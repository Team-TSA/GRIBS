"""Energy-mode integrator core (A-6 stage 2, C10; contract 164 sections 1b, 1c, 2b-4.1..2b-4.4, 3-B).

The state is [M, E_cv, x, Rt, m_out, impulse, m_gen, W_burn, H_out, Q]: the mass M and the total energy E_cv of the contents,
the burn depth x, the throat radius, the discharged and generated mass, the impulse, and three accumulators of the energy terms
(the net burn-surface term, the energy carried out by the discharge, the wall heat) that exist only so that the closure of
E_cv can be checked (energy/ledger.py).  The element inventories are B_e = M times the element fractions of the initial fill
(contract U1: the propellant and the initial fill have the same element fractions, so they stay fixed).

Each right-hand side recovers (p, T) from (E_cv, V_free, B_e) with energy/recovery.py, then
  burn rate         r = burn_rate(p)  (the legacy function of p and the geometry),
  burning mass flow mdot_burn = rho_p A_b r,   rho_p = energy_model.burn_surface.solid_density_kg_m3 (the density of the solid),
  outflow enthalpy  h = u + p V_free / M = E_cv / M + p V_free / M   (the thermodynamic density of the gas mixture is M / V_free:
                    it is taken from the STATE, never from the settings),
  nozzle flow       nozzle_state_from_chamber (the chamber state is the stagnation state, the contents at rest, H-05),
  balance           first_law_rhs with the burn-surface net term of form (a) and h_s at T_grain_inf.
The solid density rho_p and the gas density M / V_free are two different quantities (contract 164 section 1c): they are not
compared, and each has one role, in a function of its own below.

The equilibrium provider, the element fractions, the isentropic exponent gamma_s of the nozzle flow and the search range of
the state recovery are ARGUMENTS (EnergyRunInputs); the wiring from the configuration to a CEA provider is not part of this
module (run_model_energy without inputs raises NotImplementedError).  An igniter is rejected by the configuration check.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Mapping, Optional, Tuple

import numpy as np
from scipy.integrate import solve_ivp

import gribs.config
import gribs.nozzle
import gribs.propellant.burning
import gribs.simulation.integrator
from gribs.config import Config, geometry, web_thickness
from gribs.constants import R_UNIVERSAL
from gribs.energy.balance import BalanceTerms, Flow, burn_surface_term, first_law_rhs, h0_out_from_chamber, \
    h_s_for_burn_surface
from gribs.energy.config_model import BurnSurface
from gribs.energy.initial import initial_state_from_fill
from gribs.energy.ledger import closure_residuals
from gribs.energy.nozzle_from_state import _choke_threshold_pressure, nozzle_state_from_chamber
from gribs.energy.provider import EquilibriumProvider
from gribs.energy.recovery import recover_state
from gribs.energy.state import RecoveredState

IM, IE, IX, IRT, IMO, IIMP, IMG, IWB, IHO, IQ = range(10)
_ELEMENT_SUM_TOL = 1.0e-9


@dataclass(frozen=True)
class EnergyRunInputs:
    """What an energy-mode run needs besides the configuration (the wiring from the configuration is not implemented)."""
    provider: EquilibriumProvider
    element_fractions: Mapping[str, float]      # element mass fractions of the propellant and of the initial fill (U1)
    gamma_s: float                              # isentropic exponent of the nozzle flow (the source is not decided; argument)
    T_search_K: Tuple[float, float]             # searched temperature range of the state recovery
    starts_K: Tuple[float, ...]                 # start temperatures of the state recovery

    def __post_init__(self):
        if not (math.isfinite(self.gamma_s) and self.gamma_s > 1.0):
            raise ValueError(f"gamma_s must be finite and > 1, got {self.gamma_s!r}")
        fractions = dict(self.element_fractions)
        if not fractions or any(not (math.isfinite(v) and v >= 0.0) for v in fractions.values()):
            raise ValueError("element_fractions must be a non-empty table of finite, non-negative fractions")
        if abs(math.fsum(fractions.values()) - 1.0) > _ELEMENT_SUM_TOL:
            raise ValueError(f"element_fractions must add up to 1, got {math.fsum(fractions.values())!r}")
        lo, hi = self.T_search_K
        if not (math.isfinite(lo) and math.isfinite(hi) and 0.0 < lo < hi):
            raise ValueError(f"T_search_K must satisfy 0 < T_lo < T_hi, got {self.T_search_K!r}")
        if not self.starts_K or any(not (lo < t < hi) for t in self.starts_K):
            raise ValueError(f"every start temperature must lie inside T_search_K, got {list(self.starts_K)!r}")


def _burn_mass_flow(surface: BurnSurface, A_b: float, r: float) -> float:
    """mdot_burn = rho_p A_b r [kg/s]: the ONLY use of the density of the solid (a setting, contract 164 section 1c)."""
    return surface.solid_density_kg_m3 * A_b * r


def _outflow_enthalpy(E_cv: float, M: float, p: float, V_free: float) -> float:
    """h = u + p V_free / M [J/kg], with u = E_cv / M (the contents at rest): the ONLY use of the thermodynamic density of the
    gas mixture, M / V_free, and it comes from the state (contract 164 section 1c); no setting enters."""
    return E_cv / M + p * V_free / M


class _EnergyModel:
    def __init__(self, c: Config, inputs: EnergyRunInputs):
        self.c, self.inputs, self.e = c, inputs, c.energy
        self.xw = web_thickness(c)
        self.At0 = math.pi * c.R_t0 ** 2
        self.Ae = c.eps_nozzle * self.At0
        self.h_s = h_s_for_burn_surface(self.e.burn_surface, products_enthalpy_basis=getattr(
            inputs.provider, "enthalpy_convention", "undeclared"))
        self._last = None

    def B_e(self, M: float) -> dict:
        return {k: M * v for k, v in self.inputs.element_fractions.items()}

    def recover(self, M: float, E_cv: float, V_free: float) -> RecoveredState:
        key = (M, E_cv, V_free)
        if self._last is None or self._last[0] != key:
            self._last = (key, recover_state(E_cv, V_free, self.B_e(M), self.inputs.provider, self.e.condensed,
                                             self.e.recovery, T_search_K=self.inputs.T_search_K,
                                             starts_K=self.inputs.starts_K))
        return self._last[1]

    def terms(self, y, burning: bool) -> dict:
        c, e = self.c, self.e
        M, E_cv = float(y[IM]), float(y[IE])
        if not (math.isfinite(M) and M > 0.0 and math.isfinite(E_cv)):
            raise ValueError(f"the state left its range: M = {M!r} kg, E_cv = {E_cv!r} J")
        x = min(max(float(y[IX]), 0.0), self.xw)
        Rt = max(float(y[IRT]), 1.0e-9)
        At = math.pi * Rt * Rt
        Ri, Lp, Ab, V_free = geometry(x, c)
        rec = self.recover(M, E_cv, V_free)
        p, T = rec.p_Pa, rec.T_K
        if burning and Ab > 0.0:
            r = gribs.propellant.burning.burn_rate(p, c, Ab, Ri)
        else:
            r, Ab = 0.0, 0.0
        mdot_burn = _burn_mass_flow(e.burn_surface, Ab, r)
        R_eff = R_UNIVERSAL * rec.n_gas_kmol / rec.m_total_kg
        nz = nozzle_state_from_chamber(p, T, R_eff, self.inputs.gamma_s, rec.Yg, At, self.Ae, c)
        mdot_out = nz["mdot_total"]
        H0_out = h0_out_from_chamber(_outflow_enthalpy(E_cv, M, p, V_free), 0.0)
        burn_net = burn_surface_term(mdot_burn, self.h_s, p, Ab, r, form=e.burn_surface.form)
        Q = 0.0 if e.wall.adiabatic else float(e.wall.heat_rate_W)
        dE = first_law_rhs(BalanceTerms(Qdot_in=Q, burn_surface_net=burn_net, outflows=(Flow(mdot_out, H0_out),)))
        return dict(rec=rec, p=p, T=T, x=x, At=At, Ab=Ab, r=r, V_free=V_free, mdot_burn=mdot_burn, nz=nz, mdot_out=mdot_out,
                    H0_out=H0_out, burn_net=burn_net, Q=Q, dE=dE)

    def rhs(self, burning: bool):
        c = self.c

        def rhs(t, y):
            k = self.terms(y, burning)
            dRt = c.ero_throat_c * (k["p"] / c.p_ref) ** c.ero_throat_m
            return np.array([k["mdot_burn"] - k["mdot_out"], k["dE"], k["r"], dRt, k["mdot_out"], k["nz"]["F"],
                             k["mdot_burn"], k["burn_net"], k["mdot_out"] * k["H0_out"], k["Q"]])
        return rhs

    def choke_margin(self, y) -> float:
        k = self.terms(y, False)
        eps = self.Ae / k["At"]
        return k["p"] - _choke_threshold_pressure(self.c.p_a, eps, self.inputs.gamma_s)


def _event(fn, terminal=True, direction=0.0):
    fn.terminal = terminal
    fn.direction = direction
    return fn


def _counted(fun):
    state = {"calls": 0}

    def counted(t, y):
        state["calls"] += 1
        if state["calls"] > gribs.simulation.integrator.INTEGRATION_MAX_RHS_CALLS:
            raise gribs.simulation.integrator.IntegrationBudgetError(
                f"{gribs.simulation.integrator.IntegrationBudgetError.code}: more than "
                f"{gribs.simulation.integrator.INTEGRATION_MAX_RHS_CALLS} right-hand-side evaluations (t = {t:.9g} s).")
        return fun(t, y)
    return counted, state


def run_model_energy(c: Config, inputs: Optional[EnergyRunInputs] = None) -> dict:
    """Integrate the burning phase and (optionally) the blowdown phase of an energy-mode configuration.  Without inputs
    (the wiring from the configuration to the provider) it raises NotImplementedError."""
    gribs.config.validate(c)
    if getattr(c, "energy", None) is None:
        raise ValueError("run_model_energy needs a configuration with an energy_model section")
    if inputs is None:
        raise NotImplementedError(
            "energy mode: the wiring from the energy_model configuration to the equilibrium provider (propellant formulas "
            "-> element fractions -> the CEA provider) is not implemented yet (A-6 stage 2, contract 164); the integrator "
            "core takes the provider as an argument (EnergyRunInputs).")
    gribs.simulation.integrator.require_implemented_transition_policy(c)
    model = _EnergyModel(c, inputs)
    e, At0, Ae, xw = c.energy, model.At0, model.Ae, model.xw

    _Ri, _Lp, _Ab, V_free0 = geometry(0.0, c)
    start = initial_state_from_fill(e.initial_fill.temperature_K, c.p0_init, dict(inputs.element_fractions), V_free0,
                                    inputs.provider, e.condensed)
    M0, E0 = start.mass_kg, start.E_cv
    y0 = np.array([M0, E0, 0.0, c.R_t0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    init_regime = model.terms(y0, False)["nz"]["regime"]
    if c.unchoked_policy == "stop" and init_regime != "choked":
        raise ValueError("Initial state is not choked: policy 'stop' has nothing to integrate.")
    atol = np.array([1.0e-12, 1.0e-3, c.atol_x, 1.0e-14, 1.0e-12, 1.0e-9, 1.0e-12, 1.0e-3, 1.0e-3, 1.0e-3])

    ev_burn = _event(lambda t, y: y[IX] - xw, True, +1.0)
    ev_unch = _event(lambda t, y: model.choke_margin(y), c.unchoked_policy == "stop", -1.0)
    ev_rech = _event(lambda t, y: model.choke_margin(y), False, +1.0)
    ev_amb = _event(lambda t, y: model.terms(y, False)["p"] - c.p_a * 1.0005, True, -1.0)
    names = ["burnout", "unchoked", "rechoked", "ambient"]

    f1, calls1 = _counted(model.rhs(True))
    sol1 = solve_ivp(f1, (0.0, c.t_max), y0, method=c.method, events=[ev_burn, ev_unch, ev_rech, ev_amb],
                     dense_output=True, rtol=c.rtol, atol=atol, max_step=c.max_step_burn, first_step=1.0e-7)
    if not sol1.success:
        raise RuntimeError(f"Burning-phase integration failed: {sol1.message}")
    stop1 = "time-limit"
    for nm, te in zip(names, sol1.t_events):
        if nm in ("burnout", "ambient") and len(te):
            stop1 = nm
            break
        if nm == "unchoked" and c.unchoked_policy == "stop" and len(te):
            stop1 = "choking-loss"
            break

    sol2, stop2, calls2 = None, "not-run", {"calls": 0}
    if c.blowdown and stop1 == "burnout":
        y_bo = sol1.y[:, -1].copy()
        y_bo[IX] = xw
        ev_amb2 = _event(lambda t, y: model.terms(y, False)["p"] - c.p_a * 1.0005, True, -1.0)
        ev_unch2 = _event(lambda t, y: model.choke_margin(y), c.unchoked_policy == "stop", -1.0)
        f2, calls2 = _counted(model.rhs(False))
        sol2 = solve_ivp(f2, (sol1.t[-1], sol1.t[-1] + c.blowdown_tmax), y_bo, method=c.method,
                         events=[ev_amb2, ev_unch2], dense_output=True, rtol=c.rtol, atol=atol, max_step=c.max_step_blow)
        if not sol2.success:
            raise RuntimeError(f"Blowdown integration failed: {sol2.message}")
        stop2 = "ambient" if len(sol2.t_events[0]) else "time-limit"
        if c.unchoked_policy == "stop" and len(sol2.t_events[1]):
            stop2 = "choking-loss"

    return dict(sol_burn=sol1, sol_blow=sol2, stop_burn=stop1, stop_blow=stop2,
                integration_budget={"rhs_calls": calls1["calls"] + calls2["calls"]},
                Ae=Ae, At0=At0, M0=M0, E0=E0, init_regime=init_regime,
                t_events={nm: [float(v) for v in te] for nm, te in zip(names, sol1.t_events)},
                model=model)


def closure_history(run: dict) -> list:
    """The closure residuals at every node of the solution of a run (burning phase, then blowdown)."""
    model: _EnergyModel = run["model"]
    out = []
    for sol in (run["sol_burn"], run["sol_blow"]):
        if sol is None:
            continue
        for j in range(sol.t.size):
            y = sol.y[:, j]
            M = float(y[IM])
            _Ri, _Lp, _Ab, V_free = geometry(min(max(float(y[IX]), 0.0), model.xw), model.c)
            rec = model.recover(M, float(y[IE]), V_free)
            out.append(closure_residuals(
                M=M, M0=run["M0"], m_gen=float(y[IMG]), m_out=float(y[IMO]), B_e=model.B_e(M), E_cv=float(y[IE]),
                E0=run["E0"], Q_cum=float(y[IQ]), W_burn_cum=float(y[IWB]), H_out_cum=float(y[IHO]),
                V_gas=rec.V_gas_m3, V_condensed=rec.V_condensed_m3, V_free=V_free))
    return out
