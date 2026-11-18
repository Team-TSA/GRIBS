"""ODE system and integration (A-6 unit 5). Moved verbatim from section 7 of the entry: state indices, igniter
flow, right-hand side, events, the integration budget (INTEGRATION_MAX_*), the solver restart at the
property-table limits, run_model and require_implemented_transition_policy.

Tests patch run_model and the INTEGRATION_MAX_* limits here, e.g. "gribs.simulation.integrator.run_model";
the entry re-exports every name, so a patch on the entry would not reach this module. Names of other package
modules that tests patch are read as <module>.<name> (ADR decision 1).
"""
from __future__ import annotations

import math
import time as _integration_time
from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp

import gribs.config
import gribs.nozzle
import gribs.propellant.burning
import gribs.thermo
from gribs.config import Config, geometry, web_thickness
from gribs.constants import BACKEND_REFERENCE_CONSTANT
from gribs.nozzle import choke_margin
from gribs.thermo import theta_and_deriv


IP, IX, IRT, IMO, IIMP, IMG, IMIGN = range(7)


def igniter_mdot(t: float, c: Config) -> float:
    return c.ign_mdot if (c.ign_time > 0.0 and t <= c.ign_time) else 0.0


def igniter_injected_mass(c: Config, simulation_end_s: float) -> float:
    """Return the igniter mass admitted over the simulated time interval."""
    duration = min(max(float(simulation_end_s), 0.0), max(c.ign_time, 0.0))
    return max(c.ign_mdot, 0.0) * duration


def make_rhs(c: Config, Ae: float, burning: bool) -> Callable:
    xw = web_thickness(c)
    p_floor = 500.0

    def rhs(t, y):
        p0 = max(float(y[IP]), p_floor)
        x = min(max(float(y[IX]), 0.0), xw)
        Rt = max(float(y[IRT]), 1.0e-9)
        At = math.pi * Rt * Rt

        Ri, Lp, Ab, Vg = geometry(x, c)
        if burning and Ab > 0.0:
            r = gribs.propellant.burning.burn_rate(p0, c, Ab, Ri)
        else:
            r, Ab = 0.0, 0.0
        # TOTAL generated product mass flow (gas + condensed) [kg/s]: every
        # kilogram of burned propellant enters mt, regardless of its phase.
        m_gen_total = c.rho_p * r * Ab + igniter_mdot(t, c)

        # TOTAL nozzle outflow (gas + condensed), consistent mass basis.
        nz = gribs.nozzle.nozzle_state(p0, At, Ae, c)
        psi, dpsi = gribs.thermo.psi_and_deriv(p0, c)
        denom = Vg * (1.0 - p0 * dpsi / psi)
        if denom <= 0.0:
            raise ValueError("Non-physical factor (1 - p0*Psi'/Psi) <= 0.")

        dp = (psi * (m_gen_total - nz["mdot_total"]) - p0 * Ab * r) / denom
        dRt = c.ero_throat_c * (p0 / c.p_ref) ** c.ero_throat_m
        return np.array([
            dp,
            r,
            dRt,
            nz["mdot_total"],
            nz["F"],
            m_gen_total,
            igniter_mdot(t, c),
        ])

    return rhs


def _event(fn, terminal=True, direction=0.0):
    fn.terminal = terminal
    fn.direction = direction
    return fn



# ==============================================================================
# Integration budget and solver restart at the property-table limits
# ==============================================================================
#: Safety limits of one run_model() call (burning + blowdown phases together).
#: A run that exceeds one of them stops with IntegrationBudgetError
#: (E_INTEGRATION_BUDGET_EXCEEDED), reported by the CLI as an integration
#: failure (exit code 4).  The values are fixed (not configurable): normal runs
#: use a few thousand right-hand-side evaluations, a few seconds and at most a
#: few restarts.
INTEGRATION_MAX_RHS_CALLS = 2_000_000
INTEGRATION_MAX_WALL_S = 300.0
INTEGRATION_MAX_TABLE_RESTARTS = 10_000
_INTEGRATION_RESTART_FIRST_STEP_S = 1.0e-7


class IntegrationBudgetError(RuntimeError):
    """The ODE integration exceeded a safety limit (E_INTEGRATION_BUDGET_EXCEEDED)."""

    code = "E_INTEGRATION_BUDGET_EXCEEDED"


class _IntegrationBudget:
    def __init__(self):
        self.rhs_calls = 0
        self.t_start = _integration_time.monotonic()

    def wrap(self, fun):
        def counted(t, y):
            self.rhs_calls += 1
            if self.rhs_calls > INTEGRATION_MAX_RHS_CALLS:
                raise IntegrationBudgetError(
                    f"{IntegrationBudgetError.code}: more than {INTEGRATION_MAX_RHS_CALLS} "
                    f"right-hand-side evaluations (t = {t:.9g} s, p0 = {y[IP]:.9g} Pa).")
            if (self.rhs_calls & 1023) == 0 and _integration_time.monotonic() - self.t_start > INTEGRATION_MAX_WALL_S:
                raise IntegrationBudgetError(
                    f"{IntegrationBudgetError.code}: wall time above {INTEGRATION_MAX_WALL_S:g} s "
                    f"(t = {t:.9g} s, p0 = {y[IP]:.9g} Pa).")
            return fun(t, y)
        return counted


class _SegmentedSolution:
    """Solutions of consecutive solve_ivp segments joined into one object.

    Offers the attributes used by the post-processing (t, y, sol(), t_events,
    success, status, message, nfev); sol() dispatches each time to the segment
    that contains it (a restart time belongs to the earlier segment).
    """

    def __init__(self, segments, t_events, restarts):
        self.segments = segments
        self.t = np.concatenate([segments[0].t] + [s.t[1:] for s in segments[1:]])
        self.y = np.concatenate([segments[0].y] + [s.y[:, 1:] for s in segments[1:]], axis=1)
        self._ends = np.array([float(s.t[-1]) for s in segments])
        self.t_events = [np.asarray(te, dtype=float) for te in t_events]
        self.success = all(bool(s.success) for s in segments)
        self.status = segments[-1].status
        self.message = segments[-1].message
        self.nfev = int(sum(int(s.nfev) for s in segments))
        self.table_restarts = restarts

    def sol(self, t):
        scalar = np.ndim(t) == 0
        tt = np.atleast_1d(np.asarray(t, dtype=float))
        out = np.empty((self.y.shape[0], tt.size))
        k = np.searchsorted(self._ends[:-1], tt, side="left")
        for i, s in enumerate(self.segments):
            m = k == i
            if m.any():
                out[:, m] = s.sol(tt[m])
        return out[:, 0] if scalar else out


def _segmented_solve(fun, t_span, y0, c: Config, budget: _IntegrationBudget, events=(), **kw):
    """solve_ivp with a restart of the solver at the property-table limits.

    With property_policy = "clamp" (table backends only; reference_constant has
    no table limit) the pressure derivative of the chamber properties jumps to
    zero at p_fit_min / p_fit_max, so dp/dt is discontinuous there.  LSODA can
    then shrink its step to ~1e-10 s and, depending on its history, never
    recover (the run does not end).  Each crossing of a limit therefore ends
    the segment (terminal event) and a fresh solver continues from the crossing
    state; right after a restart the limit just crossed is watched only for the
    crossing back, so it cannot fire again at once.
    """
    events = list(events)
    nu = len(events)
    active = (c.property_policy == "clamp" and c.thermo_backend != BACKEND_REFERENCE_CONSTANT)
    t, t1 = float(t_span[0]), float(t_span[1])
    y = np.asarray(y0, dtype=float).copy()
    segments, t_events, restarts = [], [[] for _ in range(nu)], []
    last = None
    while True:
        limits = []
        if active:
            for name, pb in (("p_fit_max", float(c.p_fit_max)), ("p_fit_min", float(c.p_fit_min))):
                f = (lambda level: (lambda tt, yy: yy[IP] - level))(pb)
                f.terminal = True
                f.direction = -last[1] if (last is not None and last[0] == name) else 0.0
                limits.append((name, pb, f))
        kw_seg = dict(kw)
        if segments:
            kw_seg["first_step"] = min(_INTEGRATION_RESTART_FIRST_STEP_S, max(t1 - t, 1.0e-12))
        sol = solve_ivp(budget.wrap(fun), (t, t1), y, events=events + [lim[2] for lim in limits], **kw_seg)
        segments.append(sol)
        for i in range(nu):
            t_events[i].extend(float(v) for v in sol.t_events[i])
        if not sol.success or sol.status != 1:
            break
        if any(len(sol.t_events[i]) and getattr(events[i], "terminal", False) for i in range(nu)):
            break
        hit = None
        for j, (name, pb, _f) in enumerate(limits):
            if len(sol.t_events[nu + j]):
                p_prev = float(sol.y[IP, -2]) if sol.y.shape[1] > 1 else float(y[IP])
                hit = (name, 1.0 if p_prev < pb else -1.0)
                break
        if hit is None:
            break
        if len(restarts) >= INTEGRATION_MAX_TABLE_RESTARTS:
            raise IntegrationBudgetError(
                f"{IntegrationBudgetError.code}: more than {INTEGRATION_MAX_TABLE_RESTARTS} solver restarts "
                f"at the property-table limits (t = {float(sol.t[-1]):.9g} s).")
        t, y, last = float(sol.t[-1]), sol.y[:, -1].copy(), hit
        restarts.append({"t_s": t, "limit": hit[0], "direction": "up" if hit[1] > 0 else "down",
                         "p0_Pa": float(y[IP])})
        if t >= t1:
            break
    return _SegmentedSolution(segments, t_events, restarts)


def run_model(c: Config) -> dict:
    """Integrate the burning phase and (optionally) the blowdown phase."""
    gribs.config.validate(c)
    if getattr(c, "energy", None) is not None:
        raise NotImplementedError(
            "energy_model is valid configuration syntax but energy mode is not implemented yet (A-6 stage 2, "
            "contract 164). Remove the energy_model section to run the current model.")

    require_implemented_transition_policy(c)

    At0 = math.pi * c.R_t0 ** 2
    Ae = c.eps_nozzle * At0                      # geometric exit area, fixed
    xw = web_thickness(c)

    init_regime = gribs.nozzle.nozzle_state(c.p0_init, At0, Ae, c)["regime"]
    if c.unchoked_policy == "stop" and init_regime != "choked":
        raise ValueError("Initial state is not choked: policy 'stop' has "
                         "nothing to integrate.")

    theta0, _ = theta_and_deriv(c.p0_init, c)
    psi0, _ = gribs.thermo.psi_and_deriv(c.p0_init, c)
    y0 = np.array([
        c.p0_init,
        0.0,
        c.R_t0,
        0.0,
        0.0,
        0.0,
        0.0,
    ])
    # Initial chamber inventory as TOTAL product mass:
    #   mt0 = p0*Vg0/Psi(p0)   (two-phase mode; Psi = R_eff*T0)
    #   m0  = p0*Vg0/Theta(p0) (single-phase legacy mode)
    # ASSUMPTION: the initial chamber gas is equilibrium combustion products of
    # the main propellant at p0_init.  An initial fill of air or igniter gas of
    # a different composition is not represented exactly by this
    # single-composition model.
    m_total0 = c.p0_init * c.V_g0 / psi0

    atol = np.array([
        c.atol_p,
        c.atol_x,
        1e-14,
        1e-12,
        1e-9,
        1e-12,
        1e-12,
    ])

    ev_burn = _event(lambda t, y: y[IX] - xw, True, +1.0)
    ev_unch = _event(lambda t, y: choke_margin(max(y[IP], 1.0),
                                               math.pi * max(y[IRT], 1e-9) ** 2,
                                               Ae, c),
                     c.unchoked_policy == "stop", -1.0)
    ev_rech = _event(lambda t, y: choke_margin(max(y[IP], 1.0),
                                               math.pi * max(y[IRT], 1e-9) ** 2,
                                               Ae, c), False, +1.0)
    ev_amb = _event(lambda t, y: y[IP] - c.p_a * 1.0005, True, -1.0)
    names = ["burnout", "unchoked", "rechoked", "ambient"]

    budget = _IntegrationBudget()
    sol1 = _segmented_solve(make_rhs(c, Ae, True), (0.0, c.t_max), y0, c, budget,
                     method=c.method, events=[ev_burn, ev_unch, ev_rech, ev_amb],
                     dense_output=True, rtol=c.rtol, atol=atol,
                     max_step=c.max_step_burn, first_step=1.0e-7)
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

    sol2, stop2 = None, "not-run"
    if c.blowdown and stop1 == "burnout":
        y_bo = sol1.y[:, -1].copy()
        y_bo[IX] = xw
        ev_amb2 = _event(lambda t, y: y[IP] - c.p_a * 1.0005, True, -1.0)
        ev_unch2 = _event(lambda t, y: choke_margin(max(y[IP], 1.0),
                                                    math.pi * max(y[IRT], 1e-9) ** 2,
                                                    Ae, c),
                          c.unchoked_policy == "stop", -1.0)
        sol2 = _segmented_solve(make_rhs(c, Ae, False),
                         (sol1.t[-1], sol1.t[-1] + c.blowdown_tmax), y_bo, c, budget,
                         method=c.method, events=[ev_amb2, ev_unch2],
                         dense_output=True, rtol=c.rtol, atol=atol,
                         max_step=c.max_step_blow)
        if not sol2.success:
            raise RuntimeError(f"Blowdown integration failed: {sol2.message}")
        stop2 = "ambient" if len(sol2.t_events[0]) else "time-limit"
        if c.unchoked_policy == "stop" and len(sol2.t_events[1]):
            stop2 = "choking-loss"

    return dict(sol_burn=sol1, sol_blow=sol2, stop_burn=stop1, stop_blow=stop2,
                integration_budget={"rhs_calls": budget.rhs_calls,
                                    "table_restarts_burn": list(sol1.table_restarts),
                                    "table_restarts_blowdown": (list(sol2.table_restarts)
                                                                if sol2 is not None else [])},
                Ae=Ae, At0=At0, m_total0=m_total0,
                init_regime=init_regime,
                t_events={nm: [float(v) for v in te]
                          for nm, te in zip(names, sol1.t_events)})


def require_implemented_transition_policy(c: Config) -> None:
    """Reject valid-but-not-yet-implemented nozzle transition models."""
    if c.transition_policy == "shock":
        raise NotImplementedError(
            "nozzle.transition_policy = 'shock' is valid configuration syntax "
            "but the internal-shock nozzle transition model is not implemented. "
            "Use transition_policy = 'jump' for the current direct-switch model."
        )
