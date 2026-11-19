"""A-6 stage 2, C10: the energy-mode integrator core, the closure residuals, the igniter rejection and the two densities
(contract 164 sections 1b, 1c, 2b-4.1..2b-4.4; PRE_REGISTRATION_C10 and REV1).

The chemistry is a SYNTHETIC ideal-gas provider (constant cv, one gas, no condensed phase); the values are assumed inputs for the
software contract, not data of any propellant.  No assignment to attributes, no monkeypatch, no parametrize (the patch-site
census of test_a6_patch_targets.py stays as it is): a configuration with changed fields is a Config subclass that carries the
energy settings as a class attribute, and the provider is an argument.
"""
import ast
import functools
import json
import math
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pytest

import gribs.config
import gribs.config.loader
import gribs.simulation.energy_integrator as ei
import gribs.simulation.integrator
from gribs.config import Config, geometry
from gribs.energy.balance import ABSOLUTE_ENTHALPY_BASIS
from gribs.energy.config_model import EnergyModeSettings
from gribs.energy.energy_validation import IGNITER_REJECTED_CODE
from gribs.energy.exceptions import StateRecoveryError
from gribs.energy.ledger import ClosureError, ClosureResiduals, check_closure, closure_residuals
from gribs.energy.provider import TPState
from gribs.errors import ConfigurationError

ROOT = Path(__file__).resolve().parents[1]
G1 = json.loads((ROOT / "cases" / "golden" / "g1_reference_constant.json").read_text(encoding="utf-8"))

SOLID_H = {"reference_J_kg": -2.0e6, "reference_temperature_K": 293.0, "cp_J_kgK": 1300.0}
SURFACE = {"form": "a", "grain_far_temperature_K": 293.0, "solid_enthalpy": SOLID_H}
FILL = {"temperature_K": 300.0, "composition": "main_propellant_products", "element_tolerance": 1.0e-6}
CONDENSED = {"volume": "external_density", "phases": [], "transitions": []}
RECOVERY = {"phase_fraction_tolerance": 1.0e-6, "root_tolerance": 1.0e-10, "plateau_edge_tolerance": 1.0e-6, "temperature_search_K": [200.0, 8000.0], "start_temperatures_K": [600.0, 3000.0]}
ENERGY = {"initial_fill": FILL, "wall": {"adiabatic": True}, "burn_surface": SURFACE, "condensed": CONDENSED,
          "recovery": RECOVERY}
FRACTIONS = {"C": 0.3, "H": 0.7}
CLOSURE_TOL = 1.0e-9


class IdealProvider:
    """A synthetic gas: u = U_F + cv T, M = 24 kg/kmol, gamma = 1.2 (assumed inputs; no condensed phase)."""
    enthalpy_convention = ABSOLUTE_ENTHALPY_BASIS
    MW = 24.0
    GAMMA = 1.2
    CV = 8314.51 / 24.0 / 0.2
    U_F = -8.236e6

    def tp(self, T, p, B_e):
        return TPState(u_J_per_kg=self.U_F + self.CV * T, n_gas_per_kg=1.0 / self.MW, Yg=1.0, condensed_mass_fraction={},
                       converged=True, gamma_s=self.GAMMA)


class FailingProvider(IdealProvider):
    """The equilibrium stops converging after a number of calls (a counter in the object, not a patch)."""
    def __init__(self, calls_before_failure):
        self.budget = [calls_before_failure]

    def tp(self, T, p, B_e):
        self.budget[0] -= 1
        if self.budget[0] < 0:
            return TPState(u_J_per_kg=0.0, n_gas_per_kg=1.0 / self.MW, Yg=1.0, converged=False)
        return super().tp(T, p, B_e)


def document(energy=ENERGY, **sections):
    return {**G1, "two_phase_model": {**G1["two_phase_model"], "mode": "homogeneous_equilibrium"}, "energy_model": energy,
            **sections}


def config(energy=ENERGY, sections=None, **changes):
    c = gribs.config.loader.configuration_from_document(document(energy, **(sections or {}))).config
    sub = type("ConfigWithEnergy", (Config,), {"energy": c.energy})
    return sub(**{**asdict(c), "t_max": 0.05, "rtol": 1.0e-8, **changes})


def inputs(provider=None, **changes):
    values = dict(provider=provider or IdealProvider(), element_fractions=FRACTIONS)
    return ei.EnergyRunInputs(**{**values, **changes})


@functools.lru_cache(maxsize=None)
def base_run():
    return ei.run_model_energy(config(), inputs())


@functools.lru_cache(maxsize=None)
def base_history():
    return tuple(ei.closure_history(base_run()))


def worst(name):
    return max(abs(getattr(r, name)) for r in base_history())


def state_vector(c, M, T):
    """A state at x = 0 of mass M and temperature T of the synthetic gas (E_cv = M u)."""
    return np.array([M, M * (IdealProvider.U_F + IdealProvider.CV * T), 0.0, c.R_t0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])


@pytest.mark.fast
def test_e01_the_mass_closure_residual_is_small_along_the_run():
    assert len(base_history()) > 50 and base_run()["stop_burn"] == "time-limit"
    assert worst("R_M_normalised") <= CLOSURE_TOL
    assert worst("R_M") <= CLOSURE_TOL * max(base_run()["sol_burn"].y[ei.IM])            # the absolute value, kept apart


@pytest.mark.fast
def test_e02_the_element_closure_residual_is_small_along_the_run():
    assert worst("R_B_normalised") <= 1.0e-12 and worst("R_B") <= 1.0e-12 * max(base_run()["sol_burn"].y[ei.IM])


@pytest.mark.fast
def test_e03_the_energy_closure_residual_is_small_along_the_run():
    assert worst("R_E_normalised") <= CLOSURE_TOL
    assert worst("R_E") <= CLOSURE_TOL * max(abs(base_run()["sol_burn"].y[ei.IE]))
    assert base_run()["sol_burn"].y[ei.IHO, -1] != 0.0                 # the discharge carried energy out (H0 < 0: formation enthalpy)


@pytest.mark.fast
def test_e04_the_volume_closure_residual_is_small_along_the_run():
    assert worst("R_V_normalised") <= CLOSURE_TOL and worst("R_V") <= CLOSURE_TOL * geometry(0.0, config())[3] * 100.0


@pytest.mark.fast
def test_e05_the_burning_mass_flow_is_rho_p_times_the_swept_volume():
    c, run = config(), base_run()
    y = run["sol_burn"].y[:, -1]
    swept = geometry(float(y[ei.IX]), c)[3] - geometry(0.0, c)[3]                          # the free volume gained by the burn
    assert y[ei.IMG] == pytest.approx(c.rho_p * swept, rel=1.0e-6, abs=0.0)


@pytest.mark.fast
def test_e06_the_outflow_enthalpy_comes_from_the_state_not_from_the_settings():
    c1 = config()
    c2 = config(rho_p=3.0 * c1.rho_p)                                      # the one density of the solid: Config.rho_p (B-15)
    y = state_vector(c1, 4.0e-4, 1500.0)
    k1, k2 = ei._EnergyModel(c1, inputs()).terms(y, True), ei._EnergyModel(c2, inputs()).terms(y, True)
    V_free = geometry(0.0, c1)[3]
    h_state = y[ei.IE] / y[ei.IM] + k1["p"] * V_free / y[ei.IM]                          # u + p V_free / M, written out here
    assert k1["H0_out"] == pytest.approx(h_state, rel=1.0e-13, abs=0.0) and k2["H0_out"] == k1["H0_out"]
    assert k2["mdot_burn"] == pytest.approx(3.0 * k1["mdot_burn"], rel=1.0e-13, abs=0.0)      # the solid density scales the burn only
    assert k1["mdot_burn"] == pytest.approx(c1.rho_p * k1["Ab"] * k1["r"], rel=1.0e-13, abs=0.0)


def _attributes_in(function_name):
    tree = ast.parse(Path(ei.__file__).read_text(encoding="utf-8"))
    found = [n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == function_name]
    assert len(found) == 1
    return {a.attr for a in ast.walk(found[0]) if isinstance(a, ast.Attribute)}, [a.arg for a in found[0].args.args]


@pytest.mark.fast
def test_e07_the_two_densities_are_not_mixed_in_the_settings_or_in_the_code():
    keys = set()
    for section in (SURFACE, CONDENSED, FILL, RECOVERY, SOLID_H):
        keys |= set(section)
    keys |= {k for ph in CONDENSED["phases"] for k in ph}
    assert not {k for k in keys if "density" in k or k.startswith("rho")}       # B-15: no density key for the solid ...
    assert not any("gas" in k or "rho_th" in k or "thermodynamic" in k for k in keys)   # ... and none for the gas mixture
    assert isinstance(config().energy, EnergyModeSettings) and config().rho_p > 0.0     # the one rho_p of the solid is Config.rho_p
    burn_attrs, burn_args = _attributes_in("_burn_mass_flow")
    enthalpy_attrs, enthalpy_args = _attributes_in("_outflow_enthalpy")
    assert "rho_p" in burn_attrs and burn_args == ["c", "A_b", "r"]               # the solid density: in the burn flow ...
    assert "rho_p" not in enthalpy_attrs and enthalpy_args == ["E_cv", "M", "p", "V_free"]   # ... and not in h
    source = Path(ei.__file__).read_text(encoding="utf-8")
    uses = [n for n in ast.walk(ast.parse(source)) if isinstance(n, ast.Attribute) and n.attr == "rho_p"]
    assert "solid_density_kg_m3" not in source and len(uses) == 1                # one use of rho_p in the whole module
    assert "rho_s" not in source.replace("rho_p", "")                            # the symbol of the old mixed use is not used


@pytest.mark.fast
def test_e08_an_igniter_is_rejected_in_energy_mode_with_an_explicit_code():
    with_igniter = {"igniter": {"mass_flow_kg_s": 0.001, "duration_s": 0.1}}
    c = config(sections=with_igniter)
    for call in (lambda: gribs.config.validate(c), lambda: gribs.simulation.integrator.run_model(c)):
        with pytest.raises(ConfigurationError, match=IGNITER_REJECTED_CODE):
            call()
    assert IGNITER_REJECTED_CODE == "E_ENERGY_MODE_IGNITER_REJECTED"
    for harmless in ({"mass_flow_kg_s": 0.001, "duration_s": 0.0}, {"mass_flow_kg_s": 0.0, "duration_s": 0.1}):
        gribs.config.validate(config(sections={"igniter": harmless}))              # no igniter flow: accepted
    legacy = gribs.config.loader.configuration_from_document({**G1, **with_igniter}).config
    gribs.config.validate(legacy)                                                  # a legacy configuration keeps its igniter


@pytest.mark.fast
def test_e09_the_wall_heat_is_accumulated_and_the_energy_still_closes():
    heated = {**ENERGY, "wall": {"adiabatic": False, "heat_rate_W": 1000.0}}
    run = ei.run_model_energy(config(heated), inputs())
    t_end, q_end = run["sol_burn"].t[-1], run["sol_burn"].y[ei.IQ, -1]
    assert q_end == pytest.approx(1000.0 * t_end, rel=1.0e-9, abs=0.0)
    assert max(abs(r.R_E_normalised) for r in ei.closure_history(run)) <= CLOSURE_TOL
    assert base_run()["sol_burn"].y[ei.IQ, -1] == 0.0                              # an adiabatic wall adds nothing


@pytest.mark.fast
def test_e10_a_closed_adiabatic_chamber_does_not_change():
    c = config()
    y = state_vector(c, 1.0e-6, 300.0)                                             # p about 1e4 Pa, below the ambient pressure
    k = ei._EnergyModel(c, inputs()).terms(y, False)
    assert k["p"] < c.p_a and k["nz"]["regime"] == "no-flow"
    dy = ei._EnergyModel(c, inputs()).rhs(False)(0.0, y)
    assert not np.any(dy)                                                          # M, E_cv, m_out, H_out, Q ... all stay


@pytest.mark.fast
def test_e11_a_choked_discharge_removes_mass_and_the_enthalpy_of_the_state():
    c = config()
    model = ei._EnergyModel(c, inputs())
    y = state_vector(c, 4.0e-4, 2500.0)
    k, dy = model.terms(y, False), model.rhs(False)(0.0, y)
    assert k["nz"]["regime"] == "choked" and k["mdot_out"] > 0.0
    assert dy[ei.IM] == -k["mdot_out"] and dy[ei.IMO] == k["mdot_out"]
    h_state = y[ei.IE] / y[ei.IM] + k["p"] * geometry(0.0, c)[3] / y[ei.IM]
    assert dy[ei.IE] == pytest.approx(-k["mdot_out"] * h_state, rel=1.0e-13, abs=0.0)       # no burning, adiabatic, fixed walls
    assert dy[ei.IHO] == pytest.approx(k["mdot_out"] * h_state, rel=1.0e-13, abs=0.0)


@pytest.mark.fast
def test_e12_the_discharge_coefficient_scales_the_mass_flow_only():
    y = state_vector(config(), 4.0e-4, 2500.0)
    full = ei._EnergyModel(config(), inputs()).rhs(False)(0.0, y)
    half = ei._EnergyModel(config(Cd=0.5), inputs()).rhs(False)(0.0, y)
    assert half[ei.IMO] == pytest.approx(0.5 * full[ei.IMO], rel=1.0e-12, abs=0.0)
    assert half[ei.IHO] / half[ei.IMO] == pytest.approx(full[ei.IHO] / full[ei.IMO], rel=1.0e-13, abs=0.0)   # H0 per kg unchanged


@pytest.mark.fast
def test_e13_a_failing_equilibrium_stops_the_run_it_is_not_swallowed():
    with pytest.raises(StateRecoveryError):
        ei.run_model_energy(config(), inputs(FailingProvider(calls_before_failure=400)))


@pytest.mark.fast
def test_e14_burnout_and_the_blowdown_end_the_phases():
    short_web = dict(R_p=0.0045, t_max=1.0, blowdown=True, blowdown_tmax=0.3, max_step_blow=0.005)
    c = config(**short_web)
    run = ei.run_model_energy(c, inputs())
    xw = gribs.config.web_thickness(c)
    assert run["stop_burn"] == "burnout" and run["t_events"]["burnout"] and run["sol_burn"].y[ei.IX, -1] == pytest.approx(xw, abs=1e-9)
    assert run["sol_blow"] is not None and run["stop_blow"] in ("ambient", "time-limit", "choking-loss")
    assert run["sol_blow"].y[ei.IM, -1] < run["sol_blow"].y[ei.IM, 0]               # the blowdown empties the chamber
    assert max(abs(r.R_M_normalised) for r in ei.closure_history(run)) <= CLOSURE_TOL


@pytest.mark.fast
def test_e15_a_closure_residual_above_the_tolerance_raises():
    good = dict(M=2.0, M0=1.0, m_gen=1.5, m_out=0.5, B_e={"C": 0.6, "H": 1.4}, E_cv=-5.0, E0=-4.0, Q_cum=0.5, W_burn_cum=2.0,
                H_out_cum=3.5, V_gas=0.3, V_condensed=0.2, V_free=0.5)
    check_closure(closure_residuals(**good), tolerance=1.0e-12)
    for broken in (dict(m_out=0.6), dict(B_e={"C": 0.6, "H": 1.5}), dict(E_cv=-5.1), dict(V_gas=0.31)):
        with pytest.raises(ClosureError):
            check_closure(closure_residuals(**{**good, **broken}), tolerance=1.0e-3)
    for bad in (0.0, -1.0, math.nan):
        with pytest.raises(ClosureError):
            check_closure(closure_residuals(**good), tolerance=bad)
    with pytest.raises(ClosureError):
        closure_residuals(**{**good, "M": math.nan})
    assert isinstance(closure_residuals(**good), ClosureResiduals)


@pytest.mark.fast
def test_e16_a_run_needs_valid_inputs_a_gamma_above_one_from_the_state_and_a_declared_basis():
    for bad in (dict(element_fractions={}), dict(element_fractions={"C": 0.5}), dict(element_fractions={"C": 1.5, "H": -0.5})):
        with pytest.raises(ValueError):
            inputs(**bad)

    class NoGamma(IdealProvider):
        GAMMA = None

    class LowGamma(IdealProvider):
        GAMMA = 1.0

    class NanGamma(IdealProvider):
        GAMMA = math.nan
    for provider in (NoGamma(), LowGamma(), NanGamma()):                           # B-14: gamma_s > 1 or the run stops
        with pytest.raises(StateRecoveryError, match="gamma_s"):
            ei.run_model_energy(config(), inputs(provider))

    class NoBasis:
        def tp(self, T, p, B_e):
            return IdealProvider().tp(T, p, B_e)
    with pytest.raises(ValueError, match="absolute"):
        ei.run_model_energy(config(), inputs(NoBasis()))                            # h_s needs the basis of the products
    legacy = gribs.config.loader.configuration_from_document(G1).config
    with pytest.raises(ValueError, match="energy_model"):
        ei.run_model_energy(legacy, inputs())
    assert base_run()["provenance"]["gamma_s_min"] == base_run()["provenance"]["gamma_s_max"] == 1.2   # recorded in the provenance
