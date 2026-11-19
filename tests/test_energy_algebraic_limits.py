"""A-6 stage 2, C11: the algebraic limits of energy mode and the wiring (contract 164 section 2b-4.8; PRE_REGISTRATION_C11).

The limit: with the analytic backend (reference_constant: an ideal gas of constant gamma and molecular weight, one gas phase,
no condensed phase) the configuration is wired by energy/wiring.py to ReferenceConstantProvider, whose enthalpy at the reference
temperature is that of the solid at its far temperature (h(T_ref) = h_s).  Then a QUASI-STEADY state of the legacy model
(dp/dt = 0 at a frozen burn depth) is a quasi-steady state of energy mode with T = T_ref: for dp/dt = dT/dt = 0 with the free
volume growing by A_b r, the mass equation is the legacy one (Psi = R T_ref) and the energy equation of form (a) reduces to
h(T) = h_s.  The tests check this, the wiring and its provenance; they use the wired provider, no monkeypatch, no assignment, no
parametrize.  A dynamic orbit with a condensed volume (V_c != 0) is NOT compared with the legacy model (A-05).
"""
import json
import math
from pathlib import Path

import numpy as np
import pytest
from scipy.optimize import brentq

import gribs.config
import gribs.config.loader
import gribs.simulation.energy_integrator as ei
import gribs.simulation.integrator
import gribs.thermo
from gribs.config import geometry
from gribs.constants import R_UNIVERSAL
from gribs.energy.provider import TPState
from gribs.energy.wiring import EnergyRunInputs, build_energy_inputs
from gribs.errors import ConfigurationError

ROOT = Path(__file__).resolve().parents[1]
G1 = json.loads((ROOT / "cases" / "golden" / "g1_reference_constant.json").read_text(encoding="utf-8"))

FORMULAS = {"NH4CLO4(I)": {"N": 1, "H": 4, "Cl": 1, "O": 4}, "C4H6,butadiene": {"C": 4, "H": 6}}
ENERGY = {
    "initial_fill": {"temperature_K": 300.0, "composition": "main_propellant_products", "element_tolerance": 1.0e-6},
    "wall": {"adiabatic": True},
    "burn_surface": {"form": "a", "grain_far_temperature_K": 293.0,
                     "solid_enthalpy": {"reference_J_kg": -2.0e6, "reference_temperature_K": 293.0, "cp_J_kgK": 1300.0}},
    "condensed": {"volume": "external_density",
                  "phases": [{"species": "Al2O3(a)", "phase": "solid", "density_kg_m3": 3990.0}], "transitions": []},
    "recovery": {"phase_fraction_tolerance": 1.0e-6, "root_tolerance": 1.0e-10, "plateau_edge_tolerance": 1.0e-6,
                 "temperature_search_K": [200.0, 8000.0], "start_temperatures_K": [600.0, 3000.0]},
    "propellant_formulas": FORMULAS,
}
X_FROZEN = 2.0e-3                  # a burn depth inside the web (m): the geometry is frozen there
STEADY_RATE_TOL = 1.0e-6           # |d ln p / dt| and |d ln T / dt| of the quasi-steady state [1/s] (PRE_REGISTRATION_C11)
HALF_STEP_S = 1.0e-9


def document(energy=ENERGY):
    return {**G1, "two_phase_model": {**G1["two_phase_model"], "mode": "homogeneous_equilibrium"}, "energy_model": energy}


def energy_config(energy=ENERGY):
    return gribs.config.loader.configuration_from_document(document(energy)).config


def legacy_steady_pressure(c_legacy, x):
    """The pressure at which the legacy model has dp/dt = 0 at the frozen burn depth x (a root of its own right-hand side)."""
    Ae = c_legacy.eps_nozzle * math.pi * c_legacy.R_t0 ** 2
    rhs = gribs.simulation.integrator.make_rhs(c_legacy, Ae, True)
    return brentq(lambda p: float(rhs(0.0, np.array([p, x, c_legacy.R_t0, 0.0, 0.0, 0.0, 0.0]))[0]), 1.0e6, 3.0e7, xtol=1.0e-9,
                  rtol=1.0e-14)


def steady_state_vector(c, provider, p, T, x):
    """The energy-mode state vector at (p, T) of the ideal gas and the burn depth x."""
    V = geometry(x, c)[3]
    M = p * V / (provider.R * T)
    return np.array([M, M * provider.tp(T, p, {}).u_J_per_kg, x, c.R_t0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])


def rates_of_p_and_T(model, y, h=HALF_STEP_S):
    """d ln p / dt and d ln T / dt of the recovered state along the right-hand side (centred difference over the step h)."""
    dy = model.rhs(True)(0.0, y)
    out = []
    for sign in (1.0, -1.0):
        z = y + sign * h * dy
        _Ri, _Lp, _Ab, V = geometry(float(z[ei.IX]), model.c)
        out.append(model.recover(float(z[ei.IM]), float(z[ei.IE]), V))
    plus, minus = out
    return math.log(plus.p_Pa / minus.p_Pa) / (2.0 * h), math.log(plus.T_K / minus.T_K) / (2.0 * h)


@pytest.fixture
def legacy_c(tmp_path):
    """The legacy configuration of the same case, with its thermochemistry initialised (the state guard of conftest.py
    restores the backend after each test)."""
    c = gribs.config.loader.configuration_from_document(G1).config
    gribs.thermo.initialize_thermochemistry(c, tmp_path)
    return c


@pytest.mark.fast
class TestAlgebraicLimit:
    def test_a01_the_wired_provider_gives_the_ideal_gas_law_when_there_is_no_condensed_volume(self):
        c = energy_config()
        inputs = build_energy_inputs(c)
        assert isinstance(inputs, EnergyRunInputs) and abs(math.fsum(inputs.element_fractions.values()) - 1.0) < 1.0e-12
        prov = inputs.provider
        assert prov.R == pytest.approx(R_UNIVERSAL / c.reference_molecular_weight_kg_kmol, rel=1.0e-15, abs=0.0)
        model = ei._EnergyModel(c, inputs)
        for p, T in ((2.0e6, 2500.0), (1.5e7, 3000.0), (8.0e6, 3400.0)):
            y = steady_state_vector(c, prov, p, T, 1.0e-3)
            rec = model.recover(float(y[ei.IM]), float(y[ei.IE]), geometry(1.0e-3, c)[3])
            assert rec.V_condensed_m3 == 0.0 and rec.Yg == 1.0
            assert rec.p_Pa == pytest.approx(p, rel=1.0e-9, abs=0.0) and rec.T_K == pytest.approx(T, rel=1.0e-9, abs=0.0)
            assert rec.p_Pa * (geometry(1.0e-3, c)[3]) == pytest.approx(rec.m_total_kg * prov.R * rec.T_K, rel=1.0e-12, abs=0.0)
            assert rec.gamma_s == c.reference_gamma
        prov_info = inputs.provenance                                               # the choices of the wiring are recorded
        assert prov_info["provider"] == "ReferenceConstantProvider" and "gamma" in prov_info["gamma_s_definition"]
        assert prov_info["solid_density_source"].startswith("Config.rho_p") and "H-17" in prov_info["atomic_mass_table"]
        assert prov_info["reference_enthalpy_J_kg"] == ENERGY["burn_surface"]["solid_enthalpy"]["reference_J_kg"]
        with pytest.raises(ConfigurationError, match="propellant_formulas"):          # the formulas are required at run time
            build_energy_inputs(energy_config({k: v for k, v in ENERGY.items() if k != "propellant_formulas"}))

    def test_a02_a_quasi_steady_state_of_the_legacy_model_is_one_of_energy_mode(self, legacy_c):
        c = energy_config()
        p_star = legacy_steady_pressure(legacy_c, X_FROZEN)
        model = ei._EnergyModel(c, build_energy_inputs(c))
        y = steady_state_vector(c, model.inputs.provider, p_star, c.reference_temperature_K, X_FROZEN)
        d_ln_p, d_ln_T = rates_of_p_and_T(model, y)
        assert abs(d_ln_p) <= STEADY_RATE_TOL and abs(d_ln_T) <= STEADY_RATE_TOL
        k = model.terms(y, True)
        assert k["T"] == pytest.approx(c.reference_temperature_K, rel=1.0e-12, abs=0.0)
        # off the steady state the rates are not small: the test can fail (the legacy pressure 5 % lower)
        y_off = steady_state_vector(c, model.inputs.provider, 0.95 * p_star, c.reference_temperature_K, X_FROZEN)
        assert abs(rates_of_p_and_T(model, y_off)[0]) > 1.0e3 * STEADY_RATE_TOL

    def test_a03_the_isochoric_limit_with_a_wall_heat_gives_dT_dt_equal_q_over_m_cv(self):
        heated = {**ENERGY, "wall": {"adiabatic": False, "heat_rate_W": 500.0}}
        c = energy_config(heated)
        model = ei._EnergyModel(c, build_energy_inputs(c))
        y = steady_state_vector(c, model.inputs.provider, 4.0e4, 300.0, 0.0)           # below the ambient pressure: no flow
        k = model.terms(y, False)
        assert k["nz"]["regime"] == "no-flow" and k["mdot_out"] == 0.0
        dy = model.rhs(False)(0.0, y)                                                  # no burning: a closed chamber with a wall heat
        assert dy[ei.IM] == 0.0 and dy[ei.IE] == 500.0
        h = 1.0e-6
        plus = model.recover(float(y[ei.IM]), float(y[ei.IE]) + h * dy[ei.IE], geometry(0.0, c)[3])
        T_rate = (plus.T_K - 300.0) / h
        assert T_rate == pytest.approx(500.0 / (y[ei.IM] * model.inputs.provider.cv), rel=1.0e-6, abs=0.0)
        assert plus.p_Pa / plus.T_K == pytest.approx(4.0e4 / 300.0, rel=1.0e-9, abs=0.0)         # p/T is constant at fixed V and M

    def test_a04_the_first_law_identity_holds_for_the_right_hand_side(self):
        c = energy_config()
        model = ei._EnergyModel(c, build_energy_inputs(c))
        for p, T, x in ((6.0e6, 2800.0, 1.0e-3), (1.2e7, 3000.0, 3.0e-3), (2.5e6, 2400.0, 5.0e-3)):
            y = steady_state_vector(c, model.inputs.provider, p, T, x)
            k, dy = model.terms(y, True), model.rhs(True)(0.0, y)
            V = geometry(x, c)[3]
            h_chamber = y[ei.IE] / y[ei.IM] + k["p"] * V / y[ei.IM]                         # u + p V / M, written out here
            h_s = ENERGY["burn_surface"]["solid_enthalpy"]["reference_J_kg"]               # at T_grain_inf = its reference temperature
            expected = k["mdot_burn"] * h_s - k["p"] * k["Ab"] * k["r"] - k["mdot_out"] * h_chamber
            assert dy[ei.IE] == pytest.approx(expected, rel=1.0e-12, abs=0.0)
            assert dy[ei.IM] == pytest.approx(k["mdot_burn"] - k["mdot_out"], rel=1.0e-14, abs=0.0)

    def test_a05_a_condensed_volume_is_not_compared_with_the_legacy_model(self, legacy_c):
        class WithCondensed:
            enthalpy_convention = build_energy_inputs(energy_config()).provider.enthalpy_convention
            base = build_energy_inputs(energy_config()).provider

            def tp(self, T, p, B_e):
                s = self.base.tp(T, p, B_e)
                return TPState(u_J_per_kg=s.u_J_per_kg, n_gas_per_kg=s.n_gas_per_kg, Yg=0.7,
                               condensed_mass_fraction={"Al2O3(a)": 0.3}, gamma_s=s.gamma_s)
        c = energy_config()
        base_inputs = build_energy_inputs(c)
        model = ei._EnergyModel(c, EnergyRunInputs(provider=WithCondensed(), element_fractions=base_inputs.element_fractions))
        p_star = legacy_steady_pressure(legacy_c, X_FROZEN)
        y = steady_state_vector(c, base_inputs.provider, p_star, c.reference_temperature_K, X_FROZEN)
        d_ln_p, d_ln_T = rates_of_p_and_T(model, y)
        assert max(abs(d_ln_p), abs(d_ln_T)) > 1.0e3 * STEADY_RATE_TOL           # the legacy steady state is not one with V_c != 0
        for module in (ei, gribs.energy.wiring, gribs.energy.ledger):
            assert not [n for n in dir(module) if "legacy" in n.lower() and "compare" in n.lower()]   # no function compares them
        assert gribs.config.Config().two_phase_condensed_volume == "neglected"          # legacy neglects the condensed volume
