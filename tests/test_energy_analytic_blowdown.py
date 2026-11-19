"""A-6 stage 2, C15: the energy-mode right-hand side against an INDEPENDENT closed-form solution (PRE_REGISTRATION_C15, F4-01, F4-02).

A closed vessel of an ideal gas discharges through a choked nozzle (no burning, no erosion, no wall heat).  The solution is
written out here from the textbook relations, not from the code under test:
    k = Cd At G' sqrt(R T0) / V,  G' = sqrt(g) (2/(g+1))^((g+1)/(2(g-1))),  tau = 2 / ((g-1) k),
    T/T0 = (1 + t/tau)^-2,  p/p0 = (T/T0)^(g/(g-1)),  M/M0 = (T/T0)^(1/(g-1)).
The provider is a synthetic ideal gas of the given gamma (assumed inputs, not data).  The tolerance of 1e-6 is 6 orders above
the rtol of 1e-12 of this integration (the global error of DOP853 is about 1e-10 or less) and 3 orders below the effect of
a wrong exponent or gamma (1e-3 or more).  No assignment and no monkeypatch.
"""
import json
import math
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import solve_ivp

import gribs.config
import gribs.config.loader
import gribs.simulation.energy_integrator as ei
from gribs.config import Config, geometry
from gribs.energy.balance import ABSOLUTE_ENTHALPY_BASIS
from gribs.energy.provider import TPState

ROOT = Path(__file__).resolve().parents[1]
G1 = json.loads((ROOT / "cases" / "golden" / "g1_reference_constant.json").read_text(encoding="utf-8"))
ENERGY = {"initial_fill": {"temperature_K": 300.0, "composition": "main_propellant_products", "element_tolerance": 1.0e-6},
          "wall": {"adiabatic": True},
          "burn_surface": {"form": "a", "grain_far_temperature_K": 293.0,
                           "solid_enthalpy": {"reference_J_kg": -2.0e6, "reference_temperature_K": 293.0, "cp_J_kgK": 1300.0}},
          "condensed": {"volume": "external_density", "phases": [], "transitions": []},
          "recovery": {"phase_fraction_tolerance": 1.0e-6, "root_tolerance": 1.0e-10, "plateau_edge_tolerance": 1.0e-6,
                       "temperature_search_K": [200.0, 8000.0], "start_temperatures_K": [600.0, 3000.0]}}
FRACTIONS = {"C": 0.3, "H": 0.7}
MW, U_F, RU = 24.0, -8.236e6, 8314.51
T0, P0 = 3000.0, 5.0e6
TOLERANCE = 1.0e-6
TIMES_OVER_TAU = (0.02, 0.05, 0.1, 0.15, 0.2)


class GammaProvider:
    """A synthetic ideal gas of the given gamma: u = U_F + cv T, M = 24 kg/kmol (assumed inputs; no condensed phase)."""
    enthalpy_convention = ABSOLUTE_ENTHALPY_BASIS

    def __init__(self, gamma):
        self.gamma, self.cv = gamma, RU / MW / (gamma - 1.0)

    def tp(self, T, p, B_e):
        return TPState(u_J_per_kg=U_F + self.cv * T, n_gas_per_kg=1.0 / MW, Yg=1.0, condensed_mass_fraction={},
                       converged=True, gamma_s=self.gamma)


def case_config():
    c = gribs.config.loader.configuration_from_document(
        {**G1, "two_phase_model": {**G1["two_phase_model"], "mode": "homogeneous_equilibrium"}, "energy_model": ENERGY}).config
    sub = type("ConfigWithEnergy", (Config,), {"energy": c.energy})
    return sub(**{**asdict(c), "t_max": 0.05, "rtol": 1.0e-8})


def discharge(gamma):
    """The numerical states and the closed-form states of the discharge at the times of TIMES_OVER_TAU * tau."""
    c = case_config()
    assert c.ero_throat_c == 0.0                                     # no erosion: the throat area is constant
    provider = GammaProvider(gamma)
    model = ei._EnergyModel(c, ei.EnergyRunInputs(provider=provider, element_fractions=FRACTIONS))
    V = geometry(model.xw, c)[3]                                     # the free volume at burnout
    Rg = RU / MW
    M0 = P0 * V / (Rg * T0)
    y0 = np.array([M0, M0 * (U_F + provider.cv * T0), model.xw, c.R_t0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    E0 = abs(y0[ei.IE])
    atol = 1.0e-12 * np.array([M0, E0, 1.0, c.R_t0, M0, 1.0, M0, E0, E0, E0])
    Gp = math.sqrt(gamma) * (2.0 / (gamma + 1.0)) ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))
    k = c.Cd * math.pi * c.R_t0 ** 2 * Gp * math.sqrt(Rg * T0) / V
    tau = 2.0 / ((gamma - 1.0) * k)
    times = [f * tau for f in TIMES_OVER_TAU]
    sol = solve_ivp(model.rhs(False), (0.0, times[-1]), y0, method="DOP853", t_eval=times, rtol=1.0e-12, atol=atol)
    assert sol.success, sol.message
    rows = []
    for t, y in zip(sol.t, sol.y.T):
        rec = model.recover(float(y[ei.IM]), float(y[ei.IE]), V)
        theta = (1.0 + t / tau) ** -2
        rows.append(dict(t_over_tau=t / tau, p=rec.p_Pa, T=rec.T_K, M=float(y[ei.IM]),
                         p_exact=P0 * theta ** (gamma / (gamma - 1.0)), T_exact=T0 * theta,
                         M_exact=M0 * theta ** (1.0 / (gamma - 1.0)), p0=P0, T0=T0, M0=M0,
                         choked=rec.p_Pa > c.p_a * ((gamma + 1.0) / 2.0) ** (gamma / (gamma - 1.0))))
    return rows


def check(gamma):
    rows = discharge(gamma)
    assert len(rows) == len(TIMES_OVER_TAU) and all(r["choked"] for r in rows)
    for r in rows:
        for name in ("p", "T", "M"):
            assert r[name] == pytest.approx(r[name + "_exact"], rel=TOLERANCE, abs=0.0), (gamma, r["t_over_tau"], name)
        isentropic = (r["T"] / r["T0"]) ** (gamma / (gamma - 1.0))
        assert r["p"] / r["p0"] == pytest.approx(isentropic, rel=TOLERANCE, abs=0.0), (gamma, r["t_over_tau"])


@pytest.mark.fast
def test_f401_the_choked_discharge_of_an_ideal_gas_with_gamma_1_2_matches_the_closed_form():
    check(1.2)


@pytest.mark.fast
def test_f402_the_choked_discharge_of_an_ideal_gas_with_gamma_1_4_matches_the_closed_form():
    check(1.4)
