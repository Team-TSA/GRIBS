"""A-6 stage 2, C15: a whole run through the melting plateau, upward in the burning and downward in the blowdown
(PRE_REGISTRATION_C15, F4-03, F4-04).

The chemistry is a SYNTHETIC provider (assumed inputs, not data): a gas (mass fraction 0.7, 24 kg/kmol, gamma 1.2, constant cv)
and one condensed species OX, solid below T_m = 2000 K and liquid above, with constant cp and L = 1.0e6 J/kg.  On the plateau
T = T_m exactly and the phase fraction phi moves between 0 and 1.  Checked: the samples on the plateau (T == T_m, 0 < phi < 1,
the direction of phi), the closure residuals at every sample, and (e) the time spent on the plateau against a value written out
here from the synthetic data (the melted mass times L over the net energy inflow), not from the recovery or from the energy
balance of the code under test.  No assignment and no monkeypatch.
"""
import json
import math
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pytest

import gribs.config
import gribs.config.loader
import gribs.simulation.energy_integrator as ei
from gribs.config import Config, geometry
from gribs.energy.balance import ABSOLUTE_ENTHALPY_BASIS
from gribs.energy.ledger import closure_residuals
from gribs.energy.provider import TPState

ROOT = Path(__file__).resolve().parents[1]
G1 = json.loads((ROOT / "cases" / "golden" / "g1_reference_constant.json").read_text(encoding="utf-8"))
RU = 8314.51
MW, YG, YC, GAMMA = 24.0, 0.7, 0.3, 1.2
T_M, CP, TREF, H_SOLID_REF, H_LIQUID_REF = 2000.0, 1300.0, 298.15, -1.6e7, -1.5e7
LATENT = H_LIQUID_REF - H_SOLID_REF                                     # L = 1.0e6 J/kg (same cp, same reference temperature)
H_SURFACE, T_FLAME = -2.0e6, 3000.0                                      # the burn surface enthalpy; the designed flame temperature
CV_G, R_G = RU / MW / (GAMMA - 1.0), RU / MW
CLOSURE_LIMIT = 1.0e-6                                                   # = energy_postprocess.CLOSURE_TOLERANCE (the run's own limit)
N_SAMPLES, MIN_PLATEAU_SAMPLES, MIN_PHI_SPAN, DWELL_RELATIVE_TOLERANCE, PHI_NOISE = 2000, 5, 0.5, 0.02, 1.0e-12


def h_solid(T):
    return H_SOLID_REF + CP * (T - TREF)


def h_liquid(T):
    return H_LIQUID_REF + CP * (T - TREF)


U_G0 = (H_SURFACE - YC * h_liquid(T_FLAME)) / YG - CV_G * (T_FLAME - 300.0) - R_G * T_FLAME   # h(T_flame) = h_surface, gas part


class PlateauProvider:
    enthalpy_convention = ABSOLUTE_ENTHALPY_BASIS

    def tp(self, T, p, B_e):
        liquid = T > T_M
        u = YG * (U_G0 + CV_G * (T - 300.0)) + YC * (h_liquid(T) if liquid else h_solid(T))
        return TPState(u_J_per_kg=u, n_gas_per_kg=YG / MW, Yg=YG, condensed_mass_fraction={"OX(L)" if liquid else "OX(s)": YC},
                       converged=True, gamma_s=GAMMA)


def phase_data(species, kind, density, enthalpy):
    return {"species": species, "phase": kind, "density_kg_m3": density, "cp_J_kgK": CP, "enthalpy_reference_J_kg": enthalpy,
            "reference_temperature_K": TREF}


ENERGY = {"initial_fill": {"temperature_K": 300.0, "composition": "main_propellant_products", "element_tolerance": 1.0e-6},
          "wall": {"adiabatic": True},
          "burn_surface": {"form": "a", "grain_far_temperature_K": 293.0,
                           "solid_enthalpy": {"reference_J_kg": H_SURFACE, "reference_temperature_K": 293.0, "cp_J_kgK": 1300.0}},
          "condensed": {"volume": "external_density",
                        "phases": [phase_data("OX(s)", "solid", 3990.0, H_SOLID_REF), phase_data("OX(L)", "liquid", 3000.0, H_LIQUID_REF)],
                        "transitions": [{"solid_species": "OX(s)", "liquid_species": "OX(L)", "melting_temperature_K": T_M}],
                        "latent_heat_check_tolerance_rel": 1.0e-3},
          "recovery": {"phase_fraction_tolerance": 1.0e-6, "root_tolerance": 1.0e-10, "plateau_edge_tolerance": 1.0e-6,
                       "temperature_search_K": [200.0, 8000.0], "start_temperatures_K": [600.0, 3000.0]}}
FRACTIONS = {"C": 0.3, "H": 0.7}


def run_case():
    c0 = gribs.config.loader.configuration_from_document(
        {**G1, "two_phase_model": {**G1["two_phase_model"], "mode": "homogeneous_equilibrium"}, "energy_model": ENERGY}).config
    sub = type("ConfigWithEnergy", (Config,), {"energy": c0.energy})
    c = sub(**{**asdict(c0), "rtol": 1.0e-8, "R_p": 0.0045, "t_max": 1.0, "blowdown": True, "blowdown_tmax": 0.3,
               "max_step_blow": 0.005})
    return c, ei.run_model_energy(c, ei.EnergyRunInputs(provider=PlateauProvider(), element_fractions=FRACTIONS))


def samples(c, run, sol, burning):
    """N_SAMPLES equally spaced samples of one phase of the run: the recovered state, the terms and the closure residuals."""
    model = run["model"]
    out = []
    for t in np.linspace(float(sol.t[0]), float(sol.t[-1]), N_SAMPLES):
        y = sol.sol(t)
        M, E = float(y[ei.IM]), float(y[ei.IE])
        V = geometry(min(max(float(y[ei.IX]), 0.0), model.xw), c)[3]
        rec = model.recover(M, E, V)
        k = model.terms(y, burning)
        res = closure_residuals(M=M, M0=run["M0"], m_gen=float(y[ei.IMG]), m_out=float(y[ei.IMO]), B_e=model.B_e(M), E_cv=E,
                                E0=run["E0"], Q_cum=float(y[ei.IQ]), W_burn_cum=float(y[ei.IWB]), H_out_cum=float(y[ei.IHO]),
                                V_gas=rec.V_gas_m3, V_condensed=rec.V_condensed_m3, V_free=V)
        out.append(dict(t=t, M=M, T=rec.T_K, p=rec.p_Pa, phi=rec.phi, on_plateau=rec.on_plateau, V=V, mdot_burn=k["mdot_burn"],
                        mdot_out=k["mdot_out"], closure=max(abs(res.R_M_normalised), abs(res.R_B_normalised),
                                                           abs(res.R_E_normalised), abs(res.R_V_normalised))))
    return out


def phi_rate(s):
    """d phi / dt on the plateau, written from the synthetic data: M Yc L phi' = mdot_burn (h_s - u_m(phi)) - mdot_out p V / M
    (T = T_m, so E = M u_m(phi); dE/dt = mdot_burn h_s - mdot_out (u + p V / M); dM/dt = mdot_burn - mdot_out)."""
    u_m = YG * (U_G0 + CV_G * (T_M - 300.0)) + YC * (h_solid(T_M) + s["phi"] * LATENT)
    return (s["mdot_burn"] * (H_SURFACE - u_m) - s["mdot_out"] * s["p"] * s["V"] / s["M"]) / (s["M"] * YC * LATENT)


def plateau_report(rows, upward):
    """The measurements of one phase: the plateau samples, T == T_m, phi 0 < phi < 1 and its direction, the closure, and the
    time on the plateau against the value from the synthetic data."""
    idx = [i for i, s in enumerate(rows) if s["on_plateau"] and s["phi"] is not None and 0.0 < s["phi"] < 1.0]
    out = dict(n_plateau=len(idx), closure_max=max(s["closure"] for s in rows))
    if not idx:
        return out
    plat = rows[idx[0]: idx[-1] + 1]
    out.update(contiguous=idx == list(range(idx[0], idx[-1] + 1)), T_all_equal=all(rows[i]["T"] == T_M for i in idx),
               phi_first=rows[idx[0]]["phi"], phi_last=rows[idx[-1]]["phi"],
               monotone=all((b["phi"] - a["phi"]) * (1.0 if upward else -1.0) >= -PHI_NOISE for a, b in zip(plat, plat[1:])))
    dt = plat[1]["t"] - plat[0]["t"] if len(plat) > 1 else math.nan
    rates = [phi_rate(s) for s in plat]
    integral = float(sum(0.5 * (a + b) * dt for a, b in zip(rates, rates[1:])))           # trapezoid rule over the plateau samples
    dwell_obs = plat[-1]["t"] - plat[0]["t"]
    dphi = out["phi_last"] - out["phi_first"]
    out.update(integral=integral, dwell_obs=dwell_obs, dwell_pred=abs(dphi) * dwell_obs / abs(integral) if integral else math.nan)
    return out


def assert_plateau(rep, upward):
    assert rep["n_plateau"] >= MIN_PLATEAU_SAMPLES, rep                                        # (a)
    assert rep["contiguous"] and rep["T_all_equal"], rep                                       # (b)
    assert rep["monotone"], rep                                                                # (c)
    assert rep["closure_max"] <= CLOSURE_LIMIT, rep                                            # (d)
    assert abs(rep["phi_last"] - rep["phi_first"]) >= MIN_PHI_SPAN and (rep["integral"] > 0.0) == upward, rep
    assert abs(rep["dwell_obs"] - rep["dwell_pred"]) <= DWELL_RELATIVE_TOLERANCE * rep["dwell_pred"], rep   # (e)


def reports():
    c, run = run_case()
    assert run["stop_burn"] == "burnout" and run["sol_blow"] is not None
    return (plateau_report(samples(c, run, run["sol_burn"], True), True),
            plateau_report(samples(c, run, run["sol_blow"], False), False))


@pytest.mark.integration
def test_f403_the_burning_phase_crosses_the_melting_plateau_upward_and_stays_on_it_as_long_as_the_data_say():
    assert_plateau(reports()[0], True)


@pytest.mark.integration
def test_f404_the_blowdown_crosses_the_melting_plateau_downward_and_stays_on_it_as_long_as_the_data_say():
    assert_plateau(reports()[1], False)
