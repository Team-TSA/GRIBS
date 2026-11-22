"""N1d (report 246; approval records 79-83): energy mode takes the volume of the condensed phases as 0 (V0, the convention of CEA,
RP-1311 p.4).  Tests N-1 .. N-9 of report 244 section 6 / report 246 section 6, on SYNTHETIC providers (assumed inputs, not data).

N-1 the pressure is the gas pressure and the gas fills V_free (a condensed species with a nominal volume does not count)
N-2 a run without a condensed phase: alpha = 0 and no condensed volume
N-3 alpha = (sum Y/rho)/(V_free/M + sum Y/rho) at a node (hand formula), between the two phase values on the plateau
N-4 rho_c enters alpha only: doubling rho_c changes no other column and no other value of the summary (bit-identical)
N-5 a condensed species without a declared density stops the run (the behaviour kept by decision J-3)
N-6 the summary entry of alpha agrees with the column; no warning or pass/fail depends on alpha
N-7 E_cv/M + p V_free/M = u + R_eff T (within the volume tolerance), the V0 form of the outflow enthalpy
N-9a plateau: p does not depend on phi and the recovered phi is the original one (a provider whose two ends are the same: n- = n+,
     no sensible heat, no pressure term in u: the conditions under which the linear interpolation is exact)
N-9b plateau: with a composition that differs between the two ends the error of the linear interpolation is first order in d
     (the ratio of the errors at d and d/2 is in [1.8, 2.2]: decision of the 83rd record, item 3 (b))
The auxiliary residual r_h = (h - u - p V_free/M)/h is in test_energy_cea_tp.py (it needs the CEA provider).
"""
import dataclasses
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
import gribs.simulation.energy_postprocess as ep
from gribs.config import Config
from gribs.constants import R_UNIVERSAL
from gribs.energy.balance import ABSOLUTE_ENTHALPY_BASIS
from gribs.energy.config_model import CondensedData, CondensedPhaseSpec, PhaseTransitionSpec, RecoverySettings
from gribs.energy.exceptions import StateRecoveryError
from gribs.energy.initial import initial_state_from_fill, InitialFillError
from gribs.energy.provider import TPState
from gribs.energy.recovery import _Solver, plateau_phase_fraction, recover_state

ROOT = Path(__file__).resolve().parents[1]
G1 = json.loads((ROOT / "cases" / "golden" / "g1_reference_constant.json").read_text(encoding="utf-8"))
RU = R_UNIVERSAL
T_M, CP, TREF, H_S_REF, H_L_REF = 2000.0, 1300.0, 298.15, -1.6e7, -1.5e7
RHO_S, RHO_L = 3990.0, 3000.0
YG, YC, MG, CV_G, U_G0 = 0.7, 0.3, 24.0, 900.0, -3.0e6
VOLUME_TOL = 1.0e-8                                                         # recovery._VOLUME_RESIDUAL_TOL


def h_s(T):
    return H_S_REF + CP * (T - TREF)


def h_l(T):
    return H_L_REF + CP * (T - TREF)


SOLID = CondensedPhaseSpec("OX(s)", "solid", RHO_S, CP, H_S_REF, TREF)
LIQUID = CondensedPhaseSpec("OX(L)", "liquid", RHO_L, CP, H_L_REF, TREF)
TRANSITION = PhaseTransitionSpec("OX(s)", "OX(L)", T_M)
CONDENSED = CondensedData("neglected", (SOLID, LIQUID), (TRANSITION,), 1.0e-3)
RECOVERY = RecoverySettings(phase_fraction_tolerance=1.0e-6, root_tolerance=1.0e-10, plateau_edge_tolerance=1.0e-6)
B_E = {"X": 1.0}
SEARCH, STARTS = (300.0, 6000.0), (600.0, 1500.0, 3500.0)
V_FREE = 1.0 * (YG / MG) * RU * T_M / 3.0e6                                  # 3 MPa at T_m, mass 1 kg


@dataclasses.dataclass(frozen=True)
class Melting:
    """Ideal gas (n = YG/MG, constant cv) plus a condensed species OX, solid up to T_m and liquid above."""
    def tp(self, T, p, B_e):
        liquid = T > T_M
        u = YG * (U_G0 + CV_G * (T - 300.0)) + YC * (h_l(T) if liquid else h_s(T))
        return TPState(u_J_per_kg=u, n_gas_per_kg=YG / MG, Yg=YG, condensed_mass_fraction={"OX(L)" if liquid else "OX(s)": YC})


@dataclasses.dataclass(frozen=True)
class IdealGas:
    def tp(self, T, p, B_e):
        return TPState(u_J_per_kg=U_G0 + CV_G * (T - 300.0), n_gas_per_kg=1.0 / MG, Yg=1.0)


@pytest.mark.fast
def test_v01_n1_the_pressure_is_the_gas_pressure_and_the_gas_fills_v_free():
    T = 3000.0
    E = Melting().tp(T, 1.0e6, B_E).u_J_per_kg                              # M = 1 kg
    rec = recover_state(E, V_FREE, B_E, Melting(), CONDENSED, RECOVERY, T_search_K=SEARCH, starts_K=STARTS)
    assert rec.T_K == pytest.approx(T, rel=1.0e-9) and not rec.on_plateau
    assert rec.V_gas_m3 == pytest.approx(V_FREE, rel=VOLUME_TOL)                     # V0: the gas fills V_free
    assert rec.p_Pa == pytest.approx(1.0 * (YG / MG) * RU * rec.T_K / V_FREE, rel=1.0e-12)   # p V_gas = M n Ru T
    nominal = YC / RHO_L                                                    # M sum Y/rho of the liquid
    assert rec.V_condensed_diagnostic_m3 == pytest.approx(nominal, rel=1.0e-12) and rec.V_condensed_diagnostic_m3 > 0.0
    assert rec.V_gas_m3 + rec.V_condensed_diagnostic_m3 != pytest.approx(V_FREE, rel=1.0e-6)    # the nominal volume is not counted


@pytest.mark.fast
def test_v02_n2_without_a_condensed_phase_the_ideal_gas_law_holds_and_there_is_no_condensed_volume():
    none = CondensedData("neglected", (), ())
    T = 2500.0
    E = IdealGas().tp(T, 1.0e6, B_E).u_J_per_kg
    rec = recover_state(E, V_FREE, B_E, IdealGas(), none, RECOVERY, T_search_K=SEARCH, starts_K=STARTS)
    assert rec.T_K == pytest.approx(T, rel=1.0e-9)
    assert rec.p_Pa == pytest.approx((1.0 / MG) * RU * T / V_FREE, rel=1.0e-9)
    assert rec.V_condensed_diagnostic_m3 == 0.0 and rec.V_gas_m3 == pytest.approx(V_FREE, rel=VOLUME_TOL)


@pytest.mark.fast
def test_v05_n5_a_condensed_species_without_a_declared_density_stops():
    class Undeclared:
        def tp(self, T, p, B_e):
            return TPState(u_J_per_kg=-1.0e6, n_gas_per_kg=YG / MG, Yg=YG, condensed_mass_fraction={"ZZ(s)": YC})
    with pytest.raises(InitialFillError, match="no density for the condensed species"):
        initial_state_from_fill(300.0, 2.0e5, {"C": 0.5, "H": 0.5}, 0.02, Undeclared(), CONDENSED)
    E = Undeclared().tp(1000.0, 1.0e6, B_E).u_J_per_kg
    with pytest.raises(StateRecoveryError, match="no density for the condensed species"):
        recover_state(E, V_FREE, B_E, Undeclared(), CONDENSED, RECOVERY, T_search_K=SEARCH, starts_K=STARTS)


# ---- N-3, N-4, N-6: a whole run with a condensed phase ------------------------------------------------------------------------
class PlateauProvider:
    enthalpy_convention = ABSOLUTE_ENTHALPY_BASIS
    R_G = RU / 24.0
    CV = RU / 24.0 / 0.2
    GAMMA = 1.2
    T_FLAME, H_SURFACE = 3000.0, -2.0e6
    U0 = (H_SURFACE - YC * h_l(T_FLAME)) / YG - CV * (T_FLAME - 300.0) - R_G * T_FLAME

    S_REF, P_REF = 6000.0, 1.0e6

    @classmethod
    def entropy(cls, T, p):
        """The specific entropy per kg of TOTAL products (added for the numerical nozzle flow, decisions 116 and 117): the ideal gas
        (constant cp, R ln p) plus the incompressible condensed phase (constant cp) with the jump L / T_m at T_m."""
        s_gas = YG * ((cls.CV + cls.R_G) * math.log(T / 300.0) - cls.R_G * math.log(p / cls.P_REF))
        s_cond = YC * (CP * math.log(T / TREF) + ((H_L_REF - H_S_REF) / T_M if T > T_M else 0.0))
        return cls.S_REF + s_gas + s_cond

    def tp(self, T, p, B_e):
        liquid = T > T_M
        u = YG * (self.U0 + self.CV * (T - 300.0)) + YC * (h_l(T) if liquid else h_s(T))
        return TPState(u_J_per_kg=u, n_gas_per_kg=YG / 24.0, Yg=YG, condensed_mass_fraction={"OX(L)" if liquid else "OX(s)": YC},
                       converged=True, gamma_s=self.GAMMA, entropy_J_per_kgK=self.entropy(T, p))


def _phase(species, kind, density, enthalpy):
    return {"species": species, "phase": kind, "density_kg_m3": density, "cp_J_kgK": CP, "enthalpy_reference_J_kg": enthalpy,
            "reference_temperature_K": TREF}


def _energy_doc(scale):
    return {"initial_fill": {"temperature_K": 300.0, "composition": "main_propellant_products", "element_tolerance": 1.0e-6},
            "wall": {"adiabatic": True},
            "burn_surface": {"form": "a", "grain_far_temperature_K": 293.0,
                             "solid_enthalpy": {"reference_J_kg": -2.0e6, "reference_temperature_K": 293.0, "cp_J_kgK": 1300.0}},
            "condensed": {"volume": "neglected",
                          "phases": [_phase("OX(s)", "solid", RHO_S * scale, H_S_REF), _phase("OX(L)", "liquid", RHO_L * scale, H_L_REF)],
                          "transitions": [{"solid_species": "OX(s)", "liquid_species": "OX(L)", "melting_temperature_K": T_M}],
                          "latent_heat_check_tolerance_rel": 1.0e-3},
            "recovery": {"phase_fraction_tolerance": 1.0e-6, "root_tolerance": 1.0e-10, "plateau_edge_tolerance": 1.0e-6,
                         "temperature_search_K": [200.0, 8000.0], "start_temperatures_K": [600.0, 3000.0]}}


@functools.lru_cache(maxsize=None)
def run_with_density_scale(scale):
    c0 = gribs.config.loader.configuration_from_document(
        {**G1, "two_phase_model": {**G1["two_phase_model"], "mode": "homogeneous_equilibrium"},
         "energy_model": _energy_doc(scale)}).config
    sub = type("ConfigWithEnergy", (Config,), {"energy": c0.energy})
    c = sub(**{**asdict(c0), "rtol": 1.0e-8, "R_p": 0.0045, "t_max": 1.0, "blowdown": True, "blowdown_tmax": 0.3,
               "max_step_blow": 0.005})
    run = ei.run_model_energy(c, ei.EnergyRunInputs(provider=PlateauProvider(), element_fractions={"C": 0.3, "H": 0.7}))
    history = ep._history_and_closure(run)
    return c, run, history


@pytest.mark.integration
def test_v03_n3_alpha_is_the_hand_formula_at_each_node():
    _c, _run, h = run_with_density_scale(1.0)
    alpha, T, M, V = h["condensed_volume_fraction"], h["temperature_K"], h["mass_kg"], h["free_volume_m3"]
    assert alpha.shape == T.shape and np.all(alpha > 0.0)
    for i in range(T.size):
        a_s = (M[i] * YC / RHO_S) / (V[i] + M[i] * YC / RHO_S)
        a_l = (M[i] * YC / RHO_L) / (V[i] + M[i] * YC / RHO_L)
        if T[i] < T_M:
            assert alpha[i] == pytest.approx(a_s, rel=1.0e-9)
        elif T[i] > T_M:
            assert alpha[i] == pytest.approx(a_l, rel=1.0e-9)
        else:
            assert min(a_s, a_l) * (1.0 - 1.0e-9) <= alpha[i] <= max(a_s, a_l) * (1.0 + 1.0e-9)
    assert (T < T_M).any() and (T > T_M).any()                               # the run passes through both phases


@pytest.mark.integration
def test_v04_n4_rho_c_enters_alpha_only():
    _c, run1, h1 = run_with_density_scale(1.0)
    _c, run2, h2 = run_with_density_scale(2.0)
    assert list(h1) == list(h2) == list(ep.HISTORY_COLUMNS)
    for name in ep.HISTORY_COLUMNS:
        if name != "condensed_volume_fraction":
            assert np.array_equal(h1[name], h2[name]), name                  # bit-identical
    assert not np.array_equal(h1["condensed_volume_fraction"], h2["condensed_volume_fraction"])
    assert np.allclose(h2["condensed_volume_fraction"] / h1["condensed_volume_fraction"], 0.5, rtol=1.0e-2)   # alpha ~ 1/rho (alpha << 1)
    assert run1["M0"] == run2["M0"] and run1["E0"] == run2["E0"]


@pytest.mark.integration
def test_v06_n6_the_summary_entry_agrees_with_the_column_and_nothing_depends_on_it():
    c, run, h = run_with_density_scale(1.0)
    d = ep._alpha_diagnostic(h, c, run)
    a, t = h["condensed_volume_fraction"], h["time_s"]
    assert d["alpha_initial"] == a[0] and d["alpha_final"] == a[-1] and d["alpha_max"] == a.max()
    assert d["alpha_max_time_s"] == t[int(np.argmax(a))]
    assert d["densities_kg_m3"] == {"OX(s)": RHO_S, "OX(L)": RHO_L} and d["species_seen"] == ["OX(L)", "OX(s)"]
    assert "diagnostic" in d["definition"] and "not an error estimate" in d["definition"] and d["convention"].startswith("V0")
    assert not any(k in d for k in ("threshold", "warning", "passed", "limit"))
    # the closure residuals (the pass/fail quantities) do not contain alpha: they are the same four residuals as before
    assert [n for n in ep.HISTORY_COLUMNS if n.startswith("R_")] == ["R_M_kg", "R_B_kg", "R_E_J", "R_V_m3", "R_M_normalised",
                                                                   "R_B_normalised", "R_E_normalised", "R_V_normalised"]


@pytest.mark.fast
def test_v07_n7_the_outflow_enthalpy_is_u_plus_r_eff_t_in_the_v0_convention():
    for T in (1500.0, 3000.0):                                             # solid and liquid condensed phase
        E = Melting().tp(T, 1.0e6, B_E).u_J_per_kg
        rec = recover_state(E, V_FREE, B_E, Melting(), CONDENSED, RECOVERY, T_search_K=SEARCH, starts_K=STARTS)
        M = rec.m_total_kg
        R_eff = RU * rec.n_gas_kmol / M
        u = E / M
        h0 = ei._outflow_enthalpy(E, M, rec.p_Pa, V_FREE)
        assert h0 == pytest.approx(u + R_eff * rec.T_K, rel=VOLUME_TOL)       # H0_out = h of the equilibrium state (V0)
        assert rec.p_Pa * V_FREE / M == pytest.approx(R_eff * rec.T_K, rel=VOLUME_TOL)


# ---- N-9: the linear interpolation of the phase fraction on the plateau --------------------------------------------------------
@dataclasses.dataclass(frozen=True)
class FlatMelting:
    """n- = n+ (n constant), no sensible heat (cv = 0, cp = 0), no pressure term in u: the two ends of the plateau are the exact
    energies of the all-solid and the all-liquid state at T_m, so the linear interpolation is exact."""
    def tp(self, T, p, B_e):
        liquid = T > T_M
        u = YG * U_G0 + YC * (H_L_REF if liquid else H_S_REF)
        return TPState(u_J_per_kg=u, n_gas_per_kg=YG / MG, Yg=YG, condensed_mass_fraction={"OX(L)" if liquid else "OX(s)": YC})


@pytest.mark.fast
def test_v09a_n9a_on_the_plateau_p_does_not_depend_on_phi_and_the_phase_fraction_is_recovered_exactly():
    p_m = 1.0 * (YG / MG) * RU * T_M / V_FREE
    for phi in (0.1, 0.5, 0.9):
        E = 1.0 * (YG * U_G0 + YC * ((1.0 - phi) * H_S_REF + phi * H_L_REF))
        rec = recover_state(E, V_FREE, B_E, FlatMelting(), CONDENSED, RECOVERY, T_search_K=SEARCH, starts_K=STARTS)
        assert rec.on_plateau and rec.T_K == T_M
        assert rec.p_Pa == pytest.approx(p_m, rel=1.0e-9)                    # p does not depend on phi (relative 1e-9)
        assert rec.phi == pytest.approx(phi, abs=1.0e-9)                    # the recovered phi is the original one


@dataclasses.dataclass(frozen=True)
class SlopedMelting:
    """The gas composition differs between the two ends: n(T) = n0 (1 + KAPPA (T - T_m)/T_m), and u carries the energy of that
    composition change (Q_R (n/n0 - 1)), plus the sensible heat of the gas and of the condensed species."""
    def tp(self, T, p, B_e):
        n0 = YG / MG
        n = n0 * (1.0 + KAPPA * (T - T_M) / T_M)
        liquid = T > T_M
        u = YG * (U_G0 + CV_G * (T - 300.0)) + Q_R * (n / n0 - 1.0) + YC * (h_l(T) if liquid else h_s(T))
        return TPState(u_J_per_kg=u, n_gas_per_kg=n, Yg=YG, condensed_mass_fraction={"OX(L)" if liquid else "OX(s)": YC})


KAPPA, Q_R = 0.2, 2.0e5


@pytest.mark.fast
def test_v09b_n9b_with_different_compositions_at_the_two_ends_the_error_is_first_order_in_d():
    phi_true = 0.3
    # the exact state: T = T_m, n = n0 there, so the gas term Q_R (n/n0 - 1) vanishes at the plateau
    E = YG * (U_G0 + CV_G * (T_M - 300.0)) + YC * ((1.0 - phi_true) * h_s(T_M) + phi_true * h_l(T_M))
    errors = {}
    for d in (2.0e-4, 1.0e-4):
        b = _Solver(E, V_FREE, B_E, SlopedMelting(), CONDENSED, RECOVERY).plateau_bounds(TRANSITION, d)
        assert b.state_minus.n_gas_per_kg != b.state_plus.n_gas_per_kg           # the composition differs between the ends
        errors[d] = abs(plateau_phase_fraction(E, b.E_minus, b.E_plus) - phi_true)
    ratio = errors[2.0e-4] / errors[1.0e-4]
    assert errors[1.0e-4] > 1.0e-6                                             # a real error, not round-off
    assert 1.8 <= ratio <= 2.2, (errors, ratio)
