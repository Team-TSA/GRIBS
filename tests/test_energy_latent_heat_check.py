"""A-6 stage 2, C16: the latent heat of the phase data against the plateau of the provider (PRE_REGISTRATION_C16, L-01..L-07).

L_input is h_liquid(T_m) - h_solid(T_m) from the configuration's constant-cp data; L_calc is the width of the provider's own plateau
per kg of liquid, (E+ - E-) / (M dy_liq) between T_m - d and T_m + d.  No equation of the model uses L.  At t = 0 the run compares
them once (|L_calc - L_input| <= tol max(|L_input|, |L_calc|), tol = energy_model.condensed.latent_heat_check_tolerance_rel) and
rejects a mismatch with exit code 4.  The providers are synthetic (assumed inputs, not data) and are passed as arguments; L_calc is
predicted from their constants alone, before any measurement (PRE_REGISTRATION_C16, section 1-5):

    L_pred(d) = L + [c_mix 2 d + Yg A_P (p+ - p-)] / dy_liq,  c_mix = Yc cp + Yg cv,  p+- = pressure at T_m +- d.

No monkeypatch, no assignment.
"""
import dataclasses
import json
from dataclasses import asdict
from pathlib import Path

import pytest

import gribs.config
import gribs.config.loader
import gribs.simulation.energy_integrator as ei
import gribs.simulation.energy_postprocess as ep
from gribs.config import Config, geometry
from gribs.constants import R_UNIVERSAL
from gribs.energy.balance import ABSOLUTE_ENTHALPY_BASIS
from gribs.energy.config_model import CondensedData, CondensedPhaseSpec, PhaseTransitionSpec, RecoverySettings
from gribs.energy.energy_validation import LATENT_HEAT_TOLERANCE_CODE, PHASE_DATA_CODE
from gribs.energy.exceptions import LatentHeatMismatchError, StateRecoveryError
from gribs.energy.provider import TPState
from gribs.energy.recovery import latent_heat, latent_heat_from_provider
from gribs.errors import ConfigurationError, GribsError

ROOT = Path(__file__).resolve().parents[1]
G1 = json.loads((ROOT / "cases" / "golden" / "g1_reference_constant.json").read_text(encoding="utf-8"))
RU = R_UNIVERSAL

# ---- the analytic provider of L-01 / L-02 (constants as in test_energy_recovery_plateau.Melting; A_P and cv are fields) ----------
M, B_E, T_M = 1.0, {"X": 1.0}, 2327.0
YG, YC, MG, U_G0, CP, TREF, L = 0.6, 0.4, 28.0, -3.0e6, 1300.0, 298.15, 1.0e6
CV_G = 900.0
H_S_REF, H_L_REF = -1.6e7, -1.5e7                                   # L = H_L_REF - H_S_REF = 1.0e6 (same cp, same reference temperature)
RHO_S, RHO_L = 3990.0, 3000.0
V_FREE = M * (YG / MG) * RU * T_M / 3.0e6                          # V0: the gas fills V_free (N1d)
SOLID = CondensedPhaseSpec("OX(s)", "solid", RHO_S, CP, H_S_REF, TREF)
LIQUID = CondensedPhaseSpec("OX(L)", "liquid", RHO_L, CP, H_L_REF, TREF)
TRANSITION = PhaseTransitionSpec("OX(s)", "OX(L)", T_M)
CONDENSED = CondensedData("neglected", (SOLID, LIQUID), (TRANSITION,), 1.0e-3)
RECOVERY = RecoverySettings(phase_fraction_tolerance=1.0e-6, root_tolerance=1.0e-10, plateau_edge_tolerance=1.0e-6)


def h_of(phase, T):
    return phase.enthalpy_reference_J_kg + phase.cp_J_kgK * (T - phase.reference_temperature_K)


@dataclasses.dataclass(frozen=True)
class Plateau:
    """Ideal gas (constant cv, a_p p in the energy) plus one condensed species, solid below T_m and liquid above."""
    a_p: float = 0.0
    cv_g: float = CV_G

    def tp(self, T, p, B_e):
        phase = LIQUID if T > T_M else SOLID
        u = YG * (U_G0 + self.cv_g * (T - 300.0) + self.a_p * p) + YC * h_of(phase, T)
        return TPState(u_J_per_kg=u, n_gas_per_kg=YG / MG, Yg=YG, condensed_mass_fraction={phase.species: YC})


def pressure(T):
    return M * (YG / MG) * RU * T / V_FREE


def l_pred(d, a_p, cv_g=CV_G):
    """The predicted L_calc, from the constants alone (no call of the code under test)."""
    c_mix = YC * CP + YG * cv_g
    dp = pressure(T_M + d) - pressure(T_M - d)
    return L + (c_mix * 2.0 * d + YG * a_p * dp) / YC


def l_calc(provider, d_relative=None):
    return latent_heat_from_provider(TRANSITION, V_FREE, B_E, provider, CONDENSED, RECOVERY, d_relative=d_relative)


D_DEFAULT = max(1.0e-9, 10.0 * RECOVERY.root_tolerance)             # the default d_relative of the code (1e-9)


@pytest.mark.fast
def test_l01_l_calc_equals_the_prediction_and_without_pressure_dependence_equals_l():
    """Tolerances (fixed before the measurement).  Against L_pred: 1e-9 relative: round-off of E+ - E- is about 3e-15 relative
    (E about 6e6 J, the width 4e5 J), and the pressure root has a relative error of at most 1.5e-9 in ln p, which moves L_calc
    by Yg A_P p 1.5e-9 sqrt(2) / Yc = 1e-5 J/kg = 1e-11 relative; a formula without the pressure term is off by 1.09e-6, one
    without the sensible heat by 1.2e-8.  Against L (a_p = 0): 1e-7, 8 times the predicted offset 1.2333e-8 (the sensible heat
    2 c_mix d / Yc = 1.233e-2 J/kg at d = 1e-9 T_m).  N1d (V0, the condensed volume is 0): the pressure difference between the two
    ends is now 2 d p / T_m only (it was 1.09e-6 of L through the condensed volume), so the pressure term is 9.0e-6 J/kg
    (9e-12 relative); the registered predictions below are recomputed from the constants for V0 (the a_p = 0 value is unchanged)."""
    d = D_DEFAULT * T_M
    assert l_pred(d, 0.0) == pytest.approx(1000000.0123331, rel=1e-12)             # the registered predictions (see PRE_REGISTRATION_C16)
    assert l_pred(d, 1.0e-3) == pytest.approx(1000000.0123421, rel=1e-12)          # V0 (the V_c-coupled value was 1000001.0898302882)
    no_pressure = l_calc(Plateau(a_p=0.0))
    assert no_pressure == pytest.approx(l_pred(d, 0.0), rel=1.0e-9)
    assert no_pressure == pytest.approx(L, rel=1.0e-7) and no_pressure > L          # the sensible heat only raises it
    with_pressure = l_calc(Plateau(a_p=1.0e-3))
    assert with_pressure == pytest.approx(l_pred(d, 1.0e-3), rel=1.0e-9)
    assert with_pressure - no_pressure == pytest.approx(YG * 1.0e-3 * (pressure(T_M + d) - pressure(T_M - d)) / YC, rel=1.0e-2)  # the (small) pressure term
    assert latent_heat(TRANSITION, CONDENSED.phases) == pytest.approx(L, rel=1.0e-12)        # L_input of the configuration's data


@pytest.mark.fast
def test_l02_the_change_of_l_calc_with_d_is_the_predicted_slope_of_the_sensible_heat():
    """dL_calc/dd = 2 c_mix / dy_liq (a_p = 0).  c_mix = 1060 and 1600 J/(kg K) (cv 900 and 1800) give 5300 and 8000 J/(kg K).
    Tolerance 1e-6 relative on each slope between neighbouring d_relative in {1e-9, 1e-8, 1e-7, 1e-6}: round-off of one L_calc is
    about 3e-9 J/kg, the smallest change is 0.111 J/kg (1e-9 -> 1e-8), so the slope error is about 3e-8; the 1e-6 leaves 30 times."""
    ds = (1.0e-9, 1.0e-8, 1.0e-7, 1.0e-6)
    slopes = {}
    for cv in (900.0, 1800.0):
        values = [l_calc(Plateau(a_p=0.0, cv_g=cv), d_relative=dr) for dr in ds]
        predicted = 2.0 * (YC * CP + YG * cv) / YC
        pairs = [(values[i + 1] - values[i]) / ((ds[i + 1] - ds[i]) * T_M) for i in range(3)]
        for slope in pairs:
            assert slope == pytest.approx(predicted, rel=1.0e-6), (cv, pairs, predicted)
        for dr, v in zip(ds, values):
            assert v == pytest.approx(l_pred(dr * T_M, 0.0, cv), rel=1.0e-9)
        slopes[cv] = predicted
    assert slopes[900.0] == 5300.0 and slopes[1800.0] == 8000.0


# ---- the provider and the configuration of a whole run (L-03, L-04, L-07; as in test_energy_plateau_integration) ----------------
MW, RUN_YG, RUN_YC, GAMMA = 24.0, 0.7, 0.3, 1.2
RUN_CP, H_SURFACE, T_FLAME = 1300.0, -2.0e6, 3000.0
RUN_H_S, RUN_H_L = -1.6e7, -1.5e7
CV_RUN, R_RUN = RU / MW / (GAMMA - 1.0), RU / MW
RUN_T_M = 7000.0                                                    # above every temperature of the run: it never reaches T_m


def run_h(ref, T):
    return ref + RUN_CP * (T - TREF)


U_G0_RUN = (H_SURFACE - RUN_YC * run_h(RUN_H_L, T_FLAME)) / RUN_YG - CV_RUN * (T_FLAME - 300.0) - R_RUN * T_FLAME


@dataclasses.dataclass(frozen=True)
class RunPlateau:
    enthalpy_convention = ABSOLUTE_ENTHALPY_BASIS
    t_m: float = RUN_T_M

    def tp(self, T, p, B_e):
        liquid = T > self.t_m
        u = RUN_YG * (U_G0_RUN + CV_RUN * (T - 300.0)) + RUN_YC * run_h(RUN_H_L if liquid else RUN_H_S, T)
        return TPState(u_J_per_kg=u, n_gas_per_kg=RUN_YG / MW, Yg=RUN_YG,
                       condensed_mass_fraction={"OX(L)" if liquid else "OX(s)": RUN_YC}, gamma_s=GAMMA)


@dataclasses.dataclass(frozen=True)
class Counting:
    """Records the temperature of every tp call in `log` (a list passed in; appending is not an assignment)."""
    inner: object
    log: list

    @property
    def enthalpy_convention(self):
        return self.inner.enthalpy_convention

    def tp(self, T, p, B_e):
        self.log.append(T)
        return self.inner.tp(T, p, B_e)


def phase_data(species, kind, density, enthalpy):
    return {"species": species, "phase": kind, "density_kg_m3": density, "cp_J_kgK": RUN_CP, "enthalpy_reference_J_kg": enthalpy,
            "reference_temperature_K": TREF}


RECOVERY_DOC = {"phase_fraction_tolerance": 1.0e-6, "root_tolerance": 1.0e-10, "plateau_edge_tolerance": 1.0e-6,
                "temperature_search_K": [200.0, 8000.0], "start_temperatures_K": [600.0, 3000.0]}
FRACTIONS = {"C": 0.3, "H": 0.7}


def energy_doc(tol, t_m=RUN_T_M, h_liquid=RUN_H_L):
    condensed = {"volume": "neglected",
                 "phases": [phase_data("OX(s)", "solid", 3990.0, RUN_H_S), phase_data("OX(L)", "liquid", 3000.0, h_liquid)],
                 "transitions": [{"solid_species": "OX(s)", "liquid_species": "OX(L)", "melting_temperature_K": t_m}]}
    if tol is not None:
        condensed["latent_heat_check_tolerance_rel"] = tol
    return {"initial_fill": {"temperature_K": 300.0, "composition": "main_propellant_products", "element_tolerance": 1.0e-6},
            "wall": {"adiabatic": True},
            "burn_surface": {"form": "a", "grain_far_temperature_K": 293.0,
                             "solid_enthalpy": {"reference_J_kg": H_SURFACE, "reference_temperature_K": 293.0, "cp_J_kgK": 1300.0}},
            "condensed": condensed, "recovery": RECOVERY_DOC}


def run_config(tol, h_liquid=RUN_H_L):
    c0 = gribs.config.loader.configuration_from_document(
        {**G1, "two_phase_model": {**G1["two_phase_model"], "mode": "homogeneous_equilibrium"},
         "energy_model": energy_doc(tol, h_liquid=h_liquid)}).config
    sub = type("ConfigWithEnergy", (Config,), {"energy": c0.energy})
    return sub(**{**asdict(c0), "rtol": 1.0e-8, "R_p": 0.0045, "t_max": 0.01, "blowdown": False})


def run(tol, provider=None, h_liquid=RUN_H_L):
    c = run_config(tol, h_liquid)
    return c, ei.run_model_energy(c, ei.EnergyRunInputs(provider=provider or RunPlateau(), element_fractions=FRACTIONS))


@pytest.mark.integration
def test_l03_a_mismatch_of_l_input_and_l_calc_is_rejected_at_t0_with_exit_code_4(tmp_path, capsys):
    # (a) consistent data: accepted; L_calc is the prediction of this provider (no pressure term: 1e-9 relative)
    c, ok = run(1.0e-3)
    (rec,) = ok["latent_heat_check"]
    d = D_DEFAULT * RUN_T_M
    c_mix = RUN_YC * RUN_CP + RUN_YG * CV_RUN
    assert rec["L_calc_J_per_kg"] == pytest.approx(1.0e6 + 2.0 * c_mix * d / RUN_YC, rel=1.0e-9)
    assert rec["L_input_J_per_kg"] == pytest.approx(1.0e6, rel=1.0e-12) and rec["accepted"] is True
    assert rec["relative_difference"] < 1.0e-3 and rec["tolerance_rel"] == 1.0e-3
    # (b) h_liquid = -1.48e7 gives L_input = 1.2e6 while the provider's plateau is 1.0e6: the relative difference is (1.2e6 - L_calc) / 1.2e6 = 1/6
    with pytest.raises(LatentHeatMismatchError) as info:
        run(0.1, h_liquid=-1.48e7)
    err = info.value
    assert err.code == "E_ENERGY_LATENT_HEAT_MISMATCH" and isinstance(err, StateRecoveryError) and isinstance(err, GribsError)
    text = str(err)
    for part in ("E_ENERGY_LATENT_HEAT_MISMATCH", "L_input", "L_calc", "latent_heat_check_tolerance_rel", "OX(s)", "OX(L)", "7000"):
        assert part in text, (part, text)
    _, edge = run(0.2, h_liquid=-1.48e7)                                    # 1/6 < 0.2: accepted
    assert edge["latent_heat_check"][0]["relative_difference"] == pytest.approx(1.0 / 6.0, rel=1.0e-6)
    # (c) the command line side: a mismatch is exit code 4 with the exception's name
    assert ep.EnergyOutcome(error=err).write_and_report(tmp_path, {}, {}) == 4
    assert "LatentHeatMismatchError" in capsys.readouterr().err


@pytest.mark.integration
def test_l03b_l_calc_that_cannot_be_computed_is_exit_code_4_and_the_config_check_does_not_run_it(tmp_path, capsys):
    """The analytic reference_constant backend has no condensed phase, so with a declared transition y_liq does not increase across
    T_m and L_calc is undefined (a real condition, through run_energy_case).  The configuration check (gribs.config.validate, which
    --validate-config runs) accepts the same configuration: the check needs the provider."""
    doc = {**G1, "two_phase_model": {**G1["two_phase_model"], "mode": "homogeneous_equilibrium"},
           "energy_model": {**energy_doc(0.5, t_m=2327.0), "propellant_formulas": {
               "NH4CLO4(I)": {"N": 1, "H": 4, "Cl": 1, "O": 4}, "C4H6,butadiene": {"C": 4, "H": 6}}}}
    doc["output"]["directory"] = str(tmp_path / "out")
    path = tmp_path / "case.json"
    path.write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")
    loaded = gribs.config.loader.load_json_configuration(path)
    gribs.config.validate(loaded.config)                                   # accepted: nothing is run here
    outcome = ep.run_energy_case(loaded)
    assert isinstance(outcome.error, StateRecoveryError) and not isinstance(outcome.error, LatentHeatMismatchError)
    assert "no liquid increase" in str(outcome.error) and outcome.unexpected is False
    assert outcome.write_and_report(tmp_path / "out", {}, {}) == 4
    assert "Simulation failed: StateRecoveryError" in capsys.readouterr().err


@pytest.mark.integration
def test_l04_the_check_runs_once_at_t0_and_costs_what_a_direct_call_costs():
    """A run that never reaches T_m (T_m = 7000 K): every provider call at T_m +- d is the check's.  Their number equals that of
    one direct call of latent_heat_from_provider (same V_free, B_e) with another counting provider; a check in the right-hand side
    would multiply it by the number of evaluations (the run makes far more than ten)."""
    d = D_DEFAULT * RUN_T_M
    log = []
    c, result = run(1.0e-3, provider=Counting(RunPlateau(), log))
    near_run = [T for T in log if abs(T - RUN_T_M) <= 2.0 * d]
    assert result["integration_budget"]["rhs_calls"] > 10 and len(log) > len(near_run)
    assert not [T for T in log if T >= RUN_T_M - 100.0 and abs(T - RUN_T_M) > 2.0 * d]      # nothing but the check comes near T_m
    direct = []
    cfg = c.energy
    V_free0 = geometry(0.0, c)[3]
    B_e = {k: result["M0"] * v for k, v in FRACTIONS.items()}
    latent_heat_from_provider(cfg.condensed.transitions[0], V_free0, B_e, Counting(RunPlateau(), direct), cfg.condensed, cfg.recovery)
    near_direct = [T for T in direct if abs(T - RUN_T_M) <= 2.0 * d]
    assert len(near_direct) > 0 and len(near_run) == len(near_direct), (len(near_run), len(near_direct))


# ---- the configuration (L-05, L-06) ------------------------------------------------------------------------------------------
def load(energy):
    c = gribs.config.loader.configuration_from_document(
        {**G1, "two_phase_model": {**G1["two_phase_model"], "mode": "homogeneous_equilibrium"}, "energy_model": energy}).config
    gribs.config.validate(c)
    return c


def rejected(energy, code):
    with pytest.raises(ConfigurationError) as info:
        load(energy)
    assert str(info.value).startswith(code), str(info.value)
    return str(info.value)


def without_keys(doc, phase_index, *keys):
    phases = [dict(p) for p in doc["condensed"]["phases"]]
    for k in keys:
        phases[phase_index].pop(k)
    return {**doc, "condensed": {**doc["condensed"], "phases": phases}}


@pytest.mark.fast
def test_l05_the_latent_heat_tolerance_is_required_with_a_transition_and_rejected_without_one():
    assert load(energy_doc(1.0e-3)).energy.condensed.latent_heat_check_tolerance_rel == 1.0e-3
    rejected(energy_doc(None), LATENT_HEAT_TOLERANCE_CODE)                             # a transition, no key
    no_transition = {**energy_doc(None), "condensed": {"volume": "neglected", "phases": [], "transitions": []}}
    assert load(no_transition).energy.condensed.latent_heat_check_tolerance_rel is None
    unused = {**no_transition, "condensed": {**no_transition["condensed"], "latent_heat_check_tolerance_rel": 1.0e-3}}
    assert "no phase transition is declared" in rejected(unused, LATENT_HEAT_TOLERANCE_CODE)
    for bad in (0.0, -1.0e-3):
        with pytest.raises(ConfigurationError, match="latent_heat_check_tolerance_rel: must be > 0"):
            load(energy_doc(bad))


@pytest.mark.fast
def test_l06_cp_enthalpy_and_reference_temperature_are_required_for_the_phases_of_a_transition_and_rejected_for_others():
    keys = ("cp_J_kgK", "enthalpy_reference_J_kg", "reference_temperature_K")
    base = energy_doc(1.0e-3)
    for i in (0, 1):                                                       # solid and liquid of the transition
        for k in keys:
            assert k in rejected(without_keys(base, i, k), PHASE_DATA_CODE)
    other = {"species": "CO2(a)", "phase": "solid", "density_kg_m3": 1500.0}
    alone = {**base, "condensed": {**base["condensed"], "phases": base["condensed"]["phases"] + [other]}}
    assert len(load(alone).energy.condensed.phases) == 3                  # a phase of no transition with density only: accepted
    for k, v in (("cp_J_kgK", 800.0), ("enthalpy_reference_J_kg", -9.0e6), ("reference_temperature_K", 298.15)):
        with_one = {**base, "condensed": {**base["condensed"], "phases": base["condensed"]["phases"] + [{**other, k: v}]}}
        assert k in rejected(with_one, PHASE_DATA_CODE)


@pytest.mark.integration
def test_l07_the_run_records_l_calc_l_input_and_the_relative_difference_for_each_transition():
    _, result = run(1.0e-3)
    (rec,) = result["latent_heat_check"]
    assert set(rec) == {"solid_species", "liquid_species", "melting_temperature_K", "L_calc_J_per_kg", "L_input_J_per_kg",
                        "relative_difference", "tolerance_rel", "accepted"}
    assert (rec["solid_species"], rec["liquid_species"], rec["melting_temperature_K"]) == ("OX(s)", "OX(L)", RUN_T_M)
    assert rec["relative_difference"] == pytest.approx(abs(rec["L_calc_J_per_kg"] - rec["L_input_J_per_kg"])
                                                       / max(rec["L_calc_J_per_kg"], rec["L_input_J_per_kg"]), rel=1.0e-12)
    assert json.loads(json.dumps(result["latent_heat_check"])) == result["latent_heat_check"]      # what summary.json holds
