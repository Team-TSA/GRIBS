"""A-6 stage 2, C7: the CEA TP provider and its cache (contract 164 sections 3-A, 5; F-21, F-24; PRE_REGISTRATION_C7).

These tests need the official NASA CEA package: without it the whole module is skipped, and that is recorded as "not
executable (reason: cea not installed)", never as passed.  The propellant is the one of the existing golden case g5 (AP / Al /
HTPB-like); the tests check the PROVIDER's thermodynamic output and its bookkeeping, not an energy-mode run of any propellant.
No assignment, no monkeypatch.
"""
import json
import math
import os
import warnings
from pathlib import Path

import pytest

import gribs.thermo

cea = pytest.importorskip("cea")

from gribs.energy.cea_tp import TP_CACHE_PREFIX, CEATPProvider, TPCache  # noqa: E402
from gribs.energy.exceptions import StateRecoveryError  # noqa: E402
from gribs.energy.provider import TPState  # noqa: E402

NAMES = ["NH4CLO4(I)", "AL(cr)", "C4H6,butadiene"]
WT = [68.0, 18.0, 14.0]
ATOMIC = {"N": 14.0067, "H": 1.00794, "Cl": 35.453, "O": 15.9994, "Al": 26.9815, "C": 12.0107}
FORMULAS = {"NH4CLO4(I)": {"N": 1, "H": 4, "Cl": 1, "O": 4}, "AL(cr)": {"Al": 1}, "C4H6,butadiene": {"C": 4, "H": 6}}


def element_fractions():
    out = {e: 0.0 for e in ATOMIC}
    for name, wt in zip(NAMES, WT):
        mass = sum(ATOMIC[e] * n for e, n in FORMULAS[name].items())
        for e, n in FORMULAS[name].items():
            out[e] += wt / 100.0 * ATOMIC[e] * n / mass
    return out


ELEMENTS = element_fractions()
B_E = {e: 2.0 * f for e, f in ELEMENTS.items()}           # 2 kg of the main-propellant products


def provider(cache=None):
    return CEATPProvider(NAMES, WT, ELEMENTS, cache=cache)


def hp_temperature_and_enthalpy(p_Pa):
    """The HP equilibrium temperature and the reactant enthalpy h_c [J/kg], computed with the CEA objects directly, as the
    legacy backend does (the HP input is h_c / cea.R)."""
    import numpy as np
    module = gribs.thermo.import_official_cea_module()
    gribs.thermo.pin_bundled_cea_database(module, False)
    reactants = module.Mixture(NAMES)
    products = module.Mixture(NAMES, products_from_reactants=True)
    solver = module.EqSolver(products, reactants=reactants, transport=False, ions=False, trace=1.0e-10, insert=[])
    solution = module.EqSolution(solver)
    weights = np.array(WT, dtype=float)
    h_c = float(reactants.calc_property(module.ENTHALPY, weights, np.array([298.15] * 3)))
    solver.solve(solution, module.HP, h_c / float(module.R), p_Pa / 1.0e5, weights)
    assert bool(solution.converged)
    return float(solution.T), h_c


@pytest.mark.fast
def test_p01_the_provider_returns_the_state_and_the_derived_quantities_and_rejects_another_composition():
    prov = provider()
    T, p = 3300.0, 2.0e6
    state = prov.tp(T, p, B_E)
    assert isinstance(state, TPState) and state.converged and state.u_J_per_kg < 0.0
    assert 0.0 < state.Yg < 1.0 and state.n_gas_per_kg > 0.0 and math.isclose(sum(state.condensed_mass_fraction.values()), 1.0 - state.Yg, abs_tol=1e-9)
    d = prov.derived(T, p, B_E)
    assert set(d) == {"n_gas_kmol", "m_total_kg", "m_gas_kg", "Yg", "M_CEA", "M_gas", "R_eff", "R_gas"}
    assert d["m_total_kg"] == pytest.approx(2.0) and d["n_gas_kmol"] == pytest.approx(2.0 * state.n_gas_per_kg)
    assert d["m_gas_kg"] == pytest.approx(2.0 * state.Yg)
    other = {**B_E, "C": B_E["C"] * 1.1}                          # a composition that is not the provider's
    with pytest.raises(StateRecoveryError, match="outside stage 2"):
        prov.tp(T, p, other)


@pytest.mark.fast
def test_p02_the_f21_numbers_of_the_g5_case_at_2_mpa():
    T_hp, _ = hp_temperature_and_enthalpy(2.0e6)
    d = provider().derived(T_hp, 2.0e6, B_E)
    assert d["R_eff"] == pytest.approx(309.31, rel=1e-3)           # F-21 (g5, 2 MPa): R_eff = 309.31
    assert d["R_gas"] == pytest.approx(437.73, rel=1e-3)           # Ru / M_gas = 437.73
    assert d["Yg"] == pytest.approx(0.7066, rel=2e-3)


@pytest.mark.fast
def test_p03_m_cea_and_m_gas_are_two_different_fields():
    d = provider().derived(3300.0, 2.0e6, B_E)
    assert "M_CEA" in d and "M_gas" in d and d["M_CEA"] != d["M_gas"]
    assert d["M_gas"] == pytest.approx(d["Yg"] * d["M_CEA"], rel=1e-12)
    assert d["M_CEA"] * (d["n_gas_kmol"] / d["m_total_kg"]) == pytest.approx(1.0, rel=1e-9)     # M_CEA = 1 / n_gas


@pytest.mark.fast
def test_p04_r_eff_is_yg_times_r_gas_for_the_cea_output():
    d = provider().derived(3300.0, 2.0e6, B_E)
    assert d["R_eff"] == pytest.approx(d["Yg"] * d["R_gas"], rel=1e-12)
    assert d["R_eff"] == pytest.approx(8314.51 / d["M_CEA"], rel=1e-9)


@pytest.mark.fast
def test_p05_the_tp_enthalpy_at_the_hp_temperature_is_the_reactant_enthalpy():
    T_hp, h_c = hp_temperature_and_enthalpy(2.0e6)
    detail = provider().tp_detail(T_hp, 2.0e6, B_E)
    assert detail.h_J_per_kg == pytest.approx(h_c, rel=1e-6)       # the kJ -> J conversion and the absolute basis agree
    assert detail.u_J_per_kg < detail.h_J_per_kg                   # u = h - p v, with p v > 0
    assert provider().enthalpy_convention == "absolute_formation_included"


@pytest.mark.fast
def test_p05b_n7_r_h_the_enthalpy_of_the_provider_is_u_plus_the_v0_pv_term():
    """r_h = (h - u - p V_free/M)/h with V_free the V0 volume of the gas (M n Ru T/p): the provider's h agrees with u + n Ru T,
    so that H0_out = E_cv/M + p V_free/M is the enthalpy of the equilibrium state (N1d, decision of the 83rd record, item 1 (a)).
    A test-only quantity: it is not in any output and not in any equation.  Tolerance = the volume tolerance of the recovery,
    1e-8 (a choice of the N1d commit)."""
    mass = math.fsum(B_E.values())
    for T, p in ((3300.0, 2.0e6), (2326.0, 2.0e6), (2328.0, 2.0e6), (3000.0, 1.0e6)):
        d = provider().tp_detail(T, p, B_E)
        v_free = mass * d.n_gas_per_kg * 8314.51 * T / p
        r_h = (d.h_J_per_kg - d.u_J_per_kg - p * v_free / mass) / d.h_J_per_kg
        assert abs(r_h) <= 1.0e-8, (T, p, r_h)


@pytest.mark.fast
def test_p06_the_cache_key_depends_on_the_weights_the_element_order_and_the_enthalpy_convention():
    identity = provider().cache_identity()
    base = TPCache.key(identity, 3300.0, 2.0e6)
    assert TPCache.key(identity, 3300.0, 2.0e6) == base                                       # deterministic
    assert TPCache.key({**identity, "reactant_wt_percent": [68.0, 18.5, 13.5]}, 3300.0, 2.0e6) != base
    assert TPCache.key({**identity, "element_order": list(reversed(identity["element_order"]))}, 3300.0, 2.0e6) != base
    assert TPCache.key({**identity, "enthalpy_convention": "sensible_only"}, 3300.0, 2.0e6) != base
    assert TPCache.key(identity, 3300.0, 2.1e6) != base and TPCache.key(identity, 3301.0, 2.0e6) != base


@pytest.mark.fast
def test_p07_the_p_only_cache_is_neither_read_nor_written(tmp_path):
    legacy = tmp_path / "thermo_cea_python_0123456789abcdef01234567.json"
    legacy.write_text("{not json", encoding="utf-8")
    before = legacy.read_bytes()
    os.chmod(legacy, 0o000)                                        # reading it would raise PermissionError
    try:
        prov = provider(cache=TPCache(tmp_path))
        prov.tp(3300.0, 2.0e6, B_E)
        prov.tp(3300.0, 2.0e6, B_E)
    finally:
        os.chmod(legacy, 0o600)
    assert legacy.read_bytes() == before
    assert sorted(p.name for p in tmp_path.iterdir() if p.name != legacy.name)
    assert all(p.name.startswith(TP_CACHE_PREFIX) for p in tmp_path.iterdir() if p.name != legacy.name)


@pytest.mark.fast
def test_p08_a_tp_cache_round_trip_returns_the_same_state(tmp_path):
    first_cache, second_cache = TPCache(tmp_path), TPCache(tmp_path)
    a = provider(cache=first_cache).tp(3300.0, 2.0e6, B_E)
    assert first_cache.stores == 1 and first_cache.hits == 0
    b = provider(cache=second_cache).tp(3300.0, 2.0e6, B_E)          # a new provider reads the file the first one wrote
    assert second_cache.hits == 1 and second_cache.stores == 0
    assert a == b                                                    # identical, bit for bit (JSON floats round-trip)
    assert a.gamma_s > 1.0 and b.gamma_s == a.gamma_s                # the equilibrium gamma_s is part of the state (C11, B-14)
    (path,) = list(tmp_path.glob(TP_CACHE_PREFIX + "*.json"))
    doc = json.loads(path.read_text(encoding="utf-8"))
    assert doc["schema"] == "gribs-energy-tp-cache-3" and doc["key"].startswith(path.stem[len(TP_CACHE_PREFIX):])
    path.write_text(path.read_text(encoding="utf-8").replace('"Yg":', '"Yg":0.5,"Yg_":', 1), encoding="utf-8")   # corrupt the file
    with pytest.raises(StateRecoveryError):
        provider(cache=TPCache(tmp_path)).tp(3300.0, 2.0e6, B_E)


@pytest.mark.fast
def test_p09_a_cea_equilibrium_that_does_not_converge_is_an_error():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")                              # CEA warns "CEA_NOT_CONVERGED"
        with pytest.raises(StateRecoveryError, match="did not converge"):
            provider().tp(5.0, 2.0e6, B_E)                           # 5 K: no equilibrium is found (last_error = 8)


# --- C18 (report 198): the provider solves every state on a new EqSolution; the result does not depend on the history ---------
_Y_T_M_MINUS, _Y_T_M_PLUS = 2326.999998, 2327.000002                      # T_m = 2327 K of Al2O3 -+ d (the points of report 198)
_Y_PGRID = [1.5e6 + i * 1.5e5 for i in range(11)]
_Y_POINTS = ([(_Y_T_M_MINUS, p) for p in _Y_PGRID] + [(2200.0, p) for p in _Y_PGRID] + [(_Y_T_M_PLUS, p) for p in (1.5e6, 2.25e6, 3.0e6)])


def _y_close(a, b, tol):
    """The four values of a TP state agree to the relative tolerance tol."""
    for name in ("gamma_s", "u_J_per_kg", "Yg", "n_gas_per_kg"):
        x, y = getattr(a, name), getattr(b, name)
        assert abs(x - y) <= tol * max(abs(x), abs(y)), (name, x, y)


@pytest.mark.fast
def test_y01_the_states_do_not_depend_on_the_order_in_which_the_provider_solved_before():
    """Report 198: a reused EqSolution failed at the change from the solid to the liquid Al2O3 and, after a failure, gave a
    gamma_s below 1.  The 25 states of the report (T_m -+ d and 2200 K, solid side first, then the liquid) are solved in order,
    in the reverse order and each on a new provider: all succeed, gamma_s > 1 (B-14), and the values agree to 1e-7 (the old
    behaviour gave a failure or a difference of up to 0.19 in gamma_s; the remaining difference of the solver is of the order
    of 1e-10)."""
    reference = [provider().tp(T, p, B_E) for T, p in _Y_POINTS]
    forward_provider, reverse_provider = provider(), provider()
    forward = [forward_provider.tp(T, p, B_E) for T, p in _Y_POINTS]
    reverse = [reverse_provider.tp(T, p, B_E) for T, p in reversed(_Y_POINTS)][::-1]
    assert len(reference) == len(forward) == len(reverse) == 25
    for state in reference + forward + reverse:
        assert state.gamma_s > 1.0
    for ref, fwd, rev in zip(reference, forward, reverse):
        _y_close(fwd, ref, 1.0e-7)
        _y_close(rev, ref, 1.0e-7)


@pytest.mark.fast
def test_y02_the_state_after_the_change_from_the_solid_to_the_liquid_is_the_state_of_a_new_provider():
    """The point after the change from the solid (T_m - d, 3e6 Pa) to the liquid (T_m + d, 1.5e6 Pa) gave gamma_s = 0.9994
    with the reused solution (report 198; the change itself failed).  Now it is the value of a new provider, and above 1."""
    used = provider()
    used.tp(_Y_T_M_MINUS, 3.0e6, B_E)
    used.tp(_Y_T_M_PLUS, 1.5e6, B_E)                                     # raised CEA_FORTRAN_ABORT with the reused solution
    after = used.tp(_Y_T_M_PLUS, 2.25e6, B_E)
    assert after.gamma_s > 1.0
    _y_close(after, provider().tp(_Y_T_M_PLUS, 2.25e6, B_E), 1.0e-7)


@pytest.mark.fast
def test_y03_one_state_after_two_different_histories_is_the_same_to_1e_9():
    """The same (T, p) = (2200 K, 1.5e6 Pa) on a provider that solved it first and on one that solved a liquid state before."""
    first, second = provider(), provider()
    second.tp(_Y_T_M_PLUS, 3.0e6, B_E)
    _y_close(second.tp(2200.0, 1.5e6, B_E), first.tp(2200.0, 1.5e6, B_E), 1.0e-9)


_AL2O3_S = {"species": "AL2O3(a)", "phase": "solid", "density_kg_m3": 3990.0, "cp_J_kgK": 1300.0,
            "enthalpy_reference_J_kg": -1.6e7, "reference_temperature_K": 298.15}
_AL2O3_L = {"species": "AL2O3(L)", "phase": "liquid", "density_kg_m3": 3000.0, "cp_J_kgK": 1300.0,
            "enthalpy_reference_J_kg": -1.5e7, "reference_temperature_K": 298.15}


def write_g5_energy_case(tmp_path, *, with_transition):
    """The g5 energy case of p10/p11 (assumed condensed-phase and solid data).  with_transition=False: the one condensed phase
    AL2O3(L) and no transition (F-1 R-3); True: AL2O3(a), AL2O3(L), the declared transition at 2327 K and the tolerance 0.5."""
    import gribs.config.loader  # noqa: F401
    if with_transition:
        condensed = {"volume": "neglected", "latent_heat_check_tolerance_rel": 0.5,
                     "transitions": [{"solid_species": "AL2O3(a)", "liquid_species": "AL2O3(L)", "melting_temperature_K": 2327.0}],
                     "phases": [dict(_AL2O3_S), dict(_AL2O3_L)]}
    else:
        condensed = {"volume": "neglected", "transitions": [], "phases": [{k: v for k, v in _AL2O3_L.items()
                                                               if k not in ("cp_J_kgK", "enthalpy_reference_J_kg", "reference_temperature_K")}]}
    root = Path(__file__).resolve().parents[1]
    doc = json.loads((root / "cases" / "golden" / "g5_high_condensed_ap_al.json").read_text(encoding="utf-8"))
    doc["output"]["directory"] = str(tmp_path / "g5_energy_out")
    doc["solver"].update(burning_time_limit_s=0.002, relative_tolerance=1.0e-7)
    doc["two_phase_model"]["mode"] = "homogeneous_equilibrium"
    doc["energy_model"] = {
        "initial_fill": {"temperature_K": 3000.0, "composition": "main_propellant_products", "element_tolerance": 1.0e-6},
        "wall": {"adiabatic": True},
        "burn_surface": {"form": "a", "grain_far_temperature_K": 293.0,
                         "solid_enthalpy": {"reference_J_kg": -1.4e6, "reference_temperature_K": 293.0, "cp_J_kgK": 1300.0}},
        "condensed": condensed,
        "recovery": {"phase_fraction_tolerance": 1.0e-6, "root_tolerance": 1.0e-10, "plateau_edge_tolerance": 1.0e-6,
                     "temperature_search_K": [500.0, 5000.0], "start_temperatures_K": [3000.0]},
        "propellant_formulas": {"NH4CLO4(I)": {"N": 1, "H": 4, "Cl": 1, "O": 4}, "AL(cr)": {"Al": 1},
                                "C4H6,butadiene": {"C": 4, "H": 6}},
    }
    path = tmp_path / "g5_energy.json"
    path.write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")
    return path


@pytest.mark.integration
def test_p10_the_wired_cea_provider_runs_an_energy_case_to_completion_with_gamma_above_one(tmp_path):
    """C12 smoke run (B-18 (3)), in the range above T_m where Al2O3 is liquid: the energy case of the golden g5 propellant with the
    CEA provider, through run_energy_case, with ONE condensed phase (AL2O3(L)) and no transition.  It passes on completion,
    gamma_s > 1 in every state and the closure residuals only; nothing is compared with the legacy model.  The run stays at
    T = 3000-3245 K (above T_m = 2327 K), where liquid Al2O3 is present and the CEA gamma_s is above 1.  If T fell below T_m a solid
    would appear and the run would stop with exit code 4 (a condensed species not in the phases, the gap of 191 section 6).
    Synthetic condensed-phase and solid data (assumed inputs, not data of a propellant).  A run through the melting plateau with the CEA
    provider is not tested (p11 only checks the latent heat at t = 0; report 198)."""
    import gribs.config.loader
    from gribs.simulation.energy_postprocess import CLOSURE_TOLERANCE, run_energy_case
    path = write_g5_energy_case(tmp_path, with_transition=False)
    outcome = run_energy_case(gribs.config.loader.load_json_configuration(path))
    assert outcome.error is None
    gamma = outcome.energy_model["gamma_s"]
    assert 1.0 < gamma["min"] <= gamma["max"] and outcome.energy_model["provenance"]["provider"] == "CEATPProvider"
    assert all(g > 1.0 for g in outcome.history["gamma_s"])
    assert outcome.results["closure"]["passed"] is True and outcome.run_passed is True
    assert all(v <= CLOSURE_TOLERANCE for v in outcome.results["closure"]["max_normalised"].values())
    # C12b: a provider that cannot be constructed is a backend initialisation failure (exit code 3), not a run failure
    from dataclasses import asdict
    from gribs.config import Config
    from gribs.energy.wiring import build_energy_inputs
    from gribs.errors import ThermochemistryError
    c = gribs.config.loader.load_json_configuration(path).config
    sub = type("ConfigWithEnergy", (Config,), {"energy": c.energy, "reactants_spec_list": c.reactants_spec_list})
    broken = sub(**{**asdict(c), "cea_py_products_from_reactants": False, "cea_py_product_species": ("NO_SUCH_SPECIES(zz)",)})
    with pytest.raises(ThermochemistryError):
        build_energy_inputs(broken)


@pytest.mark.integration
def test_p11_the_cea_path_passes_the_latent_heat_check_at_t0(tmp_path):
    """C18 (report 198; it replaces the C16 test that fixed a stop at t = 0): with the transition AL2O3(a) -> AL2O3(L) declared
    (assumed latent heat 1.0e6 J/kg, tolerance 0.5) the t = 0 latent heat check needs the states at T_m +- d, and the CEA provider
    now gives them (every state is solved on a new EqSolution).  The run has the settings of p10 (3000 K, 0.002 s): it does NOT
    cross T_m, so this test confirms only that the check passes with the CEA provider and that the run completes above T_m with
    gamma_s > 1 and the closure residuals.  A run through the melting plateau with the CEA provider is not tested.  L_calc is not
    asserted (the condensed-phase data are assumed inputs, not data of a propellant)."""
    import gribs.config.loader
    from gribs.simulation.energy_postprocess import CLOSURE_TOLERANCE, run_energy_case
    outcome = run_energy_case(gribs.config.loader.load_json_configuration(write_g5_energy_case(tmp_path, with_transition=True)))
    assert outcome.error is None and outcome.unexpected is False and outcome.run_passed is True
    assert outcome.results["closure"]["passed"] is True
    assert all(v <= CLOSURE_TOLERANCE for v in outcome.results["closure"]["max_normalised"].values())
    assert all(g > 1.0 for g in outcome.history["gamma_s"])
    (record,) = outcome.energy_model["latent_heat_check"]
    assert record["L_input_J_per_kg"] == 1.0e6 and math.isfinite(record["L_calc_J_per_kg"])
    assert record["accepted"] is True and record["relative_difference"] <= record["tolerance_rel"] == 0.5
