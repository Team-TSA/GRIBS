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
    assert doc["schema"] == "gribs-energy-tp-cache-2" and doc["key"].startswith(path.stem[len(TP_CACHE_PREFIX):])
    path.write_text(path.read_text(encoding="utf-8").replace('"Yg":', '"Yg":0.5,"Yg_":', 1), encoding="utf-8")   # corrupt the file
    with pytest.raises(StateRecoveryError):
        provider(cache=TPCache(tmp_path)).tp(3300.0, 2.0e6, B_E)


@pytest.mark.fast
def test_p09_a_cea_equilibrium_that_does_not_converge_is_an_error():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")                              # CEA warns "CEA_NOT_CONVERGED"
        with pytest.raises(StateRecoveryError, match="did not converge"):
            provider().tp(5.0, 2.0e6, B_E)                           # 5 K: no equilibrium is found (last_error = 8)


_AL2O3_S = {"species": "AL2O3(a)", "phase": "solid", "density_kg_m3": 3990.0, "cp_J_kgK": 1300.0,
            "enthalpy_reference_J_kg": -1.6e7, "reference_temperature_K": 298.15}
_AL2O3_L = {"species": "AL2O3(L)", "phase": "liquid", "density_kg_m3": 3000.0, "cp_J_kgK": 1300.0,
            "enthalpy_reference_J_kg": -1.5e7, "reference_temperature_K": 298.15}


def write_g5_energy_case(tmp_path, *, with_transition):
    """The g5 energy case of p10/p11 (assumed condensed-phase and solid data).  with_transition=False: the one condensed phase
    AL2O3(L) and no transition (F-1 R-3); True: AL2O3(a), AL2O3(L), the declared transition at 2327 K and the tolerance 0.5."""
    import gribs.config.loader  # noqa: F401
    if with_transition:
        condensed = {"volume": "external_density", "latent_heat_check_tolerance_rel": 0.5,
                     "transitions": [{"solid_species": "AL2O3(a)", "liquid_species": "AL2O3(L)", "melting_temperature_K": 2327.0}],
                     "phases": [dict(_AL2O3_S), dict(_AL2O3_L)]}
    else:
        condensed = {"volume": "external_density", "transitions": [], "phases": [{k: v for k, v in _AL2O3_L.items()
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
    Synthetic condensed-phase and solid data (assumed inputs, not data of a propellant).  The melting plateau is NOT run with CEA
    (p11 fixes that limit, report 196)."""
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
def test_p11_the_cea_path_stops_at_t0_with_exit_code_4_when_a_melting_plateau_is_declared(tmp_path):
    """A KNOWN LIMIT (new H number, C17), fixed by this test: with the transition AL2O3(a) -> AL2O3(L) declared the t = 0 latent
    heat check needs states at T_m +- d, and the CEA path does not run them (report 196, measured at T_m +- d, 2200 K and 2500 K):
    the CEA gamma_s of the states there is at most 1 (rejected by the provider), or the TP solve of the solid side raises
    CEA_FORTRAN_ABORT.  (Above T_m, at 3000-3245 K, gamma_s is above 1 and p10 runs; the boundary temperature is not measured.)
    The run stops at t = 0 with exit code 4 and one of those two messages.  The known limit is fixed by this test: when it is
    resolved this test fails and has to be rewritten."""
    import gribs.config.loader
    from gribs.simulation.energy_postprocess import run_energy_case
    outcome = run_energy_case(gribs.config.loader.load_json_configuration(write_g5_energy_case(tmp_path, with_transition=True)))
    assert isinstance(outcome.error, StateRecoveryError) and outcome.unexpected is False
    text = str(outcome.error)
    assert ("CEA returned an isentropic exponent gamma_s" in text) or ("CEA TP solve raised" in text), text
    assert outcome.write_and_report(tmp_path / "out", {}, {}) == 4
