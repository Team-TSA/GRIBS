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
