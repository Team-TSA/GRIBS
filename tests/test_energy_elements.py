"""A-6 stage 2, C9: element mass fractions from explicit reactant formulas (contract 164; PRE_REGISTRATION_C9).

The propellant is the one of the existing golden case g5 (AP / Al / HTPB-like); the atomic masses below are written out here
independently of the module's table (they are the values of the C7 test).  The tests check the bookkeeping of the formulas, not
an energy-mode run.  No assignment, no monkeypatch, no parametrize (the patch-site census stays as it is).
"""
import json
import math
from dataclasses import asdict, fields
from pathlib import Path

import pytest

import gribs.config.loader
from gribs.config import Config
from gribs.energy.config_model import EnergyModeSettings
from gribs.energy.elements import ATOMIC_MASS_KG_PER_KMOL, check_formulas, element_mass_fractions, reactant_element_rows
from gribs.energy.initial import ElementBasis, check_initial_fill
from gribs.errors import ConfigurationError

ROOT = Path(__file__).resolve().parents[1]
G1 = json.loads((ROOT / "cases" / "golden" / "g1_reference_constant.json").read_text(encoding="utf-8"))

NAMES = ["NH4CLO4(I)", "AL(cr)", "C4H6,butadiene"]
WT = {"NH4CLO4(I)": 68.0, "AL(cr)": 18.0, "C4H6,butadiene": 14.0}
FORMULAS = {"NH4CLO4(I)": {"N": 1, "H": 4, "Cl": 1, "O": 4}, "AL(cr)": {"Al": 1}, "C4H6,butadiene": {"C": 4, "H": 6}}
ATOMIC = {"N": 14.0067, "H": 1.00794, "Cl": 35.453, "O": 15.9994, "Al": 26.9815, "C": 12.0107}

FILL = {"temperature_K": 300.0, "composition": "main_propellant_products", "element_tolerance": 1.0e-6}
SOLID_H = {"reference_J_kg": -2.0e6, "reference_temperature_K": 293.0, "cp_J_kgK": 1300.0}
SURFACE = {"form": "a", "grain_far_temperature_K": 293.0, "solid_enthalpy": SOLID_H}
CONDENSED = {"volume": "neglected",
             "phases": [{"species": "Al2O3(a)", "phase": "solid", "density_kg_m3": 3990.0}],
             "transitions": []}
RECOVERY = {"phase_fraction_tolerance": 1.0e-6, "root_tolerance": 1.0e-10, "plateau_edge_tolerance": 1.0e-6, "temperature_search_K": [200.0, 8000.0], "start_temperatures_K": [600.0, 3000.0]}
ENERGY = {"initial_fill": FILL, "wall": {"adiabatic": True}, "burn_surface": SURFACE, "condensed": CONDENSED,
          "recovery": RECOVERY}


def expected_element_fractions():
    """The element mass fractions of the mixture, written out from the atomic masses above (independent of the module)."""
    out = {e: 0.0 for e in ATOMIC}
    for name, wt in WT.items():
        mass = sum(ATOMIC[e] * n for e, n in FORMULAS[name].items())
        for e, n in FORMULAS[name].items():
            out[e] += wt / 100.0 * ATOMIC[e] * n / mass
    return out


def document(energy):
    return {**G1, "two_phase_model": {**G1["two_phase_model"], "mode": "homogeneous_equilibrium"}, "energy_model": energy}


def raises_configuration_error(call):
    with pytest.raises(ConfigurationError):
        call()


@pytest.mark.fast
def test_f1_the_element_fractions_of_g5_equal_the_independent_computation():
    got, want = element_mass_fractions(WT, FORMULAS), expected_element_fractions()
    assert set(got) == set(want)
    for element in want:
        assert got[element] == pytest.approx(want[element], rel=1.0e-12, abs=0.0)
    for element in ATOMIC:                                              # the table has the values of the C7 test
        assert ATOMIC_MASS_KG_PER_KMOL[element] == ATOMIC[element]


@pytest.mark.fast
def test_f2_every_row_and_the_mixture_sum_to_one_and_the_weight_percent_is_normalised():
    rows = reactant_element_rows(FORMULAS)
    assert set(rows) == set(FORMULAS)
    for name, row in rows.items():
        assert math.fsum(row.values()) == pytest.approx(1.0, abs=1.0e-15)
        assert set(row) == set(FORMULAS[name])
    fractions = element_mass_fractions(WT, FORMULAS)
    assert math.fsum(fractions.values()) == pytest.approx(1.0, abs=1.0e-15)
    scaled = element_mass_fractions({k: 10.0 * v for k, v in WT.items()}, FORMULAS)      # 680 / 180 / 140
    for element, value in fractions.items():
        assert scaled[element] == pytest.approx(value, rel=1.0e-14, abs=0.0)


@pytest.mark.fast
def test_f3_a_missing_or_extra_formula_and_a_bad_weight_percent_are_rejected():
    without_al = {k: v for k, v in FORMULAS.items() if k != "AL(cr)"}
    raises_configuration_error(lambda: element_mass_fractions(WT, without_al))                 # reactant without a formula
    raises_configuration_error(lambda: element_mass_fractions(WT, {**FORMULAS, "H2O": {"H": 2, "O": 1}}))   # formula without a reactant
    raises_configuration_error(lambda: element_mass_fractions({}, FORMULAS))
    for bad in (-1.0, math.nan, math.inf, True):
        raises_configuration_error(lambda bad=bad: element_mass_fractions({**WT, "AL(cr)": bad}, FORMULAS))
    raises_configuration_error(lambda: element_mass_fractions({k: 0.0 for k in WT}, FORMULAS))   # no mass


@pytest.mark.fast
def test_f4_element_symbols_are_case_sensitive_and_unknown_symbols_and_bad_counts_are_rejected():
    for wrong in ("CL", "CO", "cl", "AL", "co", "Xx"):
        raises_configuration_error(lambda wrong=wrong: check_formulas({"R": {wrong: 1}}))
    check_formulas({"R": {"Cl": 1}})
    check_formulas({"R": {"Co": 1}})
    assert ATOMIC_MASS_KG_PER_KMOL["Cl"] != ATOMIC_MASS_KG_PER_KMOL["Co"]
    monoxide = reactant_element_rows({"CO": {"C": 1, "O": 1}})["CO"]                  # carbon monoxide: C and O
    cobalt = reactant_element_rows({"Co": {"Co": 1}})["Co"]                             # cobalt: the element Co
    assert set(monoxide) == {"C", "O"} and set(cobalt) == {"Co"} and cobalt["Co"] == 1.0
    for bad in (0, -1, 0.0, math.nan, math.inf, True, "1", None):
        raises_configuration_error(lambda bad=bad: check_formulas({"R": {"Cl": bad}}))
    raises_configuration_error(lambda: check_formulas({"R": {}}))
    raises_configuration_error(lambda: check_formulas({}))
    raises_configuration_error(lambda: check_formulas({"": {"Cl": 1}}))


@pytest.mark.fast
def test_f5_the_formulas_key_is_optional_parsed_and_does_not_change_the_configuration_fields():
    legacy_keys = {f.name for f in fields(Config)}
    plain = gribs.config.loader.configuration_from_document(document(ENERGY)).config
    assert isinstance(plain.energy, EnergyModeSettings) and plain.energy.propellant_formulas is None
    assert set(asdict(plain)) == legacy_keys == set(asdict(gribs.config.loader.configuration_from_document(G1).config))
    given = gribs.config.loader.configuration_from_document(document({**ENERGY, "propellant_formulas": FORMULAS})).config
    assert given.energy.propellant_formulas == {k: {e: float(n) for e, n in v.items()} for k, v in FORMULAS.items()}
    assert set(asdict(given)) == legacy_keys and "energy" not in asdict(given) and "propellant_formulas" not in asdict(given)
    for bad in ({"R": {"CL": 1}}, {"R": {"Cl": 0}}, {"R": {}}, {}, [], {"R": [1]}):
        with pytest.raises(ConfigurationError):
            gribs.config.loader.configuration_from_document(document({**ENERGY, "propellant_formulas": bad}))


@pytest.mark.fast
def test_f6_the_rows_make_an_element_basis_that_accepts_the_fill_of_the_mixture():
    basis = ElementBasis.from_mapping(reactant_element_rows(FORMULAS))
    fill = element_mass_fractions(WT, FORMULAS)
    weights = check_initial_fill(fill, fill, basis, element_tolerance=1.0e-9)
    assert set(weights) == set(FORMULAS)
    for name, wt in WT.items():
        assert weights[name] == pytest.approx(wt / 100.0, abs=1.0e-9)
