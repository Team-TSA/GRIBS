"""A-6 stage 2, C15: the phase table and the declared phase transitions must agree (PRE_REGISTRATION_C15, F1-01..F1-09).

A configuration whose condensed phases hold a solid and a liquid of one species without the declared transition is rejected
by the configuration check (exit code 2, before any run).  "One species" is decided from the names alone: the name is split with
the regular expression of energy_validation, and two names are one species only when the bases are equal as strings (the case
is significant).  The values are assumed inputs of the software contract, not data.  No assignment and no monkeypatch.
"""
import json
from pathlib import Path

import pytest

import gribs.config
import gribs.config.loader
from gribs.energy.energy_validation import UNDECLARED_PHASE_TRANSITION_CODE
from gribs.errors import ConfigurationError

ROOT = Path(__file__).resolve().parents[1]
G1 = json.loads((ROOT / "cases" / "golden" / "g1_reference_constant.json").read_text(encoding="utf-8"))
RECOVERY = {"phase_fraction_tolerance": 1.0e-6, "root_tolerance": 1.0e-10, "plateau_edge_tolerance": 1.0e-6,
            "temperature_search_K": [200.0, 8000.0], "start_temperatures_K": [600.0, 3000.0]}
SOLID_H = {"reference_J_kg": -2.0e6, "reference_temperature_K": 293.0, "cp_J_kgK": 1300.0}


def phase(species, kind):
    return {"species": species, "phase": kind, "density_kg_m3": 3000.0, "cp_J_kgK": 1300.0,
            "enthalpy_reference_J_kg": -1.5e7, "reference_temperature_K": 298.15}


def transition(solid, liquid):
    return {"solid_species": solid, "liquid_species": liquid, "melting_temperature_K": 2327.0}


def document(phases, transitions):
    energy = {"initial_fill": {"temperature_K": 300.0, "composition": "main_propellant_products", "element_tolerance": 1.0e-6},
              "wall": {"adiabatic": True},
              "burn_surface": {"form": "a", "grain_far_temperature_K": 293.0, "solid_enthalpy": SOLID_H},
              "condensed": {"volume": "external_density", "phases": phases, "transitions": transitions},
              "recovery": RECOVERY}
    return {**G1, "two_phase_model": {**G1["two_phase_model"], "mode": "homogeneous_equilibrium"}, "energy_model": energy}


def load(phases, transitions):
    """The configuration of the document, checked by gribs.config.validate (the check of the command line, also with --validate-config)."""
    c = gribs.config.loader.configuration_from_document(document(phases, transitions)).config
    gribs.config.validate(c)
    return c


def rejected(phases, transitions):
    with pytest.raises(ConfigurationError) as info:
        load(phases, transitions)
    assert str(info.value).startswith(UNDECLARED_PHASE_TRANSITION_CODE), str(info.value)
    return str(info.value)


AL_S, AL_L = phase("AL2O3(a)", "solid"), phase("AL2O3(L)", "liquid")


@pytest.mark.fast
def test_f101_a_solid_and_a_liquid_of_one_species_without_the_transition_are_rejected_with_the_rule_in_the_message():
    message = rejected([AL_S, AL_L], [])
    assert "E_ENERGY_UNDECLARED_PHASE_TRANSITION" == UNDECLARED_PHASE_TRANSITION_CODE
    assert "^(?P<base>.+)\\((?P<mark>[A-Za-z]{1,3})\\)$" in message and "phase mark (L)" in message and "case-sensitively" in message


@pytest.mark.fast
def test_f102_the_same_two_phases_with_the_transition_are_accepted():
    c = load([AL_S, AL_L], [transition("AL2O3(a)", "AL2O3(L)")])
    assert [p.species for p in c.energy.condensed.phases] == ["AL2O3(a)", "AL2O3(L)"]


@pytest.mark.fast
def test_f103_a_transition_whose_name_differs_from_the_phase_name_only_in_case_does_not_count():
    rejected([AL_S, AL_L], [transition("Al2O3(a)", "AL2O3(L)")])


@pytest.mark.fast
def test_f104_three_phases_of_one_species_are_rejected():
    rejected([AL_S, phase("AL2O3(b)", "solid"), AL_L], [transition("AL2O3(a)", "AL2O3(L)")])


@pytest.mark.fast
def test_f105_two_solids_of_one_species_a_solid_solid_transition_are_rejected():
    rejected([AL_S, phase("AL2O3(b)", "solid")], [])


@pytest.mark.fast
def test_f106_a_phase_mark_that_does_not_match_the_phase_field_is_rejected_in_both_directions():
    rejected([phase("X(L)", "solid")], [])
    rejected([phase("Y(cr)", "liquid")], [])
    rejected([phase("Z", "liquid")], [])                                   # a liquid without a mark


@pytest.mark.fast
def test_f107_phases_of_different_species_need_no_transition():
    c = load([AL_S, phase("AL(L)", "liquid")], [])
    assert len(c.energy.condensed.phases) == 2 and c.energy.condensed.transitions == ()


@pytest.mark.fast
def test_f108_duplicated_names_and_transitions_that_are_not_a_solid_and_a_liquid_of_one_species_are_rejected():
    rejected([AL_S, AL_S], [])                                             # R-1: the same name twice
    rejected([AL_S, phase("AL(L)", "liquid")], [transition("AL2O3(a)", "AL(L)")])      # R-6: two species
    rejected([AL_S, AL_L], [transition("AL2O3(a)", "AL2O3(L)"), transition("AL2O3(a)", "NO_SUCH(L)")])   # R-6: unknown name


@pytest.mark.fast
def test_f109_names_that_differ_only_in_case_are_different_species():
    c = load([AL_S, phase("Al2O3(L)", "liquid")], [])                      # not one species: no transition is needed
    assert [p.species for p in c.energy.condensed.phases] == ["AL2O3(a)", "Al2O3(L)"]
    c = load([phase("CO(a)", "solid"), phase("Co(L)", "liquid")], [])      # CO is not Co (B-7); casefold would make them one species
    assert c.energy.condensed.transitions == ()
    rejected([AL_S, AL_L], [])                                             # the control: the same case is one species
