"""A-6 stage 2, C3: the optional energy_model section (contract 164; reports 185, 188, PRE_REGISTRATION_C3).

The values below are ASSUMED INPUTS for the software contract, not data (sources: PRE_REGISTRATION_C3.md, held items
H-01..H-08).  The tests use no assignment and no monkeypatch (the patch-site census of test_a6_patch_targets.py stays
as it is): documents are built by merging dictionaries, and a configuration with changed fields is built with
a Config subclass that carries the settings as a class attribute.
"""
import json
from dataclasses import asdict, fields, replace
from pathlib import Path

import pytest

import gribs.config
import gribs.config.loader
import gribs.simulation.integrator
from gribs.config import Config
from gribs.energy.config_model import EnergyModeSettings, WallHeat
from gribs.errors import ConfigurationError

ROOT = Path(__file__).resolve().parents[1]
G1 = json.loads((ROOT / "cases" / "golden" / "g1_reference_constant.json").read_text(encoding="utf-8"))

FILL = {"temperature_K": 300.0, "composition": "main_propellant_products", "element_tolerance": 1.0e-6}
WALL = {"adiabatic": True}
SOLID_H = {"reference_J_kg": -2.0e6, "reference_temperature_K": 293.0, "cp_J_kgK": 1300.0}
SURFACE = {"form": "a", "grain_far_temperature_K": 293.0, "solid_density_kg_m3": 1800.0, "solid_enthalpy": SOLID_H}
CONDENSED = {
    "volume": "external_density",
    "phases": [
        {"species": "Al2O3(a)", "phase": "solid", "density_kg_m3": 3990.0, "cp_J_kgK": 1300.0,
         "enthalpy_reference_J_kg": -1.6e7, "reference_temperature_K": 298.15},
        {"species": "Al2O3(L)", "phase": "liquid", "density_kg_m3": 3000.0, "cp_J_kgK": 1300.0,
         "enthalpy_reference_J_kg": -1.5e7, "reference_temperature_K": 298.15}],
    "transitions": [{"solid_species": "Al2O3(a)", "liquid_species": "Al2O3(L)", "melting_temperature_K": 2327.0}]}
RECOVERY = {"phase_fraction_tolerance": 1.0e-6, "root_tolerance": 1.0e-10, "plateau_edge_tolerance": 1.0e-6}
ENERGY = {"initial_fill": FILL, "wall": WALL, "burn_surface": SURFACE, "condensed": CONDENSED, "recovery": RECOVERY}


def without(d, *keys):
    return {k: v for k, v in d.items() if k not in keys}


def document(energy=ENERGY, mode="homogeneous_equilibrium"):
    return {**G1, "two_phase_model": {**G1["two_phase_model"], "mode": mode}, "energy_model": energy}


def load(doc):
    return gribs.config.loader.configuration_from_document(doc).config


def with_settings(c, settings):
    """The configuration c with other energy settings (a Config subclass that carries them as a class attribute)."""
    sub = type("ConfigWithEnergy", (Config,), {"energy": settings})
    return sub(**asdict(c))


def with_changes(c, **changes):
    """The configuration c with some Config fields changed and the same energy settings (no assignment: the settings
    are the class attribute of a subclass, as ``c.energy`` is an attribute of the configuration the loader returns)."""
    sub = type("ConfigWithEnergy", (Config,), {"energy": c.energy})
    return sub(**{**asdict(c), **changes})


LEGACY_KEYS = {f.name for f in fields(Config)}


@pytest.mark.fast
def test_c01_without_the_section_there_is_no_energy_attribute_field_or_summary_key():
    c = load(G1)
    assert getattr(c, "energy", None) is None
    assert "energy" not in LEGACY_KEYS
    assert "energy" not in asdict(c) and set(asdict(c)) == LEGACY_KEYS      # resolved_configuration of a legacy run


@pytest.mark.fast
def test_c02_with_the_section_the_settings_are_an_attribute_and_asdict_is_unchanged():
    configuration = gribs.config.loader.configuration_from_document(document())
    c = configuration.config
    assert isinstance(c.energy, EnergyModeSettings)
    assert c.energy.initial_fill.temperature_K == 300.0 and c.energy.burn_surface.solid_density_kg_m3 == 1800.0
    assert c.energy.condensed.volume == "external_density" and len(c.energy.condensed.phases) == 2
    assert c.energy.condensed.transitions[0].melting_temperature_K == 2327.0
    assert "energy" not in asdict(c) and set(asdict(c)) == LEGACY_KEYS == set(asdict(load(G1)))
    assert configuration.document["energy_model"] == ENERGY             # the document keeps the section


@pytest.mark.fast
def test_c03_unknown_keys_and_choices_are_rejected():
    for energy, text in (
            ({**ENERGY, "bogus": 1}, "unknown key"),
            ({**ENERGY, "wall": {**WALL, "bogus": 1}}, "unknown key"),
            ({**ENERGY, "condensed": {**CONDENSED, "volume": "neglected"}}, "not supported"),
            ({**ENERGY, "burn_surface": {**SURFACE, "form": "b"}}, "not supported")):
        with pytest.raises(ConfigurationError, match=text):
            load(document(energy))


@pytest.mark.fast
def test_c04_eta_t0_other_than_one_is_rejected_in_energy_mode_only():
    c = load(document())
    gribs.config.validate(c)                                            # eta_T0 = 1: accepted
    with pytest.raises(ConfigurationError, match="eta_T0"):
        gribs.config.validate(with_changes(c, eta_T0=0.9))
    gribs.config.validate(replace(load(G1), eta_T0=0.9))                # legacy: unchanged


@pytest.mark.fast
def test_c05_a_wall_that_is_not_adiabatic_needs_a_heat_rate():
    with pytest.raises(ConfigurationError, match="heat_rate_W"):
        load(document({**ENERGY, "wall": {"adiabatic": False}}))
    c = load(document({**ENERGY, "wall": {"adiabatic": False, "heat_rate_W": 150.0}}))
    assert c.energy.wall.heat_rate_W == 150.0
    gribs.config.validate(c)


@pytest.mark.fast
def test_c06_an_explicit_adiabatic_wall_needs_no_heat_rate():
    c = load(document())
    assert c.energy.wall.adiabatic is True and c.energy.wall.heat_rate_W is None
    gribs.config.validate(c)
    both = replace(c.energy, wall=WallHeat(adiabatic=True, heat_rate_W=5.0))
    with pytest.raises(ConfigurationError, match="must not be given"):
        gribs.config.validate(with_settings(c, both))


@pytest.mark.fast
def test_c07_the_equilibrium_pressure_diagnostic_is_rejected_in_energy_mode_only():
    with pytest.raises(ConfigurationError, match="equilibrium_pressure_diagnostic"):
        gribs.config.validate(with_changes(load(document()), equilibrium_pressure_diagnostic="sampled"))
    gribs.config.validate(replace(load(G1), equilibrium_pressure_diagnostic="sampled"))   # legacy: accepted


@pytest.mark.fast
def test_c08_an_initial_state_with_p0_only_is_rejected():
    with pytest.raises(ConfigurationError, match="temperature_K"):
        load(document({**ENERGY, "initial_fill": without(FILL, "temperature_K")}))


@pytest.mark.fast
def test_c09_the_grain_far_temperature_is_required():
    with pytest.raises(ConfigurationError, match="grain_far_temperature_K"):
        load(document({**ENERGY, "burn_surface": without(SURFACE, "grain_far_temperature_K")}))


@pytest.mark.fast
def test_c10_the_solid_density_and_enthalpy_are_required():
    for missing in ("solid_density_kg_m3", "solid_enthalpy"):
        with pytest.raises(ConfigurationError, match=missing):
            load(document({**ENERGY, "burn_surface": without(SURFACE, missing)}))


@pytest.mark.fast
def test_c11_the_three_recovery_tolerances_are_required_without_defaults():
    for key in RECOVERY:
        with pytest.raises(ConfigurationError, match=key):
            load(document({**ENERGY, "recovery": without(RECOVERY, key)}))
    with pytest.raises(ConfigurationError, match="recovery"):
        load(document(without(ENERGY, "recovery")))


@pytest.mark.fast
def test_c12_an_initial_fill_other_than_the_products_is_rejected():
    with pytest.raises(ConfigurationError, match="composition"):
        load(document({**ENERGY, "initial_fill": {**FILL, "composition": "air"}}))


@pytest.mark.fast
def test_c13_a_two_phase_mode_other_than_homogeneous_equilibrium_is_rejected():
    with pytest.raises(ConfigurationError, match="homogeneous_equilibrium"):
        gribs.config.validate(load(document(mode="single_phase_legacy")))


@pytest.mark.fast
def test_c14_run_model_is_not_implemented_for_an_energy_configuration():
    with pytest.raises(NotImplementedError, match="energy mode is not implemented"):
        gribs.simulation.integrator.run_model(load(document()))
