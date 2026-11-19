"""Parser of the optional energy_model section (A-6 stage 2, C3). Reads a JSON object; every key is required (no
defaults) except propellant_formulas (C9, optional), and unknown keys are rejected, with the messages of the other sections (ConfigurationError).

Imported lazily by gribs.config.loader.configuration_from_document, only when the section is present.
"""
from __future__ import annotations

from typing import Any, List

from gribs.config.loader import (
    _boolean, _choice, _get, _mapping, _number, _reject_unknown, _section, _string,
)
from gribs.energy.config_model import (
    BURN_SURFACE_FORMS, CONDENSED_PHASES, CONDENSED_VOLUME_CHOICES, INITIAL_FILL_COMPOSITIONS,
    BurnSurface, CondensedData, CondensedPhaseSpec, EnergyModeSettings, InitialFill, PhaseTransitionSpec,
    RecoverySettings, SolidEnthalpy, WallHeat,
)
from gribs.energy.elements import check_formulas
from gribs.errors import ConfigurationError

_KEYS = ("initial_fill", "wall", "burn_surface", "condensed", "recovery", "propellant_formulas")
_FILL_KEYS = ("temperature_K", "composition", "element_tolerance")
_WALL_KEYS = ("adiabatic", "heat_rate_W")
_SURFACE_KEYS = ("form", "grain_far_temperature_K", "solid_density_kg_m3", "solid_enthalpy")
_SOLID_H_KEYS = ("reference_J_kg", "reference_temperature_K", "cp_J_kgK")
_CONDENSED_KEYS = ("volume", "phases", "transitions")
_PHASE_KEYS = ("species", "phase", "density_kg_m3", "cp_J_kgK", "enthalpy_reference_J_kg", "reference_temperature_K")
_TRANSITION_KEYS = ("solid_species", "liquid_species", "melting_temperature_K")
_RECOVERY_KEYS = ("phase_fraction_tolerance", "root_tolerance", "plateau_edge_tolerance")


def _array_of_objects(node: dict, key: str, path: str) -> List[dict]:
    value = _get(node, key, path)
    if not isinstance(value, list):
        raise ConfigurationError(f"{path}.{key}: expected a JSON array, got {type(value).__name__}")
    return [_mapping(item, f"{path}.{key}[{i}]") for i, item in enumerate(value)]


def parse_energy_model(node: Any, path: str = "energy_model") -> EnergyModeSettings:
    root = _mapping(node, path)
    _reject_unknown(root, _KEYS, path)

    fill = _section(root, "initial_fill", path)
    _reject_unknown(fill, _FILL_KEYS, f"{path}.initial_fill")
    initial_fill = InitialFill(
        temperature_K=float(_number(fill, "temperature_K", f"{path}.initial_fill", exclusive_minimum=0.0)),
        composition=str(_choice(fill, "composition", f"{path}.initial_fill", INITIAL_FILL_COMPOSITIONS)),
        element_tolerance=float(_number(fill, "element_tolerance", f"{path}.initial_fill", exclusive_minimum=0.0)))

    wall = _section(root, "wall", path)
    _reject_unknown(wall, _WALL_KEYS, f"{path}.wall")
    adiabatic = bool(_boolean(wall, "adiabatic", f"{path}.wall"))
    heat = _number(wall, "heat_rate_W", f"{path}.wall", required=not adiabatic)
    wall_heat = WallHeat(adiabatic=adiabatic, heat_rate_W=None if heat is None else float(heat))

    surface = _section(root, "burn_surface", path)
    _reject_unknown(surface, _SURFACE_KEYS, f"{path}.burn_surface")
    solid = _mapping(_get(surface, "solid_enthalpy", f"{path}.burn_surface"), f"{path}.burn_surface.solid_enthalpy")
    _reject_unknown(solid, _SOLID_H_KEYS, f"{path}.burn_surface.solid_enthalpy")
    burn_surface = BurnSurface(
        form=str(_choice(surface, "form", f"{path}.burn_surface", BURN_SURFACE_FORMS)),
        grain_far_temperature_K=float(_number(surface, "grain_far_temperature_K", f"{path}.burn_surface",
                                              exclusive_minimum=0.0)),
        solid_density_kg_m3=float(_number(surface, "solid_density_kg_m3", f"{path}.burn_surface",
                                          exclusive_minimum=0.0)),
        solid_enthalpy=SolidEnthalpy(
            reference_J_kg=float(_number(solid, "reference_J_kg", f"{path}.burn_surface.solid_enthalpy")),
            reference_temperature_K=float(_number(solid, "reference_temperature_K",
                                                  f"{path}.burn_surface.solid_enthalpy", exclusive_minimum=0.0)),
            cp_J_kgK=float(_number(solid, "cp_J_kgK", f"{path}.burn_surface.solid_enthalpy",
                                   exclusive_minimum=0.0))))

    cond = _section(root, "condensed", path)
    _reject_unknown(cond, _CONDENSED_KEYS, f"{path}.condensed")
    phases = []
    for i, p in enumerate(_array_of_objects(cond, "phases", f"{path}.condensed")):
        pp = f"{path}.condensed.phases[{i}]"
        _reject_unknown(p, _PHASE_KEYS, pp)
        phases.append(CondensedPhaseSpec(
            species=str(_string(p, "species", pp, allow_empty=False)),
            phase=str(_choice(p, "phase", pp, CONDENSED_PHASES)),
            density_kg_m3=float(_number(p, "density_kg_m3", pp, exclusive_minimum=0.0)),
            cp_J_kgK=float(_number(p, "cp_J_kgK", pp, exclusive_minimum=0.0)),
            enthalpy_reference_J_kg=float(_number(p, "enthalpy_reference_J_kg", pp)),
            reference_temperature_K=float(_number(p, "reference_temperature_K", pp, exclusive_minimum=0.0))))
    transitions = []
    for i, t in enumerate(_array_of_objects(cond, "transitions", f"{path}.condensed")):
        tp = f"{path}.condensed.transitions[{i}]"
        _reject_unknown(t, _TRANSITION_KEYS, tp)
        transitions.append(PhaseTransitionSpec(
            solid_species=str(_string(t, "solid_species", tp, allow_empty=False)),
            liquid_species=str(_string(t, "liquid_species", tp, allow_empty=False)),
            melting_temperature_K=float(_number(t, "melting_temperature_K", tp, exclusive_minimum=0.0))))
    condensed = CondensedData(
        volume=str(_choice(cond, "volume", f"{path}.condensed", CONDENSED_VOLUME_CHOICES)),
        phases=tuple(phases), transitions=tuple(transitions))

    rec = _section(root, "recovery", path)
    _reject_unknown(rec, _RECOVERY_KEYS, f"{path}.recovery")
    recovery = RecoverySettings(
        phase_fraction_tolerance=float(_number(rec, "phase_fraction_tolerance", f"{path}.recovery",
                                               exclusive_minimum=0.0)),
        root_tolerance=float(_number(rec, "root_tolerance", f"{path}.recovery", exclusive_minimum=0.0)),
        plateau_edge_tolerance=float(_number(rec, "plateau_edge_tolerance", f"{path}.recovery",
                                             exclusive_minimum=0.0)))

    formulas = None
    if "propellant_formulas" in root:
        table = _mapping(root["propellant_formulas"], f"{path}.propellant_formulas")
        check_formulas(table, f"{path}.propellant_formulas")
        formulas = {str(name): {str(e): float(n) for e, n in formula.items()} for name, formula in table.items()}

    return EnergyModeSettings(initial_fill=initial_fill, wall=wall_heat, burn_surface=burn_surface,
                              condensed=condensed, recovery=recovery, propellant_formulas=formulas)
