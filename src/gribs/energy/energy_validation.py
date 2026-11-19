"""Cross-field checks of an energy-mode configuration (A-6 stage 2, C3; contract 164 sections 2b-4.12, 2b-FIX-2, 3-D).

Called by gribs.config.validate (lazily, and only for a configuration that has energy settings), so that a Config
built or changed directly in Python is checked as well as one read from JSON.  Imports only gribs.errors and
gribs.constants.
"""
import re

from gribs.constants import TWO_PHASE_HOMOGENEOUS
from gribs.energy.config_model import INITIAL_FILL_COMPOSITIONS
from gribs.errors import ConfigurationError

#: The error code of the rejection of an igniter in energy mode (a ConfigurationError, exit code 2 of the CLI).
IGNITER_REJECTED_CODE = "E_ENERGY_MODE_IGNITER_REJECTED"
#: The error code of a phase table that does not match the declared phase transitions (a ConfigurationError, exit code 2).
UNDECLARED_PHASE_TRANSITION_CODE = "E_ENERGY_UNDECLARED_PHASE_TRANSITION"

#: A species name is split into its base and a phase mark in brackets: "AL2O3(a)" -> ("AL2O3", "a").  A name without a mark has
#: no mark.  Two names are the same species only when their bases are equal as strings (the case is significant: CO and Co differ).
_SPECIES_MARK = re.compile(r"^(?P<base>.+)\((?P<mark>[A-Za-z]{1,3})\)$")
_LIQUID_MARK = "L"
_PHASE_RULE = (
    "Rule: a species name is split with the regular expression ^(?P<base>.+)\\((?P<mark>[A-Za-z]{1,3})\\)$ ; names with the same "
    "base (compared case-sensitively) are the same species; the phase mark (L) means liquid and phase = liquid requires the mark "
    "(L); a solid has any other 1-3 letter mark (a, b, cr, gr, II, ...) or none; a solid and a liquid of the same species need "
    "exactly one transition (solid_species, liquid_species) in energy_model.condensed.transitions; two solids of one species (a "
    "solid-solid transition) are not supported.")


def _split_species(name: str):
    m = _SPECIES_MARK.fullmatch(name)
    return (m.group("base"), m.group("mark")) if m else (name, None)


def _reject_phase_table(reason: str) -> None:
    raise ConfigurationError(f"{UNDECLARED_PHASE_TRANSITION_CODE}: energy_model.condensed: {reason}. {_PHASE_RULE}")


def _check_phase_transitions(condensed) -> None:
    """The phase table and the declared transitions must agree (rules R-1..R-6 of PRE_REGISTRATION_C15; decided from the
    configuration alone, the provider is not called)."""
    names = [p.species for p in condensed.phases]
    for name in sorted(set(names)):
        if names.count(name) > 1:
            _reject_phase_table(f"the species {name!r} is listed {names.count(name)} times in phases")
    groups = {}
    for p in condensed.phases:
        base, mark = _split_species(p.species)
        if (mark == _LIQUID_MARK) != (p.phase == "liquid"):
            _reject_phase_table(f"the phase {p.phase!r} of {p.species!r} does not match its phase mark {mark!r} "
                                f"(the mark (L) and phase = liquid go together)")
        groups.setdefault(base, []).append(p)
    pairs = [(t.solid_species, t.liquid_species) for t in condensed.transitions]
    needed = set()
    for base, members in groups.items():
        if len(members) == 1:
            continue
        solids = [p.species for p in members if p.phase == "solid"]
        liquids = [p.species for p in members if p.phase == "liquid"]
        if len(members) != 2 or len(solids) != 1 or len(liquids) != 1:
            _reject_phase_table(f"the species {base!r} has the phases {[p.species for p in members]} (more than one solid, more "
                                f"than one liquid, or more than two phases)")
        pair = (solids[0], liquids[0])
        if pairs.count(pair) != 1:
            _reject_phase_table(f"the solid {pair[0]!r} and the liquid {pair[1]!r} of one species need exactly one declared "
                                f"transition, found {pairs.count(pair)}")
        needed.add(pair)
    for pair in pairs:
        if pair not in needed:
            _reject_phase_table(f"the declared transition {pair[0]!r} -> {pair[1]!r} is not a solid and a liquid of one species "
                                f"of the phases table")


def validate_energy_mode(c) -> None:
    e = c.energy
    if c.eta_T0 != 1.0:
        raise ConfigurationError(
            f"energy_model: temperature_efficiency (eta_T0) must be 1 in energy mode, got {c.eta_T0:g}. The efficiency "
            f"multiplies T only and does not enter the energy equation; give the heat loss as energy_model.wall.")
    if c.equilibrium_pressure_diagnostic != "off":
        raise ConfigurationError(
            f"energy_model: equilibrium_pressure_diagnostic = {c.equilibrium_pressure_diagnostic!r} is not provided in "
            f"energy mode (the diagnostic evaluates the burn rate and the nozzle from p alone, which is inconsistent "
            f"with the energy-mode state T); use 'off'.")
    if c.two_phase_mode != TWO_PHASE_HOMOGENEOUS:
        raise ConfigurationError(
            f"energy_model: two_phase_model.mode must be {TWO_PHASE_HOMOGENEOUS!r} in energy mode, got "
            f"{c.two_phase_mode!r}.")
    if e.wall.adiabatic is False and e.wall.heat_rate_W is None:
        raise ConfigurationError(
            "energy_model.wall: heat_rate_W is required unless adiabatic is declared true (no default heat loss).")
    if e.wall.adiabatic is True and e.wall.heat_rate_W is not None:
        raise ConfigurationError("energy_model.wall: heat_rate_W must not be given when adiabatic is true.")
    if c.ign_mdot > 0.0 and c.ign_time > 0.0:
        raise ConfigurationError(
            f"{IGNITER_REJECTED_CODE}: an igniter (igniter.mass_flow_kg_s = {c.ign_mdot:g}, igniter.duration_s = "
            f"{c.ign_time:g}) is not supported in energy mode: the stagnation enthalpy of the igniter flow is not an input of "
            f"energy_model. Remove the igniter (a zero mass flow or a zero duration) or remove the energy_model section.")
    if e.initial_fill.composition not in INITIAL_FILL_COMPOSITIONS:
        raise ConfigurationError(
            f"energy_model.initial_fill.composition: {e.initial_fill.composition!r} is not supported; the initial fill "
            f"is limited to {list(INITIAL_FILL_COMPOSITIONS)}.")
    _check_phase_transitions(e.condensed)
