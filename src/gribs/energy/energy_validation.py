"""Cross-field checks of an energy-mode configuration (A-6 stage 2, C3; contract 164 sections 2b-4.12, 2b-FIX-2, 3-D).

Called by gribs.config.validate (lazily, and only for a configuration that has energy settings), so that a Config
built or changed directly in Python is checked as well as one read from JSON.  Imports only gribs.errors and
gribs.constants.
"""
from gribs.constants import TWO_PHASE_HOMOGENEOUS
from gribs.energy.config_model import INITIAL_FILL_COMPOSITIONS
from gribs.errors import ConfigurationError


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
    if e.initial_fill.composition not in INITIAL_FILL_COMPOSITIONS:
        raise ConfigurationError(
            f"energy_model.initial_fill.composition: {e.initial_fill.composition!r} is not supported; the initial fill "
            f"is limited to {list(INITIAL_FILL_COMPOSITIONS)}.")
