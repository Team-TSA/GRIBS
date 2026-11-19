"""Settings of the optional energy_model section (A-6 stage 2, C3; contract: report 164).

Frozen records only; no import of any gribs module (standard library only), so every layer may import them.
The parsed settings are attached to the run configuration as the attribute ``c.energy`` (like ``reactants_spec_list``),
not as a field of ``Config``: dataclasses.asdict(cfg), which the CLI and the API write as ``resolved_configuration``,
returns the fields only, so a legacy summary never shows an "energy" key (tests/test_energy_config_schema.py).

The solid density rho_p of the burning mass flow is NOT a setting of this section (C11, B-15): there is one rho_p, Config.rho_p.

The values the tests use for these inputs are assumed inputs, not data (PRE_REGISTRATION_C3.md, held items H-01..H-08).
"""
from dataclasses import dataclass
from typing import Mapping, Optional, Tuple

#: The only initial fill energy mode accepts: the equilibrium products of the main propellant (contract 2b-FIX-2).
INITIAL_FILL_COMPOSITIONS = ("main_propellant_products",)
#: The burn-surface form a configuration may select: (a) mdot*h_s - p_s*A_b*r (contract section 2 #3).
BURN_SURFACE_FORMS = ("a",)
#: The only condensed-volume treatment of energy mode (D-07): an external density per condensed phase.
CONDENSED_VOLUME_CHOICES = ("external_density",)
CONDENSED_PHASES = ("solid", "liquid")


@dataclass(frozen=True)
class InitialFill:
    temperature_K: float            # T_fill,0 of the initial charge at t = 0
    composition: str                # one of INITIAL_FILL_COMPOSITIONS
    element_tolerance: float        # tolerance of the t = 0 element check (held H-04; no default)


@dataclass(frozen=True)
class WallHeat:
    adiabatic: bool                 # explicit declaration; the default is not adiabatic
    heat_rate_W: Optional[float]    # Qdot_wall, heat into the chamber [W]; required unless adiabatic


@dataclass(frozen=True)
class SolidEnthalpy:
    reference_J_kg: float           # h_s at reference_temperature_K, on the absolute basis of the products (formation enthalpy included)
    reference_temperature_K: float
    cp_J_kgK: float                 # constant specific heat: h_s(T) = reference + cp (T - reference_temperature)


@dataclass(frozen=True)
class BurnSurface:
    form: str                       # one of BURN_SURFACE_FORMS
    grain_far_temperature_K: float  # T_grain,inf: where h_s is evaluated (contract 2b-4.3)
    solid_enthalpy: SolidEnthalpy


@dataclass(frozen=True)
class CondensedPhaseSpec:
    species: str
    phase: str                      # one of CONDENSED_PHASES
    density_kg_m3: float            # constant rho_c (D-10); required for every phase (it gives the volume)
    #: The three values of h = enthalpy_reference + cp (T - reference_temperature_K) (D-10) are required for the two phases of a
    #: declared transition and rejected for any other phase (an input that no number uses; C16, energy_validation).
    cp_J_kgK: Optional[float]
    enthalpy_reference_J_kg: Optional[float]   # h at reference_temperature_K, same basis as the products
    reference_temperature_K: Optional[float]


@dataclass(frozen=True)
class PhaseTransitionSpec:
    solid_species: str
    liquid_species: str
    melting_temperature_K: float    # T_m


@dataclass(frozen=True)
class CondensedData:
    volume: str                     # one of CONDENSED_VOLUME_CHOICES (key energy_model.condensed.volume, D-07)
    phases: Tuple[CondensedPhaseSpec, ...]
    transitions: Tuple[PhaseTransitionSpec, ...]
    #: energy_model.condensed.latent_heat_check_tolerance_rel: required (no default) when a transition is declared, rejected when
    #: none is.  The t = 0 check |L_calc - L_input| <= tol max(|L_input|, |L_calc|) uses it (C16).
    latent_heat_check_tolerance_rel: Optional[float] = None


@dataclass(frozen=True)
class RecoverySettings:
    phase_fraction_tolerance: float
    root_tolerance: float
    plateau_edge_tolerance: float


@dataclass(frozen=True)
class RecoverySearch:
    """Where the state recovery looks (C11; keys of energy_model.recovery next to the three tolerances)."""
    temperature_search_K: Tuple[float, float]   # searched temperature range (T_lo, T_hi)
    start_temperatures_K: Tuple[float, ...]     # start temperatures; every start must give the same state


@dataclass(frozen=True)
class EnergyModeSettings:
    initial_fill: InitialFill
    wall: WallHeat
    burn_surface: BurnSurface
    condensed: CondensedData
    recovery: RecoverySettings
    recovery_search: RecoverySearch
    #: Explicit formulas of the reactants, reactant name -> {element symbol: atoms} (C9; optional, None when not given).
    #: Never inferred from species names; energy mode will require it at run time (the wiring, C11).
    propellant_formulas: Optional[Mapping[str, Mapping[str, float]]] = None
