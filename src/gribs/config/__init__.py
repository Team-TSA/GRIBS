"""Runtime configuration, derived geometry and configuration checks (A-6 unit 2). Moved verbatim from the
entry: section 1 (Config, PARAM_DOC), section 2 (geometry helpers) and, from section 7, _is_nonfinite_number,
require_finite_config and validate (tests patch validate here: "gribs.config.validate")."""
from __future__ import annotations

import dataclasses
import decimal
import math
from dataclasses import dataclass
from typing import Any, Tuple

import numpy as np

from gribs.constants import KNOWN_BACKENDS, TWO_PHASE_HOMOGENEOUS


# ---- section 1: internal runtime configuration (field defaults are constructor fallbacks, not user input) ----
@dataclass
class Config:
    # ---------------- propellant / combustion ----------------
    rho_p: float = 0.0            # runtime bulk density [kg/m^3], set from the
                                  # additive-volume mixture density and the
                                  # packing fraction (never from a gas density)
    a_burn: float = 0.7e-3        # burn-rate coefficient [m/s] at p = p_ref
    n_burn: float = 0.5           # pressure exponent [-]
    p_ref: float = 1.0e6          # reference pressure of the burn law [Pa]
    sigma_p: float = 0.0          # burn-rate temperature sensitivity [1/K]
    T_grain: float = 294.0        # initial grain temperature [K]
    T_ref: float = 294.0          # reference grain temperature [K]
    ero_alpha: float = 0.0        # Lenoir-Robert erosive coefficient (0 = off)
    ero_beta: float = 53.0        # Lenoir-Robert damping constant [-]

    # ---------------- grain geometry ----------------
    R_i0: float = 3.00e-3         # initial bore radius [m]
    R_p: float = 1.20e-2          # grain outer radius [m]
    L_p0: float = 0.080           # initial grain length [m]
    V_g0: float = 6.79e-6         # initial free (gas) volume [m^3]
    n_end: float = 1.0            # number of burning end faces (0, 1 or 2)

    # ---------------- nozzle ----------------
    R_t0: float = 1.50e-3         # initial throat radius [m]
    eps_nozzle: float = 1.0       # area ratio Ae/At (1.0 = converging nozzle)
    Cd: float = 1.0               # nozzle discharge coefficient [-]
    eta_thrust: float = 1.0       # thrust (momentum) efficiency [-]
    transition_policy: str = "jump"  # 'jump' | 'shock' (shock reserved)
    use_separation: bool = True   # model overexpanded flow separation
    sep_ratio: float = 0.4        # Summerfield separation criterion pe/pa [-]
    ero_throat_c: float = 0.0     # throat erosion rate [m/s] at p_ref
    ero_throat_m: float = 1.0     # throat erosion pressure exponent [-]

    # ---------------- environment / initial state ----------------
    p_a: float = 101325.0         # ambient (back) pressure [Pa]
    p0_init: float = 0.2e6        # initial chamber pressure [Pa, absolute]

    # ---------------- igniter (optional gas injection) ----------------
    ign_mdot: float = 0.0         # igniter mass flow [kg/s]
    ign_time: float = 0.0         # igniter duration [s]

    # ---------------- thermochemistry backend (mandatory, explicit) ----------------
    thermo_backend: str = ""      # must be set to one of KNOWN_BACKENDS
    cea_pressure_points: int = 81         # logarithmic table grid points
    cea_cache_enabled: bool = True
    cea_rebuild_cache: bool = False
    cea_cache_dir: str = ""               # empty -> results/thermo_cache

    # ---- reference_constant (Phase A-5 analytic reference) options ----
    reference_temperature_K: float = 3000.0
    reference_gamma: float = 1.2
    reference_molecular_weight_kg_kmol: float = 25.0

    # ---- cea_python (official NASA CEA Python package) options ----
    cea_py_ions: bool = False
    cea_py_transport: bool = False
    cea_py_trace: float = 1.0e-10
    cea_py_products_from_reactants: bool = True
    cea_py_product_species: Tuple[str, ...] = ()
    cea_py_omit_species: Tuple[str, ...] = ()
    cea_py_insert_species: Tuple[str, ...] = ()
    cea_py_molecular_weight: str = "gas_phase_M"   # frozen compatibility token; selects EqSolution.M = M_CEA (kg-total/kmol-gas); 'total_MW' rejected
    cea_py_smooth_truncation: bool = False
    cea_py_truncation_width: float = -1.0
    cea_py_rocket_diagnostics: bool = True

    # ---- cea_legacy_executable (transitional compatibility) options ----
    cea_legacy_executable_path: str = "fcea2"
    cea_legacy_data_directory: str = ""
    cea_legacy_timeout_s: float = 120.0
    cea_legacy_batch_size: int = 8
    cea_legacy_trace: float = 1.0e-10

    # ---------------- chamber-property table / validity ----------------
    p_fit_min: float = 1.0e5      # lower validity limit of the property table [Pa]
    p_fit_max: float = 8.0e6      # upper validity limit of the property table [Pa]
    eta_T0: float = 1.0           # combustion/heat-loss efficiency on T0 [-]
    property_policy: str = "clamp"  # 'clamp' | 'extrapolate' (outside range)

    # ---------------- two-phase chamber / nozzle model ----------------
    # 'homogeneous_equilibrium' (standard): the conserved chamber mass is the
    #   TOTAL product mass mt = mg + mc; EOS p*Vg = mt*Yg(p)*Rg(p)*T0(p); the
    #   nozzle total flow comes directly from the isentropic formula
    #   evaluated with R_eff = Yg(p)*Rg(p), with no additional /Yg factor
    #   (complete entrainment, no particle slip).
    # 'single_phase_legacy': exact GRIBS v0.4.3-alpha behaviour (Yg := 1),
    #   kept ONLY for regression comparison.
    two_phase_mode: str = TWO_PHASE_HOMOGENEOUS
    # Condensed-phase volume treatment; only 'neglected' is implemented.
    two_phase_condensed_volume: str = "neglected"
    # Nozzle two-phase treatment; only 'complete' (entrainment) is implemented.
    two_phase_nozzle_entrainment: str = "complete"

    # ---------------- nozzle low-pressure policy ----------------
    unchoked_policy: str = "switch"  # 'switch' (subsonic branch) | 'stop'

    # ---------------- solver ----------------
    method: str = "LSODA"         # LSODA | BDF | Radau
    rtol: float = 1.0e-10
    atol_p: float = 1.0e-3        # [Pa]
    atol_x: float = 1.0e-13       # [m]
    t_max: float = 600.0          # integration horizon of the burning phase [s]
    max_step_burn: float = 0.02   # [s]
    blowdown: bool = True         # integrate the post-burnout blowdown
    blowdown_tmax: float = 5.0    # [s]
    max_step_blow: float = 5.0e-4  # [s]

    # ---------------- diagnostics ----------------
    # Equilibrium pressure diagnostic (p_eq): 'off' disables it (fastest),
    # 'sampled' evaluates up to 200 points during the burning phase.
    # Default is 'off' because the diagnostic is expensive (each point
    # calls equilibrium_pressure() which does a pressure sweep + brentq).
    equilibrium_pressure_diagnostic: str = "off"  # 'off' | 'sampled'


PARAM_DOC = {
    "rho_p": "bulk propellant density [kg/m^3] (computed from the reactants JSON: additive-volume ideal density x packing fraction)",
    "a_burn": "burn-rate coefficient [m/s at p_ref]",
    "n_burn": "burn-rate pressure exponent [-]",
    "p_ref": "reference pressure of the burn law [Pa]",
    "sigma_p": "burn-rate temperature sensitivity [1/K]",
    "T_grain": "initial grain temperature [K]",
    "T_ref": "reference grain temperature [K]",
    "ero_alpha": "Lenoir-Robert erosive-burning coefficient (0 disables it)",
    "ero_beta": "Lenoir-Robert damping constant [-]",
    "R_i0": "initial bore radius [m]",
    "R_p": "grain outer radius [m]",
    "L_p0": "initial grain length [m]",
    "V_g0": "initial free gas volume [m^3]",
    "n_end": "number of burning end faces (0, 1 or 2)",
    "R_t0": "initial throat radius [m]",
    "eps_nozzle": "nozzle area ratio Ae/At [-]",
    "Cd": "nozzle discharge coefficient [-]",
    "eta_thrust": "thrust (momentum) efficiency [-]",
    "use_separation": "model overexpanded flow separation (true/false)",
    "sep_ratio": "separation criterion pe/pa [-]",
    "ero_throat_c": "throat erosion rate at p_ref [m/s]",
    "ero_throat_m": "throat erosion pressure exponent [-]",
    "p_a": "ambient back pressure [Pa]",
    "p0_init": "initial chamber pressure [Pa abs]",
    "ign_mdot": "igniter mass flow [kg/s]",
    "ign_time": "igniter duration [s]",
    "thermo_backend": ("mandatory explicit choice: 'cea_python' (official NASA CEA "
                       "Python package) | 'cea_legacy_executable' (external fcea2) | "
                       "'reference_constant' (Phase A-5 analytic reference)"),
    "reference_temperature_K": "reference_constant backend temperature [K]",
    "reference_gamma": "reference_constant backend specific-heat ratio [-]",
    "reference_molecular_weight_kg_kmol": (
        "reference_constant backend molecular weight [kg/kmol]"),
    "cea_pressure_points": "number of logarithmic chamber-pressure grid points",
    "cea_cache_enabled": "reuse the reproducible chamber-property cache",
    "cea_rebuild_cache": "ignore and rebuild an existing chamber-property cache",
    "cea_cache_dir": "directory for reproducible chamber-property caches",
    "cea_py_ions": "official CEA backend: include ionized species",
    "cea_py_transport": "official CEA backend: compute transport properties",
    "cea_py_trace": "official CEA backend: trace-species threshold (<=0 uses the API default)",
    "cea_py_products_from_reactants": "official CEA backend: build the product set from the reactant elements",
    "cea_py_molecular_weight": ("'gas_phase_M' only (frozen compatibility token name). "
                                "It selects EqSolution.M = M_CEA = 1/EqSolution.n "
                                "[kg-total/kmol-gas] = M_gas/Yg, i.e. the total-mixture "
                                "mass per kmol of gas; it is NOT the gas-phase mean "
                                "molecular weight. The former 'total_MW' choice is "
                                "rejected (v0.4.4): EqSolution.MW counts condensed moles, "
                                "is tabulated as a diagnostic and never enters the EOS "
                                "or nozzle."),
    "cea_legacy_executable_path": "path to the native CEA executable (compatibility backend)",
    "cea_legacy_data_directory": "directory containing thermo.lib and trans.lib",
    "cea_legacy_timeout_s": "timeout for each external CEA batch [s]",
    "cea_legacy_batch_size": "assigned pressures per external CEA run (8 = tested legacy workaround)",
    "cea_legacy_trace": "external CEA output trace threshold",
    "p_fit_min": "lower validity limit of the property table [Pa]",
    "p_fit_max": "upper validity limit of the property table [Pa]",
    "eta_T0": "combustion efficiency applied to T0 [-]",
    "property_policy": "'clamp' or 'extrapolate' outside the table range",
    "two_phase_mode": ("'homogeneous_equilibrium' (total product mass conserved, "
                       "EOS on the gas phase via Psi = R_eff*T0) | "
                       "'single_phase_legacy' (v0.4.3 behaviour, regression only)"),
    "two_phase_condensed_volume": "condensed-phase volume treatment ('neglected')",
    "two_phase_nozzle_entrainment": ("nozzle two-phase treatment ('complete' = both "
                                     "phases leave with the same velocity)"),
    "unchoked_policy": "'switch' to subsonic branch, or 'stop'",
    "method": "ODE method: LSODA | BDF | Radau",
    "rtol": "relative tolerance",
    "atol_p": "absolute tolerance on pressure [Pa]",
    "atol_x": "absolute tolerance on burn depth [m]",
    "t_max": "integration horizon of the burning phase [s]",
    "max_step_burn": "maximum solver step, burning phase [s]",
    "blowdown": "integrate post-burnout blowdown (true/false)",
    "blowdown_tmax": "blowdown horizon [s]",
    "max_step_blow": "maximum solver step, blowdown [s]",
}



# ---- section 2: derived geometry helpers ----
def web_thickness(c: Config) -> float:
    """Burn depth at which the grain is consumed (first geometric limit)."""
    radial = c.R_p - c.R_i0
    axial = c.L_p0 / c.n_end if c.n_end > 0.0 else math.inf
    return min(radial, axial)


def initial_propellant_volume(c: Config) -> float:
    return math.pi * (c.R_p ** 2 - c.R_i0 ** 2) * c.L_p0


def propellant_mass(c: Config) -> float:
    return c.rho_p * initial_propellant_volume(c)


def geometry(x: float, c: Config):
    """Return (Ri, Lp, Ab, Vg) for burn depth x [m] (clamped to the web)."""
    xw = web_thickness(c)
    xe = min(max(x, 0.0), xw)
    Ri = min(c.R_i0 + xe, c.R_p)
    Lp = max(c.L_p0 - c.n_end * xe, 0.0)
    dA = max(c.R_p ** 2 - Ri ** 2, 0.0)
    Ab = 2.0 * math.pi * Ri * Lp + c.n_end * math.pi * dA
    Vprop = math.pi * dA * Lp
    Vg = c.V_g0 + initial_propellant_volume(c) - Vprop
    return Ri, Lp, Ab, Vg


# ---- configuration checks (formerly at the end of section 7) ----
def _is_nonfinite_number(value: Any) -> bool:
    """True only for a numeric NaN or infinity; finite values are never converted.

    Python float (and subclasses such as numpy.float64) uses math.isfinite;
    other numpy floating scalars and 0-d float arrays use numpy.isfinite in
    their own precision (a finite long double beyond the float range stays
    finite); decimal.Decimal uses is_finite().  bool and int (including huge
    integers) are always finite.  Other types are left to the existing checks.
    """
    if isinstance(value, (bool, int)):
        return False
    if isinstance(value, float):
        return not math.isfinite(value)
    if isinstance(value, decimal.Decimal):
        return not value.is_finite()
    if isinstance(value, (np.floating, np.complexfloating)):
        return not bool(np.isfinite(value))
    if isinstance(value, np.ndarray) and value.ndim == 0 and value.dtype.kind in "fc":
        return not bool(np.isfinite(value))
    return False


def require_finite_config(c: Config) -> None:
    """Reject a NaN or infinite value in any Config field (finiteness only)."""
    for field in dataclasses.fields(c):
        value = getattr(c, field.name)
        if _is_nonfinite_number(value):
            raise ValueError(f"{field.name} must be a finite number, got {value!r}.")


def validate(c: Config) -> None:
    # Every numeric Config field must be finite, also for Config objects built
    # or changed directly in Python (run_model() calls validate() first).
    require_finite_config(c)
    if c.R_i0 >= c.R_p:
        raise ValueError("R_i0 must be smaller than R_p.")
    for name in ("rho_p", "a_burn", "p_ref", "R_p", "L_p0", "V_g0", "R_t0",
                 "p_a", "p0_init", "T_grain", "T_ref"):
        if getattr(c, name) <= 0.0:
            raise ValueError(f"{name} must be positive.")
    if c.eps_nozzle < 1.0:
        raise ValueError("eps_nozzle (Ae/At) must be >= 1.")
    if c.n_end < 0.0 or c.n_end > 2.0:
        raise ValueError("n_end must lie between 0 and 2.")
    if c.property_policy not in ("clamp", "extrapolate"):
        raise ValueError("property_policy must be 'clamp' or 'extrapolate'.")
    if c.thermo_backend not in KNOWN_BACKENDS:
        raise ValueError("thermo_backend must be one of: " + ", ".join(KNOWN_BACKENDS) + ".")
    if c.p_fit_min <= 0.0 or c.p_fit_max <= c.p_fit_min:
        raise ValueError("Require 0 < p_fit_min < p_fit_max.")
    if c.unchoked_policy not in ("switch", "stop"):
        raise ValueError("unchoked_policy must be 'switch' or 'stop'.")
    if c.transition_policy not in ("jump", "shock"):
        raise ValueError("transition_policy must be 'jump' or 'shock'.")
    # A-7: eta_T0 range check (review §3-E1)
    # eta_T0 > 1 would produce temperatures above the adiabatic flame
    # temperature, which is physically inconsistent with the current model.
    if not (0.0 < c.eta_T0 <= 1.5):
        raise ValueError(
            f"temperature_efficiency (eta_T0) must lie in (0, 1.5], got "
            f"{c.eta_T0}.  Values above 1.0 are physically inconsistent "
            f"with the current model (they exceed the adiabatic flame "
            f"temperature).  Use eta_T0 = 1.0 for the standard model.")
    # A-7: pressure range vs initial/ambient pressure (review §3-E2)
    # The blowdown terminal pressure approaches p_a; if p_a is too close to
    # p_fit_min, the tail of the blowdown is dominated by clamping.
    # This is a warning, not an error, because the default config works
    # adequately with the current margin.
    if c.p_a < c.p_fit_min:
        raise ValueError(
            f"ambient pressure p_a = {c.p_a:.0f} Pa is below the "
            f"property table lower limit p_fit_min = {c.p_fit_min:.0f} Pa.  "
            f"Require p_a >= p_fit_min.")
