#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GRIBS v0.4.3-alpha
================================================================================
JSON-CONFIGURED INTERNAL BALLISTICS OF A SOLID ROCKET MOTOR
  Grain    : cylindrical bore + N burning end face(s) (outer surface inhibited)
  Nozzle   : converging (Ae = At) or converging-diverging (Ae/At > 1)
  Chamber  : unsteady mass/state equation with pressure-dependent gas properties

--------------------------------------------------------------------------------
GOVERNING MODEL  (unchanged since v0.3.1-alpha)
--------------------------------------------------------------------------------
  Equation of state in the chamber        :  p0 * Vg = m_g * Theta(p0),
                                             Theta = R(p0) * T0(p0)
  Differentiating and substituting
      dVg/dt = Ab * r ,  dm_g/dt = mdot_gen - mdot_out
  gives the pressure equation used here (valid for any nozzle regime):

      dp0/dt = [ Theta*(mdot_gen - mdot_out) - p0*Ab*r ] / [ Vg*(1 - p0*Theta'/Theta) ]

  Geometry (exactly consistent: dVg/dx = Ab is verified numerically)
      Ri = Ri0 + x ,  Lp = Lp0 - N_end*x
      Ab = 2*pi*Ri*Lp + N_end*pi*(Rp^2 - Ri^2)
      Vg = Vg0 + Vp0 - pi*(Rp^2 - Ri^2)*Lp

  Burn rate (Saint-Robert + temperature sensitivity + optional erosive burning,
  Lenoir-Robert form, solved implicitly):
      r0 = a*(p0/p_ref)^n * exp(sigma_p*(T_grain - T_ref))
      r  = r0 + alpha*G^0.8/D_h^0.2 * exp(-beta*rho_p*r/G)

  Nozzle (isentropic, discharge coefficient Cd, thrust efficiency eta_F)
      choked      :  mdot = Cd*At*p0*sqrt(g/(R*T0))*(2/(g+1))^((g+1)/(2(g-1)))
                     exit Mach from the area-Mach relation (Ae/At)
      unchoked    :  pe = pa, subsonic exit Mach from the pressure ratio
                     (this branch is *continuous* with the choked branch)
      separated   :  Summerfield criterion pe < f_sep*pa (only for Ae/At > 1)
      thrust      :  F = eta_F*mdot*ve + (pe - pa)*Ae_effective

  Throat erosion :  dRt/dt = C_ero*(p0/p_ref)^m_ero  (Ae is held fixed)

  Solid propellant density (independent of any gas property):
      rho_ideal = 1 / sum_i(w_i/rho_i)        (additive-volume mixing)
      rho_p     = rho_ideal * packing_fraction
  The CEA combustion-product (gas) density is NEVER used as rho_p.

--------------------------------------------------------------------------------
THERMOCHEMISTRY BACKENDS  (v0.4)
--------------------------------------------------------------------------------
  The chamber-gas properties R(p0), T0(p0) and gamma_s(p0) are produced behind a
  single backend interface; the internal-ballistics solver never needs to know
  how they were generated.

    cea_python              (preferred)  official NASA CEA Python package
                                         ``import cea``  (https://github.com/nasa/CEA)
                                         HP equilibrium, no external executable,
                                         no thermo.lib/trans.lib beside GRIBS,
                                         no temporary .inp/.out/.plt files.
    cea_legacy_executable   (transitional compatibility) the v0.3.1-alpha pathway
                                         that runs the external fcea2 executable and
                                         parses its .plt output.  Kept for regression
                                         comparison only.

  The backend is an explicitly configured choice - there is no default and no
  silent fallback, and the v0.3 manual-correlation backend (legacy_fit) has been
  removed in v0.4.0-alpha (see CHANGELOG.md).

  Both backends build a pressure-indexed chamber-property table, cache it as
  reproducible JSON, and interpolate it with PCHIP in ln(p) - including an
  analytic derivative of Theta = R*T used by the ODE pressure equation.

--------------------------------------------------------------------------------
NUMERICS
--------------------------------------------------------------------------------
  State y = [p0, x, Rt, m_out, Impulse, m_gen] integrated with solve_ivp
  (LSODA/BDF/Radau, dense output, tight tolerances).  Quadratures for impulse
  and integrated masses are carried as ODE states, so they inherit the solver
  error control instead of relying on post-hoc trapezoidal sums.
  Root-finding events capture: burnout, choking loss/recovery, ambient pressure.
  A set of numerical self-tests is executed before every production run.

--------------------------------------------------------------------------------
OUTPUT
--------------------------------------------------------------------------------
  results/ballistics.png   single figure, 8 panels, publication-style
  results/time_history.csv full time history
  results/summary.txt      human-readable summary (English)
  results/summary.json     machine-readable summary + the exact input set +
                           full thermochemistry provenance
  results/thermo_cache/    reproducible chamber-property tables

  All user inputs are read from gribs_config.json beside this script.
  Select thermochemistry.backend in that file.  An alternative complete
  configuration may be supplied with:
      python GRIBS_v0.4.3-alpha.py --config another_config.json
  Run numerical self-tests only with:
      python GRIBS_v0.4.3-alpha.py --selftest

--------------------------------------------------------------------------------
STATUS / LIMITATIONS
--------------------------------------------------------------------------------
  * GRIBS keeps its own unsteady internal-ballistics and nozzle-flow model.
    CEA theoretical rocket performance (c*, Cf, Isp) is reported ONLY as a
    clearly-named diagnostic and never feeds the GRIBS results.
  * CEA does not provide burn-rate coefficients; the Saint-Robert inputs remain
    empirical user data.
  * Results require independent validation before any engineering use.
================================================================================
"""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import importlib
import importlib.metadata as importlib_metadata
import json
import math
import os
import platform
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.interpolate import PchipInterpolator

# ==============================================================================
# 0. PROGRAM IDENTITY AND CONSTANTS
# ==============================================================================
PROGRAM_NAME = "GRIBS"
PROGRAM_VERSION = "0.4.3-alpha"
SCHEMA_VERSION = "0.4.3-alpha"
#: Schema string of the pre-migration configuration; used only by the explicit
#: v0.3 -> v0.4 migration helper (never for normal operation).
LEGACY_SCHEMA_VERSION = "0.3.0-alpha"

BACKEND_CEA_PYTHON = "cea_python"
BACKEND_CEA_LEGACY_EXECUTABLE = "cea_legacy_executable"
#: The complete set of selectable backends.  ``thermochemistry.backend`` is a
#: required configuration key whose value must be one of these; the program never
#: picks a backend on the user's behalf.
KNOWN_BACKENDS = (BACKEND_CEA_PYTHON, BACKEND_CEA_LEGACY_EXECUTABLE)
#: Only used by the explicit v0.3 migration helper: the v0.3 name "cea2" meant
#: "the external fcea2 executable", which is now named cea_legacy_executable.
LEGACY_BACKEND_ALIASES = {"cea2": BACKEND_CEA_LEGACY_EXECUTABLE}
#: Backends that existed before v0.4.0-alpha and have been removed.  Naming one of
#: them in a configuration is an error with an actionable message - never a silent
#: substitution by another backend.
REMOVED_BACKENDS = {
    "legacy_fit": ("the manual R/T0/gamma correlation backend (legacy_fit) was "
                   "removed in GRIBS v0.4.0-alpha. Use 'cea_python' (official NASA "
                   "CEA Python package, preferred) or 'cea_legacy_executable' "
                   "(external fcea2 executable)."),
}

CEA_PROJECT_URL = "https://github.com/nasa/CEA"
CEA_INSTALL_COMMAND = "python -m pip install cea"
#: Official NASA CEA Python API series this file was developed and validated
#: against (see VALIDATION_REPORT.md).  Older majors do not provide the
#: EqSolver/EqSolution/RocketSolver surface used here.
CEA_MIN_API_VERSION = (3, 0)
CEA_TESTED_API_VERSION = (3, 3)

#: Universal gas constant used for R = Ru/M  [J/(kmol K)].
#: Verified at runtime against ``cea.R`` for the official backend.
R_UNIVERSAL = 8314.51
G0 = 9.80665                      # standard gravity [m/s^2]
_TRAPZ = np.trapezoid if hasattr(np, "trapezoid") else np.trapz

_PHASE_TAGS = ("(cr)", "(L)", "(l)", "(gr)", "(s)", "(a)", "(b)", "(I)", "(II)",
               "(III)", "(IV)", "(V)")


# ==============================================================================
# 0b. EXCEPTION HIERARCHY
# ==============================================================================
class GribsError(Exception):
    """Base class for all GRIBS exceptions."""


class ConfigurationError(GribsError):
    """JSON configuration is missing, malformed, or inconsistent."""


class ThermochemistryError(GribsError):
    """A thermochemistry backend could not produce physical gas properties."""


class ThermochemistryCacheError(ThermochemistryError):
    """A cached chamber-property table exists but is unusable."""


class BackendUnavailableError(ThermochemistryError):
    """A selected thermochemistry backend cannot be initialized."""


# ==============================================================================
# 1. INTERNAL RUNTIME CONFIGURATION
# Values are populated from the complete JSON configuration. Field defaults are
# implementation fallbacks required by the dataclass constructor, not normal user input.
# ==============================================================================
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

    # ---- cea_python (official NASA CEA Python package) options ----
    cea_py_ions: bool = False
    cea_py_transport: bool = False
    cea_py_trace: float = 1.0e-10
    cea_py_products_from_reactants: bool = True
    cea_py_product_species: Tuple[str, ...] = ()
    cea_py_omit_species: Tuple[str, ...] = ()
    cea_py_insert_species: Tuple[str, ...] = ()
    cea_py_molecular_weight: str = "gas_phase_M"   # gas_phase_M | total_MW
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

    # ---------------- nozzle low-pressure policy ----------------
    unchoked_policy: str = "switch"  # 'switch' (subsonic branch) | 'stop'

    # ---------------- solver ----------------
    method: str = "LSODA"         # LSODA | BDF | Radau
    rtol: float = 1.0e-9
    atol_p: float = 1.0e-3        # [Pa]
    atol_x: float = 1.0e-13       # [m]
    t_max: float = 600.0          # integration horizon of the burning phase [s]
    max_step_burn: float = 0.02   # [s]
    blowdown: bool = True         # integrate the post-burnout blowdown
    blowdown_tmax: float = 5.0    # [s]
    max_step_blow: float = 5.0e-4  # [s]


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
                       "Python package) | 'cea_legacy_executable' (external fcea2)"),
    "cea_pressure_points": "number of logarithmic chamber-pressure grid points",
    "cea_cache_enabled": "reuse the reproducible chamber-property cache",
    "cea_rebuild_cache": "ignore and rebuild an existing chamber-property cache",
    "cea_cache_dir": "directory for reproducible chamber-property caches",
    "cea_py_ions": "official CEA backend: include ionized species",
    "cea_py_transport": "official CEA backend: compute transport properties",
    "cea_py_trace": "official CEA backend: trace-species threshold (<=0 uses the API default)",
    "cea_py_products_from_reactants": "official CEA backend: build the product set from the reactant elements",
    "cea_py_molecular_weight": "'gas_phase_M' (gas-phase M, as in the legacy .plt 'm' column) | 'total_MW'",
    "cea_legacy_executable_path": "path to the native CEA executable (compatibility backend)",
    "cea_legacy_data_directory": "directory containing thermo.lib and trans.lib",
    "cea_legacy_timeout_s": "timeout for each external CEA batch [s]",
    "cea_legacy_batch_size": "assigned pressures per external CEA run (8 = tested legacy workaround)",
    "cea_legacy_trace": "external CEA output trace threshold",
    "p_fit_min": "lower validity limit of the property table [Pa]",
    "p_fit_max": "upper validity limit of the property table [Pa]",
    "eta_T0": "combustion efficiency applied to T0 [-]",
    "property_policy": "'clamp' or 'extrapolate' outside the table range",
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


# ==============================================================================
# 2. DERIVED GEOMETRY HELPERS
# ==============================================================================
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


# ==============================================================================
# 3. THERMOCHEMISTRY BACKENDS
# ------------------------------------------------------------------------------
#  The internal-ballistics solver only ever asks a backend for:
#      props(p)                  -> (R [J/(kg K)], T0 [K], gamma_s [-])
#      theta_and_derivative(p)   -> (Theta = R*T0 [J/kg], dTheta/dp [J/(kg Pa)])
#      is_extrapolated(p)        -> bool
#      metadata()                -> provenance dictionary for summary.json
#
#  Backends never modify grain geometry, burn rate, nozzle flow or the ODE.
# ==============================================================================

#: Actionable message required whenever the official package is missing.
MISSING_CEA_PACKAGE_MESSAGE = f"""The official NASA CEA Python package is required for the
'{BACKEND_CEA_PYTHON}' backend.

Install it with:
    {CEA_INSTALL_COMMAND}

Official project:
    {CEA_PROJECT_URL}"""

_TABLE_REQUIRED_KEYS = ("pressure_Pa", "temperature_K", "gamma_s",
                        "molecular_weight_kg_kmol", "gas_constant_J_kgK")


@dataclass(frozen=True)
class ReactantSpec:
    """One solid-propellant constituent (mass fraction, temperature, density)."""
    name: str
    wt_percent: float
    temperature_K: float
    density_kg_m3: float


def ideal_mixture_density(reactants: Sequence[ReactantSpec]) -> float:
    """Additive-volume (ideal) mixture density  rho = 1 / sum_i(w_i/rho_i)."""
    specific_volume = sum((r.wt_percent / 100.0) / r.density_kg_m3 for r in reactants)
    if not math.isfinite(specific_volume) or specific_volume <= 0.0:
        raise ThermochemistryError(
            "The reactant densities produced an invalid mixture specific volume.")
    return 1.0 / specific_volume


def _finite_positive(value: Any) -> bool:
    try:
        v = float(value)
    except (TypeError, ValueError):
        return False
    return math.isfinite(v) and v > 0.0


# ------------------------------------------------------------------------------
# 3.1 Backend interface
# ------------------------------------------------------------------------------
class ThermochemistryBackend:
    """Solver-facing interface.  Calls made by the ODE right-hand side must be side-effect free and inexpensive:
    implementations may be queried during every right-hand-side evaluation, so the
    solver-facing methods must use interpolation or closed-form evaluation only."""

    name = "abstract"
    #: True when the backend evaluates an external CEA calculation while building
    #: the chamber-property table (never inside the ODE right-hand side).
    uses_cea = False

    #: lower/upper validity limit of the property representation [Pa]
    pmin: float = 1.0e5
    pmax: float = 8.0e6

    def props(self, pressure_pa: float) -> Tuple[float, float, float]:
        raise NotImplementedError

    def theta_and_derivative(self, pressure_pa: float) -> Tuple[float, float]:
        raise NotImplementedError

    def is_extrapolated(self, pressure_pa: float) -> bool:
        raise NotImplementedError

    def metadata(self) -> Dict[str, Any]:
        raise NotImplementedError

    # -- convenience ---------------------------------------------------------
    def theta(self, pressure_pa: float) -> float:
        return self.theta_and_derivative(pressure_pa)[0]

    @property
    def outside_range_policy(self) -> str:
        raise NotImplementedError

    def describe(self) -> str:
        return self.name


# ------------------------------------------------------------------------------
# 3.2 Shared pressure-table machinery (PCHIP in ln p + reproducible cache)
# ------------------------------------------------------------------------------
class PressureTableBackend(ThermochemistryBackend):
    """Base class for backends that tabulate chamber properties on a pressure grid.

    The table is built once, cached as JSON, and interpolated with monotonic PCHIP
    interpolators in z = ln(p).  The derivative of Theta = R*T0 is taken from the
    analytic derivative of the interpolator and converted with

        dTheta/dp = dTheta/dln(p) / p

    so that the ODE pressure equation never sees finite-difference noise and CEA is
    never called from inside solve_ivp.
    """

    cache_schema = 2
    table_stem = "thermo_table"
    #: human-readable identification of the equilibrium formulation (metadata)
    equilibrium_formulation = "unspecified"

    def __init__(self, c: Config, reactants: Sequence[ReactantSpec],
                 packing_fraction: float, cache_root: Path):
        self.c = c
        self.reactants = list(reactants)
        self.packing_fraction = float(packing_fraction)
        self.ideal_mixture_density = ideal_mixture_density(self.reactants)
        self.bulk_mixture_density = self.ideal_mixture_density * self.packing_fraction

        self.cache_root = Path(cache_root)
        self.key = self.cache_key()
        self.cache_path = self.cache_root / f"{self.table_stem}_{self.key}.json"
        self.cache_used = False
        self.cache_written = False

        if c.cea_cache_enabled and self.cache_path.is_file() and not c.cea_rebuild_cache:
            table = self._read_cache()
            self.cache_used = True
        else:
            table = self._build_table()
            if c.cea_cache_enabled:
                self.cache_root.mkdir(parents=True, exist_ok=True)
                self.cache_path.write_text(
                    json.dumps(table, indent=2, sort_keys=False) + "\n", encoding="utf-8")
                self.cache_written = True
        self.table = table
        self._load_table(table)
        self.property_validation = self._validate_table()

    # -- cache key -----------------------------------------------------------
    def cache_key_payload(self) -> Dict[str, Any]:
        """Everything that can change the tabulated properties.

        Solid constituent densities and the packing fraction affect rho_p only and
        are deliberately excluded (they are still reported in the provenance), so
        changing them does not invalidate the gas-property cache.
        """
        return {
            "cache_schema": self.cache_schema,
            "backend": self.name,
            "reactants": [[r.name, float(r.wt_percent), float(r.temperature_K)]
                          for r in self.reactants],
            "equilibrium": self.equilibrium_formulation,
            "pressure_min_Pa": float(self.c.p_fit_min),
            "pressure_max_Pa": float(self.c.p_fit_max),
            "pressure_points": int(self.c.cea_pressure_points),
            "temperature_efficiency": float(self.c.eta_T0),
            "options": self.option_signature(),
        }

    def option_signature(self) -> Dict[str, Any]:
        return {}

    def cache_key(self) -> str:
        payload = self.cache_key_payload()
        blob = json.dumps(payload, sort_keys=True, default=str).encode("utf-8")
        return hashlib.sha256(blob).hexdigest()[:24]

    # -- cache IO ------------------------------------------------------------
    def _read_cache(self) -> Dict[str, Any]:
        try:
            raw = self.cache_path.read_text(encoding="utf-8")
        except OSError as exc:
            raise ThermochemistryCacheError(
                f"Chamber-property cache could not be read: {self.cache_path}\n"
                f"Original error: {exc}\n"
                f"Set 'rebuild_cache' to true for this backend, or delete the file."
            ) from exc
        try:
            table = json.loads(raw)
        except json.JSONDecodeError as exc:
            raise ThermochemistryCacheError(
                f"Chamber-property cache is not valid JSON: {self.cache_path}\n"
                f"Original error: {exc}\n"
                f"Set 'rebuild_cache' to true for this backend, or delete the file."
            ) from exc
        if not isinstance(table, dict):
            raise ThermochemistryCacheError(
                f"Chamber-property cache is malformed (expected a JSON object): {self.cache_path}")
        missing = [k for k in _TABLE_REQUIRED_KEYS if k not in table]
        if missing:
            raise ThermochemistryCacheError(
                f"Chamber-property cache is missing {missing}: {self.cache_path}")
        cached_key = table.get("cache_key")
        if cached_key is not None and cached_key != self.key:
            raise ThermochemistryCacheError(
                f"Chamber-property cache key mismatch: {self.cache_path}\n"
                f"  cached:   {cached_key}\n"
                f"  expected: {self.key}\n"
                f"Set 'rebuild_cache' to true for this backend, or delete the file.")
        if table.get("backend") not in (None, self.name):
            raise ThermochemistryCacheError(
                f"Chamber-property cache belongs to backend {table.get('backend')!r} "
                f"but {self.name!r} was selected: {self.cache_path}")
        p = np.asarray(table["pressure_Pa"], dtype=float)
        npts = int(self.c.cea_pressure_points)
        if p.size != npts:
            raise ThermochemistryCacheError(
                f"Chamber-property cache holds {p.size} pressure points but "
                f"{npts} were configured: {self.cache_path}")
        if not np.all(np.isfinite(p)) or not np.all(np.diff(p) > 0.0):
            raise ThermochemistryCacheError(
                f"Chamber-property cache pressure grid is not strictly increasing: "
                f"{self.cache_path}")
        return table

    # -- table -> interpolators ---------------------------------------------
    def _load_table(self, table: Dict[str, Any]) -> None:
        p = np.asarray(table["pressure_Pa"], dtype=float)
        R = np.asarray(table["gas_constant_J_kgK"], dtype=float)
        M = np.asarray(table["molecular_weight_kg_kmol"], dtype=float)
        T = np.asarray(table["temperature_K"], dtype=float) * self.c.eta_T0
        g = np.asarray(table["gamma_s"], dtype=float)
        if p.ndim != 1 or R.shape != p.shape or T.shape != p.shape or g.shape != p.shape:
            raise ThermochemistryCacheError(
                "Chamber-property table arrays have inconsistent shapes.")
        if len(p) < 4 or not np.all(np.diff(p) > 0.0):
            raise ThermochemistryError(
                "The chamber-property pressure grid must contain at least four "
                "strictly increasing pressures.")
        if not (np.all(np.isfinite(p)) and np.all(p > 0.0)):
            raise ThermochemistryError(
                "The chamber-property pressure grid contains non-finite or "
                "non-positive pressures.")
        if not np.allclose(R * M, R_UNIVERSAL, rtol=1e-8):
            raise ThermochemistryError(
                "Inconsistent units in the chamber-property table: R * M differs "
                "from the universal gas constant (Ru = %.5f J/(kmol K))." % R_UNIVERSAL)
        if not (np.all(np.isfinite(T)) and np.all(T > 0.0)):
            raise ThermochemistryError(
                "The chamber-property table contains non-physical temperatures.")
        if not (np.all(np.isfinite(g)) and np.all(g > 1.0)):
            raise ThermochemistryError(
                "The chamber-property table contains non-physical gamma_s (must be > 1).")
        self.pmin = float(p[0])
        self.pmax = float(p[-1])
        z = np.log(p)
        self.Ri = PchipInterpolator(z, R, extrapolate=True)
        self.Ti = PchipInterpolator(z, T, extrapolate=True)
        self.gi = PchipInterpolator(z, g, extrapolate=True)
        self.thi = PchipInterpolator(z, R * T, extrapolate=True)
        self.dthi = self.thi.derivative()

    def _validate_table(self) -> Dict[str, Any]:
        """Numerical sanity report for the table (recorded in summary.json)."""
        p = np.asarray(self.table["pressure_Pa"], dtype=float)
        theta = np.asarray(self.table["gas_constant_J_kgK"], dtype=float) * \
            np.asarray(self.table["temperature_K"], dtype=float) * self.c.eta_T0
        # p * Theta'(p) / Theta = (dTheta/dln p) / Theta   [dimensionless]
        dz = np.diff(np.log(p))
        slope = np.diff(theta) / dz
        factor = 1.0 - slope / theta[:-1]
        report = {
            "activation_clamped": self.c.property_policy == "clamp",
            "pressure_points": int(p.size),
            "pressure_min_Pa": float(p[0]),
            "pressure_max_Pa": float(p[-1]),
            "temperature_min_K": float(np.min(self.table["temperature_K"])),
            "temperature_max_K": float(np.max(self.table["temperature_K"])),
            "gamma_min": float(np.min(self.table["gamma_s"])),
            "gamma_max": float(np.max(self.table["gamma_s"])),
            "molecular_weight_min_kg_kmol": float(np.min(self.table["molecular_weight_kg_kmol"])),
            "molecular_weight_max_kg_kmol": float(np.max(self.table["molecular_weight_kg_kmol"])),
            "theta_min_J_kg": float(np.min(theta)),
            "theta_max_J_kg": float(np.max(theta)),
            "min_pressure_factor_1_minus_pTheta_over_Theta": float(np.min(factor)),
            "all_points_finite": bool(np.all(np.isfinite(p))),
            "pressure_strictly_increasing": bool(np.all(np.diff(p) > 0.0)),
        }
        if report["min_pressure_factor_1_minus_pTheta_over_Theta"] <= 0.0:
            raise ThermochemistryError(
                "The chamber-pressure equation is non-physical over the configured "
                "property table: min(1 - p*Theta'/Theta) = "
                f"{report['min_pressure_factor_1_minus_pTheta_over_Theta']:.6g} <= 0.\n"
                "Check the thermochemistry inputs (reactants, temperatures, "
                "pressure range) before continuing.")
        return report

    # -- solver-facing API ---------------------------------------------------
    def _effective_pressure(self, p: float) -> float:
        if self.c.property_policy == "clamp":
            return min(max(p, self.pmin), self.pmax)
        return max(p, 1.0)

    def _outside_range(self, p: float) -> bool:
        return bool(p < self.pmin or p > self.pmax)

    def props(self, p: float) -> Tuple[float, float, float]:
        pe = self._effective_pressure(p)
        z = math.log(pe)
        R = float(self.Ri(z))
        T = float(self.Ti(z))
        g = float(self.gi(z))
        if not (_finite_positive(R) and _finite_positive(T) and math.isfinite(g) and g > 1.0):
            raise ThermochemistryError(
                f"Non-physical interpolated gas properties at p = {p:.6g} Pa: "
                f"R = {R:.6g}, T0 = {T:.6g}, gamma_s = {g:.6g}.")
        return R, T, g

    def theta_and_derivative(self, p: float) -> Tuple[float, float]:
        pe = self._effective_pressure(p)
        z = math.log(pe)
        theta = float(self.thi(z))
        clamped = self.c.property_policy == "clamp" and p != pe
        if clamped:
            dtheta = 0.0
        else:
            dtheta = float(self.dthi(z)) / pe
        if not (math.isfinite(theta) and theta > 0.0 and math.isfinite(dtheta)):
            raise ThermochemistryError(
                f"Non-physical Theta or dTheta/dp at p = {p:.6g} Pa "
                f"(Theta = {theta:.6g} J/kg, dTheta/dp = {dtheta:.6g}).")
        return theta, dtheta

    def is_extrapolated(self, p: float) -> bool:
        return self._outside_range(p)

    @property
    def outside_range_policy(self) -> str:
        return self.c.property_policy

    # -- provenance ----------------------------------------------------------
    def base_metadata(self) -> Dict[str, Any]:
        return {
            "backend": self.name,
            "configuration_schema_version": SCHEMA_VERSION,
            "equilibrium_formulation": self.equilibrium_formulation,
            "pressure_grid": {
                "minimum_Pa": float(self.c.p_fit_min),
                "maximum_Pa": float(self.c.p_fit_max),
                "points": int(self.c.cea_pressure_points),
                "spacing": "logarithmic (geometric)",
            },
            "outside_range_policy": self.c.property_policy,
            "temperature_efficiency": float(self.c.eta_T0),
            "interpolation": {
                "variable": "z = ln(p)",
                "method": "PCHIP (monotone cubic Hermite, scipy PchipInterpolator)",
                "quantities": ["R", "T0", "gamma_s", "Theta = R*T0"],
                "derivative": "dTheta/dp = d/dz[PCHIP(Theta)](z) / p",
            },
            "cache_key": self.key,
            "cache_file": str(self.cache_path),
            "cache_enabled": bool(self.c.cea_cache_enabled),
            "cache_used": bool(self.cache_used),
            "cache_written": bool(self.cache_written),
            "reactants": [asdict(r) for r in self.reactants],
            "mass_fractions_percent": [float(r.wt_percent) for r in self.reactants],
            "reactant_temperatures_K": [float(r.temperature_K) for r in self.reactants],
            "molecular_weight_units": "kg/kmol (numerically identical to g/mol)",
            "universal_gas_constant_J_kmolK": R_UNIVERSAL,
            "specific_gas_constant_definition": "R(p) = Ru / M(p)",
            "ideal_mixture_density_kg_m3": float(self.ideal_mixture_density),
            "packing_fraction": float(self.packing_fraction),
            "bulk_propellant_density_kg_m3": float(self.bulk_mixture_density),
            "propellant_density_model": ("rho_p = packing_fraction / sum_i(w_i/rho_i); "
                                         "independent of any CEA gas density"),
            "property_validation": dict(self.property_validation),
        }


# ------------------------------------------------------------------------------
# 3.3 cea_legacy_executable backend (transitional compatibility only)
# ------------------------------------------------------------------------------
class CEALegacyExecutableBackend(PressureTableBackend):
    """NASA CEA2 HP equilibrium through the external ``fcea2`` executable.

    This is the exact v0.3.1-alpha pathway: discover the executable and
    thermo.lib/trans.lib, write a legacy CEA input deck, run it in a temporary
    directory, parse the generated .plt file, and interpolate the result.

    It exists for regression comparison during the v0.4-alpha transition.  New
    work should select ``cea_python``.
    """

    name = BACKEND_CEA_LEGACY_EXECUTABLE
    uses_cea = True
    table_stem = "thermo_cea_legacy"
    equilibrium_formulation = ("NASA CEA HP (assigned enthalpy and pressure) equilibrium "
                               "through the external fcea2 executable; reactant enthalpy "
                               "from NAME reactants at their specified temperatures")

    def option_signature(self) -> Dict[str, Any]:
        c = self.c
        return {
            "executable": str(c.cea_legacy_executable_path),
            "data_directory": str(c.cea_legacy_data_directory),
            "batch_size": int(c.cea_legacy_batch_size),
            "trace": float(c.cea_legacy_trace),
            "timeout_s": float(c.cea_legacy_timeout_s),
        }

    # -- discovery -----------------------------------------------------------
    def _resolve_executable(self) -> str:
        p = Path(self.c.cea_legacy_executable_path).expanduser()
        candidates = [p] if p.is_absolute() else [
            Path(__file__).resolve().parent / p,
            Path.cwd() / p,
        ]
        for candidate in candidates:
            if candidate.is_file():
                candidate = candidate.resolve()
                if not os.access(candidate, os.X_OK):
                    raise BackendUnavailableError(
                        f"CEA executable exists but is not executable: {candidate}\n"
                        f"Grant execute permission with: chmod +x '{candidate}'")
                return str(candidate)
        exe = shutil.which(str(p))
        if exe:
            return exe
        searched = ", ".join(str(v) for v in candidates)
        raise BackendUnavailableError(
            f"The '{BACKEND_CEA_LEGACY_EXECUTABLE}' backend requires an external CEA "
            f"executable.\n"
            f"Executable not found: {self.c.cea_legacy_executable_path}\n"
            f"Searched: {searched}, and the system PATH.\n"
            f"Use the '{BACKEND_CEA_PYTHON}' backend (official Python package) to avoid "
            f"external executables entirely.")

    def _data_dir(self) -> Path:
        candidates = []
        if self.c.cea_legacy_data_directory:
            candidates.append(Path(self.c.cea_legacy_data_directory).expanduser())
        candidates.append(Path(self._resolve_executable()).parent)
        for d in candidates:
            if (d / "thermo.lib").is_file():
                return d.resolve()
        raise BackendUnavailableError(
            f"The '{BACKEND_CEA_LEGACY_EXECUTABLE}' backend requires the CEA "
            f"thermodynamic database.\n"
            f"thermo.lib was not found in cea_legacy_executable.data_directory "
            f"({self.c.cea_legacy_data_directory!r}) or beside the executable.\n"
            f"The '{BACKEND_CEA_PYTHON}' backend ships its own database and does not "
            f"need thermo.lib or trans.lib.")

    def _input_text(self, pressures_bar: Sequence[float]) -> str:
        lines = ["problem hp",
                 "  p,bar = " + " ".join(f"{p:.12g}" for p in pressures_bar),
                 f"  trace = {self.c.cea_legacy_trace:.6e}", "reactants"]
        # NAME reactants allow arbitrary solid-propellant formulations without
        # forcing an oxidizer/fuel partition or an O/F reinterpretation.
        for r in self.reactants:
            lines.append(f"  name = {r.name} wt% = {r.wt_percent:.12g} "
                         f"t,k = {r.temperature_K:.12g}")
        lines += ["output short", "  plot p t gam m", "end", ""]
        return "\n".join(lines)

    def _run_batch(self, pressures_pa: np.ndarray, ibatch: int):
        exe = self._resolve_executable()
        data = self._data_dir()
        with tempfile.TemporaryDirectory(prefix="gribs_cea2_") as td:
            td = Path(td)
            shutil.copy2(data / "thermo.lib", td / "thermo.lib")
            if (data / "trans.lib").is_file():
                shutil.copy2(data / "trans.lib", td / "trans.lib")
            stem = f"gribs_{ibatch:04d}"
            pbar = [p / 1e5 for p in pressures_pa]
            (td / f"{stem}.inp").write_text(self._input_text(pbar), encoding="ascii")
            try:
                cp = subprocess.run([exe], input=stem + "\n", text=True, cwd=td,
                                    capture_output=True,
                                    timeout=self.c.cea_legacy_timeout_s)
            except subprocess.TimeoutExpired as exc:
                raise ThermochemistryError(
                    f"The external CEA executable timed out after "
                    f"{self.c.cea_legacy_timeout_s:g} s (batch {ibatch}).") from exc
            if cp.returncode != 0:
                raise ThermochemistryError(
                    f"The external CEA executable failed (batch {ibatch}, "
                    f"rc={cp.returncode}): {cp.stderr[-2000:]}")
            plt_file = td / f"{stem}.plt"
            out_file = td / f"{stem}.out"
            if not plt_file.is_file():
                tail = (out_file.read_text(errors="replace")[-4000:]
                        if out_file.is_file() else cp.stdout[-4000:])
                raise ThermochemistryError(
                    "The external CEA executable did not create a plot file. "
                    "Output tail:\n" + tail)
            rows = []
            for line in plt_file.read_text(errors="replace").splitlines():
                if not line.strip() or line.lstrip().startswith("#"):
                    continue
                vals = line.replace("D", "E").split()
                if len(vals) >= 4:
                    try:
                        rows.append(tuple(map(float, vals[:4])))
                    except ValueError:
                        pass
            if len(rows) != len(pressures_pa):
                tail = (out_file.read_text(errors="replace")[-4000:]
                        if out_file.is_file() else "")
                raise ThermochemistryError(
                    f"The external CEA executable returned {len(rows)} points, "
                    f"expected {len(pressures_pa)}.\n{tail}")
            return rows

    def _build_table(self) -> Dict[str, Any]:
        p = np.geomspace(self.c.p_fit_min, self.c.p_fit_max,
                         self.c.cea_pressure_points)
        rows: List[Tuple[float, float, float, float]] = []
        # Use eight assigned-pressure points per CEA run.  Legacy/distributed
        # fcea2 builds may reserve array slots internally and return fewer than
        # ten points when ten are requested.  Eight is compatible with the
        # observed formatted-output width and does not change the pressure grid.
        batch = int(self.c.cea_legacy_batch_size)
        for i in range(0, len(p), batch):
            rows.extend(self._run_batch(p[i:i + batch], i // batch))
        arr = np.asarray(rows, dtype=float)
        order = np.argsort(arr[:, 0])
        arr = arr[order]
        pp = arr[:, 0] * 1e5
        T = arr[:, 1]
        gamma = arr[:, 2]
        M = arr[:, 3]
        R = R_UNIVERSAL / M
        if not (np.all(np.isfinite(arr)) and np.all(R > 0) and np.all(T > 0)
                and np.all(gamma > 1)):
            raise ThermochemistryError(
                "The external CEA executable produced non-physical chamber properties.")
        return {
            "schema": self.cache_schema,
            "backend": self.name,
            "cache_key": self.key,
            "source": "external NASA CEA executable, HP equilibrium, .plt output",
            "cache_role": "TRANSITIONAL COMPATIBILITY BACKEND (v0.4-alpha only)",
            "pressure_Pa": pp.tolist(),
            "temperature_K": T.tolist(),
            "gamma_s": gamma.tolist(),
            "molecular_weight_kg_kmol": M.tolist(),
            "gas_constant_J_kgK": R.tolist(),
        }

    def option_metadata(self) -> Dict[str, Any]:
        try:
            exe = self._resolve_executable()
        except BackendUnavailableError:
            exe = None
        try:
            data_dir = str(self._data_dir())
        except BackendUnavailableError:
            data_dir = None
        thermo = Path(data_dir) / "thermo.lib" if data_dir else None
        return {
            "compatibility_backend": True,
            "warning": ("Transitional backend retained from v0.3.1-alpha for regression "
                        "comparison; prefer 'cea_python' for new work."),
            "executable": exe or self.c.cea_legacy_executable_path,
            "data_directory": data_dir,
            "thermo_lib": (str(thermo) if thermo is not None and thermo.is_file() else None),
            "thermo_lib_size_bytes": (thermo.stat().st_size if thermo is not None
                                      and thermo.is_file() else None),
            "batch_size": int(self.c.cea_legacy_batch_size),
            "trace": float(self.c.cea_legacy_trace),
            "timeout_s": float(self.c.cea_legacy_timeout_s),
        }

    def metadata(self) -> Dict[str, Any]:
        md = self.base_metadata()
        md.update(self.option_metadata())
        md["cea_package"] = None
        return md


# ------------------------------------------------------------------------------
# 3.4 cea_python backend (official NASA CEA Python package)
# ------------------------------------------------------------------------------
def import_official_cea_module():
    """Import and validate the official NASA CEA Python package.

    Raises
    ------
    BackendUnavailableError
        with an actionable message when the package is missing, when a different
        package named ``cea`` is installed, or when the installed official API
        version is not supported.
    """
    try:
        module = importlib.import_module("cea")
    except Exception as exc:                      # ImportError, OSError, ...
        message = (f"{MISSING_CEA_PACKAGE_MESSAGE}\n\n"
                   f"Original import error: {type(exc).__name__}: {exc}")
        raise BackendUnavailableError(message) from exc

    version = getattr(module, "__version__", None)
    required_attrs = ("Mixture", "EqSolver", "EqSolution", "RocketSolver",
                      "RocketSolution", "HP", "ENTHALPY", "R")
    missing = [a for a in required_attrs if not hasattr(module, a)]
    if version is None or missing:
        raise BackendUnavailableError(
            f"The module importable as 'cea' is not the official NASA CEA Python "
            f"package.\n"
            f"  module path:  {getattr(module, '__file__', '<unknown>')}\n"
            f"  __version__:  {version!r}\n"
            f"  missing API:  {missing or 'none'}\n\n"
            f"{MISSING_CEA_PACKAGE_MESSAGE}")

    numbers: List[int] = []
    for part in str(version).split("."):
        digits = "".join(ch for ch in part if ch.isdigit())
        if not digits:
            break
        numbers.append(int(digits))
    numbers += [0] * (2 - len(numbers))
    if tuple(numbers[:2]) < CEA_MIN_API_VERSION:
        raise BackendUnavailableError(
            f"Unsupported official NASA CEA Python API version: {version}\n"
            f"GRIBS v{PROGRAM_VERSION} requires the CEA Python API >= "
            f"{'.'.join(str(v) for v in CEA_MIN_API_VERSION)} "
            f"(validated with {'.'.join(str(v) for v in CEA_TESTED_API_VERSION)}.x).\n"
            f"Upgrade with:\n    {CEA_INSTALL_COMMAND}")
    return module


def _cea_package_metadata() -> Dict[str, Any]:
    info: Dict[str, Any] = {"distribution_name": None, "distribution_version": None,
                            "summary": None, "author_email": None, "license": None,
                            "project_url": CEA_PROJECT_URL}
    try:
        meta = importlib_metadata.metadata("cea")
    except Exception:
        return info
    info["distribution_name"] = meta.get("Name")
    info["distribution_version"] = meta.get("Version")
    info["summary"] = meta.get("Summary")
    info["author_email"] = meta.get("Author-email")
    info["license"] = meta.get("License-Expression") or meta.get("License")
    for key in ("Project-URL", "Home-page"):
        value = meta.get(key)
        if value:
            info["project_url"] = value
            break
    return info


def _cea_database_identity(cea_module: Any) -> Dict[str, Any]:
    """Identify the thermodynamic database that ships inside the official package.

    The digest is part of the chamber-property cache key, so replacing or
    upgrading the package database can never silently reuse stale properties.
    """
    info: Dict[str, Any] = {
        "thermo_lib_path": None, "thermo_lib_size_bytes": None,
        "thermo_lib_sha256_prefix": None, "trans_lib_path": None,
        "trans_lib_size_bytes": None,
    }
    try:
        root = Path(cea_module.__file__).resolve().parent
    except Exception:                              # pragma: no cover - defensive
        return info
    for key, name in (("thermo_lib", "thermo.lib"), ("trans_lib", "trans.lib")):
        path = root / "data" / name
        if not path.is_file():
            continue
        data = path.read_bytes()
        info[f"{key}_path"] = str(path)
        info[f"{key}_size_bytes"] = len(data)
        info[f"{key}_sha256_prefix"] = hashlib.sha256(data).hexdigest()[:16]
    return info


class CEAPythonBackend(PressureTableBackend):
    """Official NASA CEA Python package backend (preferred for v0.4).

    Equilibrium formulation
    -----------------------
    * Problem type: **HP** - assigned pressure and assigned enthalpy, i.e. an
      adiabatic, constant-pressure chamber equilibrium (exactly the formulation
      used by the v0.3 ``problem hp`` input deck).  Pyrolysis/heat losses are not
      modelled by CEA; the optional GRIBS ``temperature_efficiency`` is applied
      afterwards, on the tabulated T0 only.
    * Reactant enthalpy: computed by ``Mixture.calc_property(cea.ENTHALPY, w, T)``
      with the **per-constituent** temperatures from the configuration and the
      mass fractions of the formulation.  The value is passed to the official
      solver as ``h/R`` [K] (the convention of the official examples).
    * Pressure units: Pa in the configuration, converted to **bar** for the API.
    * Temperature units: K.
    * Molecular weight: ``EqSolution.M`` by default - the gas-phase molecular
      weight in kg/kmol (= g/mol), i.e. the same quantity the legacy ``plot m``
      column provided.  ``total_MW`` (including condensed species) is selectable.
    * gamma: ``EqSolution.gamma_s`` (isentropic exponent of the equilibrium gas).
    * Specific gas constant: R = Ru/M with Ru = 8314.51 J/(kmol K) - the value of
      ``cea.R``, verified at runtime.
    * Condensed species: included in the equilibrium (GRIBS only consumes the
      gas-phase M and gamma_s; the condensed mass is reported as a diagnostic).
    * Ions: off by default (``ions`` option).  Transport: off by default.
    * Product species: built from the reactant elements
      (``Mixture(..., products_from_reactants=True)``), which reproduces the
      legacy "no explicit product list" behaviour, or given explicitly through
      ``product_species``.

    The chamber state is obtained from ``EqSolver`` only; the GRIBS unsteady
    internal-ballistics and nozzle model is untouched.  ``RocketSolver`` is used
    exclusively for clearly-named theoretical comparison values in summary.json.
    """

    name = BACKEND_CEA_PYTHON
    uses_cea = True
    table_stem = "thermo_cea_python"
    equilibrium_formulation = ("NASA CEA HP equilibrium (assigned pressure and enthalpy) "
                               "at each tabulated chamber pressure; reactant enthalpy from "
                               "the configured constituent temperatures")

    def __init__(self, c: Config, reactants: Sequence[ReactantSpec],
                 packing_fraction: float, cache_root: Path):
        # the option guards run before the base-class table build, so the
        # configuration must already be reachable
        self.c = c
        self.reactants = list(reactants)
        self.cea = import_official_cea_module()
        self.cea_metadata = _cea_package_metadata()
        self.cea_database = _cea_database_identity(self.cea)
        self._check_gas_constant()
        self._check_option_consistency()
        super().__init__(c, reactants, packing_fraction, cache_root)

    # -- guards --------------------------------------------------------------
    def _check_gas_constant(self) -> None:
        cea_R = float(self.cea.R)
        if not math.isfinite(cea_R) or abs(cea_R - R_UNIVERSAL) > 1e-6 * R_UNIVERSAL:
            raise ThermochemistryError(
                f"Inconsistent units: cea.R = {cea_R} J/(kmol K) but GRIBS expects "
                f"{R_UNIVERSAL} J/(kmol K).\n"
                f"Check the installed official CEA package version "
                f"({getattr(self.cea, '__version__', 'unknown')}).")

    def _check_option_consistency(self) -> None:
        c = self.c
        if c.cea_py_molecular_weight not in ("gas_phase_M", "total_MW"):
            raise ConfigurationError(
                "thermochemistry.cea_python.molecular_weight must be "
                "'gas_phase_M' or 'total_MW'.")
        if c.cea_py_products_from_reactants and c.cea_py_product_species:
            raise ConfigurationError(
                "thermochemistry.cea_python: 'product_species' must be empty when "
                "'products_from_reactants' is true.")

    def option_signature(self) -> Dict[str, Any]:
        c = self.c
        return {
            "cea_package_version": str(getattr(self.cea, "__version__", "unknown")),
            "cea_library_version": str(self.cea.lib_version()),
            "cea_module_path": str(getattr(self.cea, "__file__", "unknown")),
            "ions": bool(c.cea_py_ions),
            "transport": bool(c.cea_py_transport),
            "trace": float(c.cea_py_trace),
            "products_from_reactants": bool(c.cea_py_products_from_reactants),
            "product_species": list(c.cea_py_product_species),
            "omit_species": list(c.cea_py_omit_species),
            "insert_species": list(c.cea_py_insert_species),
            "molecular_weight": str(c.cea_py_molecular_weight),
            "smooth_truncation": bool(c.cea_py_smooth_truncation),
            "truncation_width": float(c.cea_py_truncation_width),
            "equilibrium_type": "HP",
            "pressure_unit_for_api": "bar",
            "enthalpy_convention": "h/R [K]",
            "database": "package-internal thermo.lib / trans.lib",
            "database_identity": dict(self.cea_database),
        }

    # -- mixture construction -----------------------------------------------
    def _build_mixtures(self):
        cea = self.cea
        c = self.c
        names = [r.name for r in self.reactants]
        weights = np.array([r.wt_percent for r in self.reactants], dtype=float)
        temperatures = np.array([r.temperature_K for r in self.reactants], dtype=float)
        try:
            reactant_mixture = cea.Mixture(names)
        except Exception as exc:
            raise ThermochemistryError(
                "CEA rejected the reactant mixture.\n"
                f"  reactants: {names}\n"
                f"  original error: {type(exc).__name__}: {exc}\n"
                "Check for a mis-spelled species name, a missing phase designator "
                "(the official database distinguishes e.g. 'C' (gaseous carbon), "
                "'C(gr)' (graphite) and does not contain 'C(L)'), or a species that "
                "does not exist in the installed thermo.lib.") from exc

        if c.cea_py_products_from_reactants:
            try:
                product_mixture = cea.Mixture(
                    names, products_from_reactants=True,
                    omit=list(c.cea_py_omit_species), ions=bool(c.cea_py_ions))
            except Exception as exc:
                raise ThermochemistryError(
                    "CEA could not construct the product species set from the "
                    f"reactant elements (reactants: {names}).\n"
                    f"  original error: {type(exc).__name__}: {exc}\n"
                    "Set 'products_from_reactants' to false and provide an explicit "
                    "'product_species' list.") from exc
        else:
            species = list(c.cea_py_product_species)
            if not species:
                raise ConfigurationError(
                    "thermochemistry.cea_python.product_species must list at least one "
                    "species when 'products_from_reactants' is false.")
            try:
                product_mixture = cea.Mixture(
                    species, omit=list(c.cea_py_omit_species),
                    ions=bool(c.cea_py_ions))
            except Exception as exc:
                raise ThermochemistryError(
                    f"CEA could not construct the requested product mixture: {species}\n"
                    f"  original error: {type(exc).__name__}: {exc}") from exc

        kwargs: Dict[str, Any] = {
            "reactants": reactant_mixture,
            "transport": bool(c.cea_py_transport),
            "ions": bool(c.cea_py_ions),
            "trace": float(c.cea_py_trace),
            "insert": list(c.cea_py_insert_species),
        }
        if c.cea_py_smooth_truncation:
            kwargs["smooth_truncation"] = True
            kwargs["truncation_width"] = float(c.cea_py_truncation_width)
        try:
            solver = cea.EqSolver(product_mixture, **kwargs)
        except Exception as exc:
            raise ThermochemistryError(
                "The official CEA equilibrium solver could not be constructed.\n"
                f"  options: { {k: v for k, v in kwargs.items() if k != 'reactants'} }\n"
                f"  original error: {type(exc).__name__}: {exc}") from exc

        solution = cea.EqSolution(solver)
        try:
            h_c = float(reactant_mixture.calc_property(cea.ENTHALPY, weights, temperatures))
        except Exception as exc:
            raise ThermochemistryError(
                "CEA could not evaluate the reactant-mixture enthalpy.\n"
                f"  reactants: {names}\n"
                f"  mass fractions [%]: {weights.tolist()}\n"
                f"  temperatures [K]: {temperatures.tolist()}\n"
                f"  original error: {type(exc).__name__}: {exc}\n"
                "Check the constituent temperatures against the valid ranges of the "
                "database entries.") from exc
        if not math.isfinite(h_c):
            raise ThermochemistryError(
                f"The reactant-mixture enthalpy is not finite: h = {h_c!r} J/kg.")
        return reactant_mixture, product_mixture, solver, solution, weights, h_c

    # -- table ---------------------------------------------------------------
    def _build_table(self) -> Dict[str, Any]:
        cea = self.cea
        c = self.c
        _, product_mixture, solver, solution, weights, h_c = self._build_mixtures()
        h_R = h_c / float(cea.R)

        pressures = np.geomspace(c.p_fit_min, c.p_fit_max, c.cea_pressure_points)
        T_out: List[float] = []
        M_out: List[float] = []
        MW_out: List[float] = []
        g_out: List[float] = []
        rho_out: List[float] = []
        cp_out: List[float] = []
        a_out: List[float] = []
        h_out: List[float] = []
        species_top: Dict[str, float] = {}

        for p_pa in pressures:
            p_bar = float(p_pa) / 1.0e5
            try:
                solver.solve(solution, cea.HP, h_R, p_bar, weights)
            except Exception as exc:
                raise ThermochemistryError(
                    f"CEA equilibrium solve failed at p = {p_pa:.8g} Pa "
                    f"({p_bar:.8g} bar).\n"
                    f"  reactants: {[r.name for r in self.reactants]}\n"
                    f"  mass fractions [%]: {weights.tolist()}\n"
                    f"  original error: {type(exc).__name__}: {exc}") from exc
            if not bool(solution.converged):
                raise ThermochemistryError(
                    f"CEA did not converge at p = {p_pa:.8g} Pa ({p_bar:.8g} bar); "
                    f"last_error = {int(solution.last_error)}.\n"
                    "Adjust the pressure range, the reactant set or the CEA options.")

            T = float(solution.T)
            M = float(solution.M)
            MW = float(solution.MW)
            gamma = float(solution.gamma_s)
            if not _finite_positive(T):
                raise ThermochemistryError(
                    f"CEA returned a non-physical chamber temperature "
                    f"T = {T!r} K at p = {p_pa:.8g} Pa.")
            molecular_weight = M if c.cea_py_molecular_weight == "gas_phase_M" else MW
            if not _finite_positive(molecular_weight):
                raise ThermochemistryError(
                    f"CEA returned a non-physical molecular weight "
                    f"({c.cea_py_molecular_weight}) = {molecular_weight!r} kg/kmol "
                    f"at p = {p_pa:.8g} Pa.")
            if not (math.isfinite(gamma) and gamma > 1.0):
                raise ThermochemistryError(
                    f"CEA returned a non-physical isentropic exponent "
                    f"gamma_s = {gamma!r} at p = {p_pa:.8g} Pa.")

            T_out.append(T)
            M_out.append(molecular_weight)
            MW_out.append(MW)
            g_out.append(gamma)
            rho_out.append(float(solution.density))
            cp_out.append(float(solution.cp_eq))
            a_out.append(float(solution.sonic_velocity) if hasattr(solution, "sonic_velocity")
                         else float("nan"))
            h_out.append(float(solution.enthalpy))
            if p_pa == pressures[-1]:
                species_top = {k: float(v) for k, v in
                               sorted(solution.mass_fractions.items(),
                                      key=lambda kv: -kv[1])[:8]}

        R_out = [R_UNIVERSAL / m for m in M_out]
        return {
            "schema": self.cache_schema,
            "backend": self.name,
            "cache_key": self.key,
            "source": (f"official NASA CEA Python package {getattr(cea, '__version__', '?')} "
                       f"(HP equilibrium)"),
            "cea_version": str(getattr(cea, "__version__", "unknown")),
            "cea_library_version": str(cea.lib_version()),
            "cea_module_path": str(getattr(cea, "__file__", "unknown")),
            "pressure_Pa": pressures.tolist(),
            "temperature_K": T_out,
            "gamma_s": g_out,
            "molecular_weight_kg_kmol": M_out,
            "gas_constant_J_kgK": R_out,
            "reactant_enthalpy_J_kg": h_c,
            "reactant_enthalpy_over_R_K": h_R,
            "reactant_enthalpy_recomputed_kJ_kg": h_out,
            "molecular_weight_definition": c.cea_py_molecular_weight,
            "product_species_count": int(product_mixture.num_species),
            "num_gas": int(solver.num_gas),
            "num_condensed": int(solver.num_condensed),
            "num_elements": int(solver.num_elements),
            "diagnostics": {
                "density_kg_m3": rho_out,
                "cp_eq_kJ_kgK": cp_out,
                "sonic_velocity_m_s": a_out,
                "total_molecular_weight_kg_kmol": MW_out,
                "top_mass_fractions_at_max_pressure": species_top,
                "note": ("diagnostic values only; they never replace the GRIBS "
                         "internal-ballistics or nozzle results"),
            },
        }

    def option_metadata(self) -> Dict[str, Any]:
        c = self.c
        return {
            "cea_package": self.cea_metadata,
            "cea_python_module_path": str(getattr(self.cea, "__file__", "unknown")),
            "cea_python_version": str(getattr(self.cea, "__version__", "unknown")),
            "cea_library_version": str(self.cea.lib_version()),
            "ions": bool(c.cea_py_ions),
            "transport": bool(c.cea_py_transport),
            "trace": float(c.cea_py_trace),
            "products_from_reactants": bool(c.cea_py_products_from_reactants),
            "product_species": list(c.cea_py_product_species),
            "omit_species": list(c.cea_py_omit_species),
            "insert_species": list(c.cea_py_insert_species),
            "molecular_weight": str(c.cea_py_molecular_weight),
            "smooth_truncation": bool(c.cea_py_smooth_truncation),
            "truncation_width": float(c.cea_py_truncation_width),
            "product_species_count": self.table.get("product_species_count"),
            "num_gas": self.table.get("num_gas"),
            "num_condensed": self.table.get("num_condensed"),
            "num_elements": self.table.get("num_elements"),
            "reactant_enthalpy_J_kg": self.table.get("reactant_enthalpy_J_kg"),
            "reactant_enthalpy_over_R_K": self.table.get("reactant_enthalpy_over_R_K"),
            "equilibrium_constraint": "HP (assigned enthalpy and pressure)",
            "reactant_enthalpy_treatment": (
                "Mass-weighted constituent enthalpies at their individual configured "
                "temperatures; passed to CEA as h/R [K]."),
            "pressure_units": "configuration in Pa, official API in bar",
            "temperature_units": "K",
            "molecular_weight_definition": (
                "EqSolution.M (gas-phase, kg/kmol)") if c.cea_py_molecular_weight
                == "gas_phase_M" else "EqSolution.MW (incl. condensed, kg/kmol)",
            "gamma_definition": "EqSolution.gamma_s (equilibrium isentropic exponent)",
            "heat_capacity_treatment": "equilibrium (cp_eq reported as a diagnostic)",
            "condensed_species": ("included in the equilibrium; GRIBS consumes the "
                                  "gas-phase M and gamma_s only"),
            "ion_settings": "enabled" if c.cea_py_ions else "disabled",
            "transport_settings": "enabled" if c.cea_py_transport else "disabled",
            "product_species_selection": ("built from the reactant elements"
                                          if c.cea_py_products_from_reactants
                                          else "explicit list from the configuration"),
            "rocket_solver_usage": ("diagnostic only: IAC chamber state used to report "
                                    "cea_theoretical_* comparison values; the GRIBS "
                                    "nozzle model is not replaced"),
            "finite_area_combustor": False,
            "external_files_required": False,
            "database": "package-internal thermo.lib / trans.lib",
            "database_identity": dict(self.cea_database),
        }

    def metadata(self) -> Dict[str, Any]:
        md = self.base_metadata()
        md.update(self.option_metadata())
        md["table_source"] = self.table.get("source")
        md["table_diagnostics_summary"] = {
            "density_kg_m3_min_max": [float(np.nanmin(self.table["diagnostics"]["density_kg_m3"])),
                                      float(np.nanmax(self.table["diagnostics"]["density_kg_m3"]))],
            "cp_eq_kJ_kgK_min_max": [float(np.nanmin(self.table["diagnostics"]["cp_eq_kJ_kgK"])),
                                     float(np.nanmax(self.table["diagnostics"]["cp_eq_kJ_kgK"]))],
            "top_mass_fractions_at_max_pressure":
                self.table["diagnostics"]["top_mass_fractions_at_max_pressure"],
        }
        return md

    # -- theoretical rocket diagnostics (never used by the GRIBS solver) -----
    def theoretical_rocket_performance(self, chamber_pressure_Pa: float,
                                       expansion_ratio: float) -> Dict[str, Any]:
        """CEA IAC theoretical rocket performance at one chamber pressure.

        The result is returned with ``cea_theoretical_*`` names and is never fed
        back into GRIBS.  For a converging nozzle (Ae/At = 1) the throat station is
        selected, because it coincides with the geometric exit.

        The reactant mixture is mandatory: ``RocketSolver`` also accepts the
        no-reactant overload, but empirical tests with cea 3.3.4 show that the
        solver then returns ``c_star = inf`` and ``Cf = 0`` while still reporting
        ``converged = True``.  The station arrays are therefore validated below
        and a non-physical result is reported as an error, never as a number.
        """
        cea = self.cea
        c = self.c
        reactant_mixture, product_mixture, _, _, weights, h_c = self._build_mixtures()
        try:
            rocket_solver = cea.RocketSolver(product_mixture, reactants=reactant_mixture,
                                             transport=False,
                                             trace=float(c.cea_py_trace))
        except Exception as exc:                   # pragma: no cover - defensive
            raise ThermochemistryError(
                "The official CEA rocket solver could not be constructed with an "
                "explicit reactant mixture.\n"
                f"  original error: {type(exc).__name__}: {exc}") from exc
        rocket_solution = cea.RocketSolution(rocket_solver)
        eps = float(expansion_ratio)
        supar_request = eps if eps > 1.0 + 1.0e-12 else 1.0 + 1.0e-6
        pc_bar = float(chamber_pressure_Pa) / 1.0e5
        try:
            rocket_solver.solve(rocket_solution, weights, pc_bar,
                                supar=[supar_request], hc=h_c / float(cea.R), iac=True)
        except Exception as exc:
            raise ThermochemistryError(
                f"CEA theoretical rocket evaluation failed at p_c = "
                f"{chamber_pressure_Pa:.8g} Pa with Ae/At = {eps:g}.\n"
                f"  original error: {type(exc).__name__}: {exc}") from exc
        if not bool(rocket_solution.converged):
            raise ThermochemistryError(
                f"CEA theoretical rocket solution did not converge at "
                f"p_c = {chamber_pressure_Pa:.8g} Pa.")

        P = np.asarray(rocket_solution.P, dtype=float) * 1.0e5
        T = np.asarray(rocket_solution.T, dtype=float)
        mach = np.asarray(rocket_solution.Mach, dtype=float)
        ae_at = np.asarray(rocket_solution.ae_at, dtype=float)
        c_star = np.asarray(rocket_solution.c_star, dtype=float)
        cf = np.asarray(rocket_solution.coefficient_of_thrust, dtype=float)
        isp = np.asarray(rocket_solution.Isp, dtype=float)
        isp_vac = np.asarray(rocket_solution.Isp_vacuum, dtype=float)
        mw = np.asarray(rocket_solution.M, dtype=float)
        gamma = np.asarray(rocket_solution.gamma_s, dtype=float)

        # guard against a silently unusable rocket solution (see the docstring)
        for label, arr in (("P", P), ("T", T), ("Mach", mach), ("ae_at", ae_at),
                           ("c_star", c_star), ("coefficient_of_thrust", cf),
                           ("Isp", isp), ("Isp_vacuum", isp_vac), ("M", mw),
                           ("gamma_s", gamma)):
            if not np.all(np.isfinite(arr)):
                raise ThermochemistryError(
                    f"CEA returned a non-finite rocket-performance array '{label}' "
                    f"at p_c = {chamber_pressure_Pa:.8g} Pa and Ae/At = {eps:g}: "
                    f"{np.asarray(arr).tolist()}.\n"
                    "The official RocketSolver can report converged=True with an "
                    "unusable solution when it is not given an explicit reactant "
                    "mixture; verify that 'reactants' was supplied and that the "
                    "chamber pressure is inside the tabulated validity range.")
        if not (np.all(P > 0.0) and np.all(T > 0.0) and np.all(c_star > 0.0)):
            raise ThermochemistryError(
                f"CEA returned non-positive rocket-performance values at p_c = "
                f"{chamber_pressure_Pa:.8g} Pa: P = {P.tolist()}, T = {T.tolist()}, "
                f"c_star = {c_star.tolist()}.")

        # Station selection (never assume index 0 blindly): the chamber station is
        # the first station with Mach = 0 and no area ratio.
        chamber_candidates = np.flatnonzero(np.isclose(mach, 0.0, atol=1e-12))
        if chamber_candidates.size == 0:
            raise ThermochemistryError(
                "CEA theoretical rocket solution does not contain an identifiable "
                f"chamber station (Mach = {mach.tolist()}).")
        i_ch = int(chamber_candidates[0])

        if eps > 1.0 + 1.0e-12:
            i_ex = int(np.argmin(np.abs(ae_at - eps)))
            exit_kind = f"supersonic exit at Ae/At = {ae_at[i_ex]:.6g}"
        else:
            throat = np.flatnonzero(np.isclose(mach, 1.0, atol=1e-9) &
                                    np.isclose(ae_at, 1.0, atol=1e-9))
            i_ex = int(throat[0]) if throat.size else int(np.argmin(np.abs(mach - 1.0)))
            exit_kind = "throat (converging nozzle, Ae/At = 1)"

        return {
            "cea_theoretical_cstar_m_s": float(c_star[i_ch]),
            "cea_theoretical_cf": float(cf[i_ex]),
            "cea_theoretical_isp_s": float(isp[i_ex]) / G0,
            "cea_theoretical_isp_vacuum_s": float(isp_vac[i_ex]) / G0,
            "cea_theoretical_isp_exit_m_s": float(isp[i_ex]),
            "cea_theoretical_chamber_pressure_Pa": float(P[i_ch]),
            "cea_theoretical_chamber_temperature_K": float(T[i_ch]),
            "cea_theoretical_chamber_molecular_weight_kg_kmol": float(mw[i_ch]),
            "cea_theoretical_chamber_gamma_s": float(gamma[i_ch]),
            "cea_theoretical_exit_pressure_Pa": float(P[i_ex]),
            "cea_theoretical_exit_temperature_K": float(T[i_ex]),
            "cea_theoretical_exit_mach": float(mach[i_ex]),
            "cea_theoretical_exit_area_ratio": float(ae_at[i_ex]),
            "cea_theoretical_stations": int(rocket_solution.num_pts),
            "cea_theoretical_chamber_station_index": i_ch,
            "cea_theoretical_exit_station_index": i_ex,
            "cea_theoretical_exit_definition": exit_kind,
            "cea_theoretical_note": (
                "Theoretical CEA values for comparison only: GRIBS keeps its own "
                "unsteady internal-ballistics and nozzle model (cstar_eff_m_s, "
                "CF_eff, Isp_s are GRIBS results)."),
        }


# ------------------------------------------------------------------------------
# 3.5 Backend factory and solver-facing accessors
# ------------------------------------------------------------------------------
_THERMO: Optional[ThermochemistryBackend] = None


def build_backend(c: Config, reactants: Sequence[ReactantSpec],
                  packing_fraction: float, cache_root: Path) -> ThermochemistryBackend:
    """Instantiate the explicitly selected thermochemistry backend.

    There is deliberately **no default and no silent fallback**: the value comes
    from the mandatory ``thermochemistry.backend`` configuration key, and if the
    selected backend cannot initialize the error is raised with an actionable
    diagnostic instead of switching to another backend.
    """
    if c.thermo_backend == BACKEND_CEA_PYTHON:
        return CEAPythonBackend(c, reactants, packing_fraction, cache_root)
    if c.thermo_backend == BACKEND_CEA_LEGACY_EXECUTABLE:
        return CEALegacyExecutableBackend(c, reactants, packing_fraction, cache_root)
    if c.thermo_backend in REMOVED_BACKENDS:
        raise ConfigurationError(
            f"thermochemistry.backend: {c.thermo_backend!r} - "
            f"{REMOVED_BACKENDS[c.thermo_backend]}")
    raise ConfigurationError(
        f"thermochemistry.backend: unknown backend {c.thermo_backend!r}. "
        f"Set it explicitly to one of: {', '.join(KNOWN_BACKENDS)}.")


def initialize_thermochemistry(c: Config, default_cache_dir: Path) -> ThermochemistryBackend:
    """Build (or load from cache) the chamber-property backend and its density."""
    global _THERMO
    cache_root = Path(c.cea_cache_dir).expanduser() if c.cea_cache_dir else Path(default_cache_dir)
    _THERMO = build_backend(c, c.reactants_spec_list, c.packing_fraction, cache_root)
    # Additive-volume mixture density x packing fraction.  A CEA combustion-product
    # (gas) density is never used as the solid-propellant density.
    c.rho_p = _THERMO.bulk_mixture_density
    return _THERMO


def thermochemistry() -> ThermochemistryBackend:
    if _THERMO is None:
        raise RuntimeError("Thermochemistry has not been initialized.")
    return _THERMO


def gas_props(p0: float, c: Config) -> Tuple[float, float, float]:
    """Return chamber R [J/(kg K)], adiabatic equilibrium T0 [K], gamma_s [-]."""
    return thermochemistry().props(p0)


def theta_and_deriv(p0: float, c: Config) -> Tuple[float, float]:
    """Return Theta = R*T0 [J/kg] and dTheta/dp [J/(kg Pa)]."""
    return thermochemistry().theta_and_derivative(p0)


def is_extrapolated(p0: float, c: Config) -> bool:
    return thermochemistry().is_extrapolated(p0)


# ==============================================================================
# 4. BURN RATE (Saint-Robert + temperature sensitivity + erosive burning)
#    Unchanged from v0.3.1-alpha.
# ==============================================================================
def base_burn_rate(p0: float, c: Config) -> float:
    return (c.a_burn * (max(p0, 0.0) / c.p_ref) ** c.n_burn
            * math.exp(c.sigma_p * (c.T_grain - c.T_ref)))


def burn_rate(p0: float, c: Config, Ab: float, Ri: float) -> float:
    """Total regression rate including the optional erosive contribution."""
    r0 = base_burn_rate(p0, c)
    if c.ero_alpha <= 0.0 or Ri <= 0.0 or Ab <= 0.0:
        return r0
    A_port = math.pi * Ri ** 2
    D_h = 2.0 * Ri
    r = r0
    for _ in range(60):                     # damped fixed-point iteration
        G = c.rho_p * r * Ab / A_port       # port mass flux at the aft end
        if G <= 1.0e-12:
            return r0
        r_new = r0 + c.ero_alpha * G ** 0.8 / D_h ** 0.2 \
            * math.exp(-c.ero_beta * c.rho_p * r / G)
        if abs(r_new - r) <= 1.0e-13 * max(r, 1.0e-12):
            return r_new
        r = 0.5 * (r + r_new)
    return r


# ==============================================================================
# 5. NOZZLE MODEL  (unchanged from v0.3.1-alpha)
# ==============================================================================
def area_mach(M: float, g: float) -> float:
    """A/At as a function of Mach number (isentropic)."""
    return (1.0 / M) * ((2.0 / (g + 1.0))
                        * (1.0 + 0.5 * (g - 1.0) * M * M)) ** ((g + 1.0) / (2.0 * (g - 1.0)))


def mach_from_area(eps: float, g: float, supersonic: bool) -> float:
    """Invert the area-Mach relation."""
    if eps <= 1.0 + 1.0e-14:
        return 1.0
    f = lambda M: area_mach(M, g) - eps
    if supersonic:
        return brentq(f, 1.0 + 1.0e-12, 200.0, xtol=1e-14, rtol=8.9e-16, maxiter=200)
    return brentq(f, 1.0e-7, 1.0, xtol=1e-16, rtol=8.9e-16, maxiter=200)


def critical_ratio(g: float) -> float:
    """pe/p0 at the throat when choked."""
    return (2.0 / (g + 1.0)) ** (g / (g - 1.0))


def choke_threshold_pressure(p0: float, At: float, Ae: float, c: Config) -> float:
    """Chamber pressure above which the nozzle is choked (evaluated at gamma(p0))."""
    _, _, g = gas_props(p0, c)
    eps = Ae / At
    if eps <= 1.0 + 1.0e-12:
        return c.p_a / critical_ratio(g)
    Me_sub = mach_from_area(eps, g, supersonic=False)
    return c.p_a * (1.0 + 0.5 * (g - 1.0) * Me_sub ** 2) ** (g / (g - 1.0))


def choke_margin(p0: float, At: float, Ae: float, c: Config) -> float:
    """> 0 when choked, < 0 when subsonic."""
    return p0 - choke_threshold_pressure(p0, At, Ae, c)


def choke_limit_pressure(c: Config, At: float, Ae: float) -> Optional[float]:
    """Solve p0 = p_threshold(p0) (the fixed point marking the choking limit)."""
    f = lambda p: choke_margin(p, At, Ae, c)
    lo, hi = c.p_a * (1.0 + 1e-9), c.p_a * 1.0e4
    try:
        if f(lo) * f(hi) > 0.0:
            return None
        return brentq(f, lo, hi, xtol=1e-8, rtol=1e-14)
    except (ValueError, RuntimeError):
        return None


def nozzle_state(p0: float, At: float, Ae: float, c: Config) -> dict:
    """Mass flow, thrust, exit conditions and flow regime."""
    if p0 <= c.p_a:
        return dict(mdot=0.0, F=0.0, Me=0.0, pe=p0, ve=0.0,
                    regime="no-flow", A_eff=Ae)

    R, T0, g = gas_props(p0, c)
    eps = Ae / At
    sq = math.sqrt(g / (R * T0))
    p_thr = choke_threshold_pressure(p0, At, Ae, c)

    if p0 >= p_thr:                                   # ---- choked ----
        mdot = c.Cd * At * p0 * sq \
            * (2.0 / (g + 1.0)) ** ((g + 1.0) / (2.0 * (g - 1.0)))
        Me = 1.0 if eps <= 1.0 + 1.0e-12 else mach_from_area(eps, g, True)
        pe = p0 * (1.0 + 0.5 * (g - 1.0) * Me * Me) ** (-g / (g - 1.0))
        Te = T0 / (1.0 + 0.5 * (g - 1.0) * Me * Me)
        ve = Me * math.sqrt(g * R * Te)
        A_eff, pe_eff, ve_eff, regime = Ae, pe, ve, "choked"

        if c.use_separation and eps > 1.0 and pe < c.sep_ratio * c.p_a:
            p_sep = c.sep_ratio * c.p_a
            arg = (p0 / p_sep) ** ((g - 1.0) / g) - 1.0
            if arg > 0.0:
                M_sep = math.sqrt(2.0 / (g - 1.0) * arg)
                eps_sep = area_mach(M_sep, g)
                if 1.0 < eps_sep < eps:
                    Te_s = T0 / (1.0 + 0.5 * (g - 1.0) * M_sep ** 2)
                    A_eff = eps_sep * At
                    pe_eff = p_sep
                    ve_eff = M_sep * math.sqrt(g * R * Te_s)
                    Me = M_sep
                    regime = "separated"
        F = c.eta_thrust * mdot * ve_eff + (pe_eff - c.p_a) * A_eff
        return dict(mdot=mdot, F=F, Me=Me, pe=pe_eff, ve=ve_eff,
                    regime=regime, A_eff=A_eff)

    # ---- subsonic (exit pressure equals ambient); continuous with the above ----
    Me = math.sqrt(2.0 / (g - 1.0) * ((p0 / c.p_a) ** ((g - 1.0) / g) - 1.0))
    mdot = c.Cd * Ae * p0 * sq * Me \
        * (1.0 + 0.5 * (g - 1.0) * Me * Me) ** (-(g + 1.0) / (2.0 * (g - 1.0)))
    Te = T0 / (1.0 + 0.5 * (g - 1.0) * Me * Me)
    ve = Me * math.sqrt(g * R * Te)
    F = c.eta_thrust * mdot * ve
    return dict(mdot=mdot, F=F, Me=Me, pe=c.p_a, ve=ve,
                regime="subsonic", A_eff=Ae)


# ==============================================================================
# 6. QUASI-STEADY EQUILIBRIUM PRESSURE (diagnostic only)
# ==============================================================================
def equilibrium_pressure(Ab: float, Ri: float, At: float, Ae: float,
                         c: Config) -> Optional[float]:
    """Lowest pressure satisfying mdot_gen(p) = mdot_out(p); None if not found."""
    if not np.isfinite(Ab) or Ab <= 0.0:
        return None

    def bal(p):
        try:
            v = c.rho_p * burn_rate(p, c, Ab, Ri) * Ab - nozzle_state(p, At, Ae, c)["mdot"]
            return float(v) if np.isfinite(v) else np.nan
        except (ValueError, OverflowError, FloatingPointError, GribsError):
            return np.nan

    grid = np.geomspace(c.p_a * (1.0 + 1e-9), 1.0e11, 140)
    vals = np.array([bal(p) for p in grid])
    for p1, p2, f1, f2 in zip(grid[:-1], grid[1:], vals[:-1], vals[1:]):
        if not (np.isfinite(f1) and np.isfinite(f2)):
            continue
        if f1 == 0.0:
            return float(p1)
        if f1 * f2 < 0.0:
            return float(brentq(bal, p1, p2, xtol=1.0, rtol=1e-12))
    return None


# ==============================================================================
# 7. ODE SYSTEM
#    y = [p0, x, Rt, m_out, Impulse, m_gen]     (unchanged from v0.3.1-alpha)
# ==============================================================================
IP, IX, IRT, IMO, IIMP, IMG = range(6)


def igniter_mdot(t: float, c: Config) -> float:
    return c.ign_mdot if (c.ign_time > 0.0 and t <= c.ign_time) else 0.0


def make_rhs(c: Config, Ae: float, burning: bool) -> Callable:
    xw = web_thickness(c)
    p_floor = 500.0

    def rhs(t, y):
        p0 = max(float(y[IP]), p_floor)
        x = min(max(float(y[IX]), 0.0), xw)
        Rt = max(float(y[IRT]), 1.0e-9)
        At = math.pi * Rt * Rt

        Ri, Lp, Ab, Vg = geometry(x, c)
        if burning and Ab > 0.0:
            r = burn_rate(p0, c, Ab, Ri)
        else:
            r, Ab = 0.0, 0.0
        m_gen = c.rho_p * r * Ab + igniter_mdot(t, c)

        nz = nozzle_state(p0, At, Ae, c)
        theta, dtheta = theta_and_deriv(p0, c)
        denom = Vg * (1.0 - p0 * dtheta / theta)
        if denom <= 0.0:
            raise ValueError("Non-physical factor (1 - p0*Theta'/Theta) <= 0.")

        dp = (theta * (m_gen - nz["mdot"]) - p0 * Ab * r) / denom
        dRt = c.ero_throat_c * (p0 / c.p_ref) ** c.ero_throat_m
        return np.array([dp, r, dRt, nz["mdot"], nz["F"], m_gen])

    return rhs


def _event(fn, terminal=True, direction=0.0):
    fn.terminal = terminal
    fn.direction = direction
    return fn


def run_model(c: Config) -> dict:
    """Integrate the burning phase and (optionally) the blowdown phase."""
    validate(c)
    At0 = math.pi * c.R_t0 ** 2
    Ae = c.eps_nozzle * At0                      # geometric exit area, fixed
    xw = web_thickness(c)

    init_regime = nozzle_state(c.p0_init, At0, Ae, c)["regime"]
    if c.unchoked_policy == "stop" and init_regime != "choked":
        raise ValueError("Initial state is not choked: policy 'stop' has "
                         "nothing to integrate.")

    theta0, _ = theta_and_deriv(c.p0_init, c)
    y0 = np.array([c.p0_init, 0.0, c.R_t0, 0.0, 0.0, 0.0])
    m_gas0 = c.p0_init * c.V_g0 / theta0

    atol = np.array([c.atol_p, c.atol_x, 1e-14, 1e-12, 1e-9, 1e-12])

    ev_burn = _event(lambda t, y: y[IX] - xw, True, +1.0)
    ev_unch = _event(lambda t, y: choke_margin(max(y[IP], 1.0),
                                               math.pi * max(y[IRT], 1e-9) ** 2,
                                               Ae, c),
                     c.unchoked_policy == "stop", -1.0)
    ev_rech = _event(lambda t, y: choke_margin(max(y[IP], 1.0),
                                               math.pi * max(y[IRT], 1e-9) ** 2,
                                               Ae, c), False, +1.0)
    ev_amb = _event(lambda t, y: y[IP] - c.p_a * 1.0005, True, -1.0)
    names = ["burnout", "unchoked", "rechoked", "ambient"]

    sol1 = solve_ivp(make_rhs(c, Ae, True), (0.0, c.t_max), y0,
                     method=c.method, events=[ev_burn, ev_unch, ev_rech, ev_amb],
                     dense_output=True, rtol=c.rtol, atol=atol,
                     max_step=c.max_step_burn, first_step=1.0e-7)
    if not sol1.success:
        raise RuntimeError(f"Burning-phase integration failed: {sol1.message}")

    stop1 = "time-limit"
    for nm, te in zip(names, sol1.t_events):
        if nm in ("burnout", "ambient") and len(te):
            stop1 = nm
            break
        if nm == "unchoked" and c.unchoked_policy == "stop" and len(te):
            stop1 = "choking-loss"
            break

    sol2, stop2 = None, "not-run"
    if c.blowdown and stop1 == "burnout":
        y_bo = sol1.y[:, -1].copy()
        y_bo[IX] = xw
        ev_amb2 = _event(lambda t, y: y[IP] - c.p_a * 1.0005, True, -1.0)
        ev_unch2 = _event(lambda t, y: choke_margin(max(y[IP], 1.0),
                                                    math.pi * max(y[IRT], 1e-9) ** 2,
                                                    Ae, c),
                          c.unchoked_policy == "stop", -1.0)
        sol2 = solve_ivp(make_rhs(c, Ae, False),
                         (sol1.t[-1], sol1.t[-1] + c.blowdown_tmax), y_bo,
                         method=c.method, events=[ev_amb2, ev_unch2],
                         dense_output=True, rtol=c.rtol, atol=atol,
                         max_step=c.max_step_blow)
        if not sol2.success:
            raise RuntimeError(f"Blowdown integration failed: {sol2.message}")
        stop2 = "ambient" if len(sol2.t_events[0]) else "time-limit"
        if c.unchoked_policy == "stop" and len(sol2.t_events[1]):
            stop2 = "choking-loss"

    return dict(sol_burn=sol1, sol_blow=sol2, stop_burn=stop1, stop_blow=stop2,
                Ae=Ae, At0=At0, m_gas0=m_gas0, init_regime=init_regime,
                t_events={nm: [float(v) for v in te]
                          for nm, te in zip(names, sol1.t_events)})


def validate(c: Config) -> None:
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


# ==============================================================================
# 8. POST-PROCESSING  (unchanged from v0.3.1-alpha)
# ==============================================================================
def _phase_grid(sol, n_lin: int, log_head: bool) -> np.ndarray:
    t0, t1 = float(sol.t[0]), float(sol.t[-1])
    if t1 <= t0:
        return np.array([t0])
    g = [np.linspace(t0, t1, n_lin)]
    span = t1 - t0
    if log_head:
        g.append(t0 + np.geomspace(max(1e-7, span * 1e-7), span, 900))
    g.append(np.asarray(sol.t, dtype=float))
    for te in sol.t_events:
        if len(te):
            g.append(np.asarray(te, dtype=float))
    t = np.unique(np.concatenate(g))
    return t[(t >= t0) & (t <= t1)]


def sample(res: dict, c: Config) -> dict:
    Ae = res["Ae"]
    xw = web_thickness(c)
    segs = [(res["sol_burn"], "burn", 4000, True)]
    if res["sol_blow"] is not None:
        segs.append((res["sol_blow"], "blowdown", 1600, True))

    rec = {k: [] for k in
           ("t", "phase", "p0", "x", "Rt", "At", "Kn", "Ri", "Lp", "Ab", "Vg",
            "r", "mdot_gen", "mdot_out", "F", "Me", "pe", "ve", "regime",
            "R", "T0", "gamma", "m_out", "impulse", "m_gen", "m_gas_eos",
            "m_gas_bal", "extrap", "cstar", "CF")}

    first = True
    for sol, phase, n_lin, head in segs:
        tg = _phase_grid(sol, n_lin, head)
        if not first:
            tg = tg[1:]
        first = False
        ys = sol.sol(tg)
        for k, t in enumerate(tg):
            p0 = max(float(ys[IP, k]), 1.0)
            x = min(max(float(ys[IX, k]), 0.0), xw)
            Rt = float(ys[IRT, k])
            At = math.pi * Rt * Rt
            Ri, Lp, Ab, Vg = geometry(x, c)
            R, T0, g = gas_props(p0, c)
            nz = nozzle_state(p0, At, Ae, c)
            if phase == "burn":
                r = burn_rate(p0, c, Ab, Ri)
                mg = c.rho_p * r * Ab + igniter_mdot(t, c)
            else:
                r, Ab, mg = 0.0, 0.0, 0.0
            theta = R * T0
            m_eos = p0 * Vg / theta
            m_bal = res["m_gas0"] + float(ys[IMG, k]) - float(ys[IMO, k])

            rec["t"].append(float(t));          rec["phase"].append(phase)
            rec["p0"].append(p0);               rec["x"].append(x)
            rec["Rt"].append(Rt);               rec["At"].append(At)
            rec["Kn"].append(Ab / At);          rec["Ri"].append(Ri)
            rec["Lp"].append(Lp);               rec["Ab"].append(Ab)
            rec["Vg"].append(Vg);               rec["r"].append(r)
            rec["mdot_gen"].append(mg);         rec["mdot_out"].append(nz["mdot"])
            rec["F"].append(nz["F"]);           rec["Me"].append(nz["Me"])
            rec["pe"].append(nz["pe"]);         rec["ve"].append(nz["ve"])
            rec["regime"].append(nz["regime"]); rec["R"].append(R)
            rec["T0"].append(T0);               rec["gamma"].append(g)
            rec["m_out"].append(float(ys[IMO, k]))
            rec["impulse"].append(float(ys[IIMP, k]))
            rec["m_gen"].append(float(ys[IMG, k]))
            rec["m_gas_eos"].append(m_eos);     rec["m_gas_bal"].append(m_bal)
            rec["extrap"].append(is_extrapolated(p0, c))
            rec["cstar"].append(p0 * At / nz["mdot"] if nz["mdot"] > 0 else np.nan)
            rec["CF"].append(nz["F"] / (p0 * At) if p0 * At > 0 else np.nan)

    out = {k: (np.asarray(v) if k not in ("phase", "regime") else np.asarray(v, dtype=object))
           for k, v in rec.items()}

    # quasi-steady equilibrium pressure (diagnostic, burning phase only)
    t = out["t"]
    p_eq = np.full_like(t, np.nan)
    burn_idx = np.nonzero(out["phase"] == "burn")[0]
    if burn_idx.size:
        pick = burn_idx[np.unique(np.linspace(0, burn_idx.size - 1,
                                              min(200, burn_idx.size)).astype(int))]
        vals = np.array([equilibrium_pressure(out["Ab"][i], out["Ri"][i],
                                              out["At"][i], Ae, c) or np.nan
                         for i in pick])
        ok = np.isfinite(vals)
        if ok.sum() >= 2:
            p_eq[burn_idx] = np.interp(t[burn_idx], t[pick][ok], vals[ok],
                                       left=np.nan, right=np.nan)
    out["p_eq"] = p_eq
    return out


def summarize(res: dict, h: dict, c: Config) -> dict:
    t = h["t"]
    burn = h["phase"] == "burn"
    m_prop = propellant_mass(c)
    At0, Ae = res["At0"], res["Ae"]
    t_bo = float(t[burn][-1])

    i_pmax = int(np.argmax(h["p0"]))
    i_fmax = int(np.argmax(h["F"]))
    I_burn = float(h["impulse"][burn][-1])
    I_tot = float(h["impulse"][-1])
    m_out_tot = float(h["m_out"][-1])
    m_gen_tot = float(h["m_gen"][burn][-1])

    pAt_int = float(_TRAPZ(h["p0"][burn] * h["At"][burn], t[burn]))
    cstar_eff = pAt_int / m_out_tot if m_out_tot > 0 else float("nan")
    CF_eff = I_burn / pAt_int if pAt_int > 0 else float("nan")
    p_mean = float(pAt_int / (At0 * (t_bo - t[0]))) if t_bo > t[0] else float("nan")

    denom = np.maximum(np.abs(h["m_gas_eos"]), 1e-15)
    mass_err = float(np.max(np.abs(h["m_gas_bal"] - h["m_gas_eos"]) / denom))

    p_star = choke_limit_pressure(c, At0, Ae)
    mg0 = c.rho_p * base_burn_rate(c.p0_init, c) * geometry(0.0, c)[2]
    mo0 = nozzle_state(c.p0_init, At0, Ae, c)["mdot"]

    backend = thermochemistry()
    s = dict(
        web_mm=web_thickness(c) * 1e3,
        propellant_density_kg_m3=c.rho_p,
        propellant_mass_g=m_prop * 1e3,
        burn_time_s=t_bo,
        total_time_s=float(t[-1]),
        burn_stop_reason=res["stop_burn"],
        blowdown_stop_reason=res["stop_blow"],
        initial_regime=res["init_regime"],
        initial_mdot_gen_kg_s=mg0,
        initial_mdot_out_kg_s=mo0,
        choke_limit_pressure_Pa=p_star,
        p_max_Pa=float(h["p0"][i_pmax]), t_p_max_s=float(t[i_pmax]),
        p_min_burn_Pa=float(np.min(h["p0"][burn])),
        t_p_min_s=float(t[burn][int(np.argmin(h["p0"][burn]))]),
        p_mean_burn_Pa=p_mean,
        F_max_N=float(h["F"][i_fmax]), t_F_max_s=float(t[i_fmax]),
        F_mean_burn_N=float(I_burn / (t_bo - t[0])) if t_bo > t[0] else float("nan"),
        impulse_burn_Ns=I_burn,
        impulse_total_Ns=I_tot,
        impulse_blowdown_Ns=I_tot - I_burn,
        Isp_s=I_tot / (m_prop * G0) if m_prop > 0 else float("nan"),
        cstar_eff_m_s=cstar_eff,
        CF_eff=CF_eff,
        Kn_initial=float(h["Kn"][0]), Kn_max=float(np.max(h["Kn"][burn])),
        Kn_final=float(h["Kn"][burn][-1]),
        throat_radius_final_mm=float(h["Rt"][-1]) * 1e3,
        propellant_mass_balance_error=abs(m_gen_tot / m_prop - 1.0) if m_prop > 0 else float("nan"),
        gas_mass_consistency_error=mass_err,
        property_extrapolation=bool(np.any(h["extrap"])),
        property_extrapolation_points=int(np.count_nonzero(h["extrap"])),
        property_extrapolation_fraction=float(np.mean(h["extrap"])),
        property_policy=c.property_policy,
        thermochemistry_backend=backend.name,
        property_range_min_Pa=float(backend.pmin),
        property_range_max_Pa=float(backend.pmax),
        pressure_min_Pa=float(np.min(h["p0"])),
        pressure_max_Pa=float(np.max(h["p0"])),
        regimes_visited=sorted(set(h["regime"].tolist())),
        event_times=res["t_events"],
    )
    return s


# ==============================================================================
# 9. SELF-TESTS
# ==============================================================================
#: Documented tolerance of the analytic-vs-finite-difference derivative check.
#: The chamber-property tables use a C1-continuous PCHIP spline whose second
#: derivative jumps at the knots, so a centred difference with h = 1e-5 p is
#: limited to O(h) there (1e-5 relative tolerance).
DERIVATIVE_TOLERANCE_TABLE = 1.0e-5


def self_tests(c: Config, verbose=True) -> bool:
    out, ok_all = [], True

    def chk(name, ok, detail=""):
        nonlocal ok_all
        ok_all &= bool(ok)
        out.append(f"  [{'PASS' if ok else 'FAIL'}] {name:<46s} {detail}")

    backend = thermochemistry()

    # T1 analytic dTheta/dp against a centred difference (derivative validation)
    tol_deriv = DERIVATIVE_TOLERANCE_TABLE
    e = 0.0
    for p in np.geomspace(max(c.p_fit_min * 1.05, 1.1e5), c.p_fit_max * 0.95, 25):
        hstep = p * 1e-5
        fd = (theta_and_deriv(p + hstep, c)[0] - theta_and_deriv(p - hstep, c)[0]) / (2 * hstep)
        e = max(e, abs(fd / theta_and_deriv(p, c)[1] - 1.0))
    chk("analytic vs numerical dTheta/dp", e < tol_deriv,
        f"max rel. err = {e:.2e} (tol {tol_deriv:.0e})")

    # T2 geometric consistency dVg/dx = Ab
    xw = web_thickness(c)
    e = 0.0
    for x in np.linspace(0.02 * xw, 0.98 * xw, 25):
        hstep = xw * 1e-6
        fd = (geometry(x + hstep, c)[3] - geometry(x - hstep, c)[3]) / (2 * hstep)
        e = max(e, abs(fd / geometry(x, c)[2] - 1.0))
    chk("geometry consistency dVg/dx = Ab", e < 1e-7, f"max rel. err = {e:.2e}")

    # T3 area-Mach inversion round trip
    e = 0.0
    for eps in (1.5, 2.5, 6.0, 20.0):
        for g in (1.10, 1.15, 1.25):
            e = max(e, abs(area_mach(mach_from_area(eps, g, True), g) / eps - 1.0),
                    abs(area_mach(mach_from_area(eps, g, False), g) / eps - 1.0))
    chk("area-Mach inversion round trip", e < 1e-10, f"max rel. err = {e:.2e}")

    # T4 continuity of mdot and F across the choking boundary
    # Test only the nozzle geometry specified in the current configuration.
    worst_mdot = 0.0
    worst_thrust = 0.0

    eps = c.eps_nozzle

    cc = Config(**{
        **asdict(c),
        "eps_nozzle": eps,
        "use_separation": False
    })

    At0 = math.pi * cc.R_t0 ** 2
    Ae = eps * At0

    pstar = choke_limit_pressure(cc, At0, Ae)

    if pstar is None:
        chk("choked/subsonic continuity (mdot, F)", False,
            "choking-limit pressure could not be determined")
    else:
        lo = nozzle_state(pstar * (1.0 - 1e-9), At0, Ae, cc)
        hi = nozzle_state(pstar * (1.0 + 1e-9), At0, Ae, cc)
        worst_mdot = abs(lo["mdot"] / hi["mdot"] - 1.0)
        worst_thrust = abs((lo["F"] + 1e-12) / (hi["F"] + 1e-12) - 1.0)
        worst = max(worst_mdot, worst_thrust)
        chk("choked/subsonic continuity (mdot, F)", worst < 1e-6,
            f"Ae/At = {eps:.6g}, mdot jump = {worst_mdot:.2e}, "
            f"F jump = {worst_thrust:.2e}")

        # T5 closed-form check of the converging-nozzle choked thrust
        cc = Config(**{**asdict(c), "eps_nozzle": 1.0, "Cd": 1.0,
                       "eta_thrust": 1.0, "use_separation": False})
        At0 = math.pi * cc.R_t0 ** 2
        e = 0.0
        for p in np.geomspace(cc.p_a / 0.5, 50e6, 20):
            nz = nozzle_state(p, At0, At0, cc)
            _, _, g = gas_props(p, cc)
            F_ref = At0 * ((g + 1.0) * critical_ratio(g) * p - cc.p_a)
            e = max(e, abs(nz["F"] / F_ref - 1.0))
        chk("closed-form choked thrust identity", e < 1e-12,
            f"max rel. err = {e:.2e}")

        # T6 finite and physical right-hand side at t = 0
        try:
            dy = make_rhs(c, c.eps_nozzle * math.pi * c.R_t0 ** 2, True)(
                0.0, np.array([c.p0_init, 0.0, c.R_t0, 0.0, 0.0, 0.0]))
            chk("finite initial derivatives", np.all(np.isfinite(dy)),
                f"dp/dt = {dy[0]:+.3e} Pa/s, r = {dy[1]:.3e} m/s")
        except Exception as exc:                              # pragma: no cover
            chk("finite initial derivatives", False, str(exc))

    # T7 chamber-property table completeness and physicality
    try:
        pv = backend.property_validation
        ok = bool(pv.get("all_points_finite", False)) and \
            bool(pv.get("pressure_strictly_increasing", False)) and \
            float(pv.get("temperature_min_K", 0.0)) > 0.0 and \
            float(pv.get("gamma_min", 0.0)) > 1.0 and \
            float(pv.get("min_pressure_factor_1_minus_pTheta_over_Theta", -1.0)) > 0.0
        detail = (f"{pv.get('pressure_points', 'n/a')} points, "
                  f"T0 = {pv.get('temperature_min_K', float('nan')):.1f}-"
                  f"{pv.get('temperature_max_K', float('nan')):.1f} K, "
                  f"gamma = {pv.get('gamma_min', float('nan')):.4f}-"
                  f"{pv.get('gamma_max', float('nan')):.4f}")
        chk("chamber-property table physical", ok, detail)
    except Exception as exc:                                  # pragma: no cover
        chk("chamber-property table physical", False, str(exc))

    if verbose:
        print(f"Self-tests  (thermochemistry backend: {backend.name})")
        print("\n".join(out))
        print(f"  -> {'all checks passed' if ok_all else 'FAILURES DETECTED'}\n")
    return ok_all


# ==============================================================================
# 10. PLOTTING (one single, carefully designed figure)
#     Panel definitions and meanings unchanged from v0.3.1-alpha.
# ==============================================================================
C_PRESS = "#14507d"
C_THRUST = "#b3331f"
C_GEN = "#2e8b57"
C_OUT = "#7b4fa8"
C_ACC = "#d98c00"
C_GREY = "#5a5a5a"
REGIME_FACE = {"choked": "#dce9f5", "subsonic": "#fdecd2",
               "separated": "#ece2f7", "no-flow": "#ededed"}
REGIME_NAME = {"choked": "choked (sonic throat)", "subsonic": "subsonic nozzle",
               "separated": "separated (overexpanded)", "no-flow": "no outflow"}


def _style() -> None:
    plt.rcParams.update({
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "axes.facecolor": "#fdfdfe",
        "axes.edgecolor": "#54606b",
        "axes.linewidth": 0.9,
        "axes.grid": True,
        "axes.axisbelow": True,
        "grid.color": "#c5ccd3",
        "grid.linewidth": 0.55,
        "grid.alpha": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "font.family": "DejaVu Sans",
        "font.size": 9.0,
        "axes.titlesize": 10.0,
        "axes.titleweight": "bold",
        "axes.labelsize": 9.0,
        "legend.fontsize": 7.8,
        "legend.frameon": True,
        "legend.framealpha": 0.92,
        "legend.edgecolor": "#b9c2ca",
        "xtick.labelsize": 8.2,
        "ytick.labelsize": 8.2,
        "xtick.direction": "out",
        "ytick.direction": "out",
    })


def _regime_spans(t, regime):
    spans, i0 = [], 0
    for i in range(1, len(t)):
        if regime[i] != regime[i0]:
            spans.append((t[i0], t[i], regime[i0]))
            i0 = i
    spans.append((t[i0], t[-1], regime[i0]))
    return spans


def _shade(ax, spans, t_lo=None, t_hi=None, scale=1.0):
    for a, b, key in spans:
        if t_lo is not None:
            a, b = max(a, t_lo), min(b, t_hi)
        if b <= a:
            continue
        ax.axvspan(a * scale, b * scale, color=REGIME_FACE.get(key, "#eeeeee"),
                   alpha=0.85, lw=0, zorder=0)


#: Style of the key-results block and the geometry of its automatic fit.
#: The block is sized to its own axes: it grows up to KEY_RESULTS_MAX_FONTSIZE for
#: readability and shrinks (with a warning) if a configuration produces so many
#: lines that it would otherwise run over the neighbouring panels.
KEY_RESULTS_MAX_FONTSIZE = 9.6      # never grow beyond this (readability)
KEY_RESULTS_MIN_FONTSIZE = 5.2      # never shrink below this (a warning is printed)
KEY_RESULTS_LINESPACING = 1.28      # multiple of the font size
KEY_RESULTS_BOX_PAD = 0.55          # bbox pad, in font-size units
KEY_RESULTS_CHAR_WIDTH = 0.602      # DejaVu Sans Mono advance width, in em


def key_results_fontsize(fig, ax, n_lines: int, longest_line: int) -> float:
    """Point size whose text box just fits inside ``ax`` (clamped to sane limits).

    The box size is estimated analytically (line count, monospace advance width
    and the ``bbox`` padding) from the axes geometry, so the fit is decided
    before anything is drawn and the key-results block can never spill over the
    neighbouring panels - whatever the configuration produces (extra event
    lines, CEA diagnostics, long backend names, ...).  A negative return value
    means "does not fit even at the minimum size".
    """
    width_pt = fig.get_size_inches()[0] * ax.get_position().width * 72.0
    height_pt = fig.get_size_inches()[1] * ax.get_position().height * 72.0
    fs_h = height_pt / (max(1, n_lines) * KEY_RESULTS_LINESPACING + 2.0 * KEY_RESULTS_BOX_PAD + 0.8)
    fs_w = width_pt / (max(1, longest_line) * KEY_RESULTS_CHAR_WIDTH + 2.0 * KEY_RESULTS_BOX_PAD)
    return float(min(KEY_RESULTS_MAX_FONTSIZE, fs_h, fs_w))


def _draw_key_results(fig, ax, lines: Sequence[str]):
    """Draw the key-results block, fitted to (and contained in) its own axes."""
    text = "\n".join(lines)
    fontsize = key_results_fontsize(fig, ax, len(lines), max(len(l) for l in lines))
    if fontsize < KEY_RESULTS_MIN_FONTSIZE:
        print(f"WARNING: the key-results block needs a font size of "
              f"{fontsize:.2f} pt to fit its panel; clamping to "
              f"{KEY_RESULTS_MIN_FONTSIZE:.2f} pt (the block may overflow).")
        fontsize = KEY_RESULTS_MIN_FONTSIZE
    return ax.text(0.0, 1.0, text, transform=ax.transAxes, va="top", ha="left",
                   family="monospace", fontsize=fontsize,
                   linespacing=KEY_RESULTS_LINESPACING,
                   bbox=dict(boxstyle=f"round,pad={KEY_RESULTS_BOX_PAD}",
                             fc="#f4f7fa", ec="#a9bbcc", lw=1.0))


def make_figure(res, h, s, c: Config, path: Path, cea_diag: Optional[dict] = None) -> None:
    _style()
    t = h["t"]
    burn = h["phase"] == "burn"
    t_bo = s["burn_time_s"]
    spans = _regime_spans(t, h["regime"])
    has_blow = res["sol_blow"] is not None

    # The key-results box needs much more vertical room than the tallest plot
    # panel, so it gets its own full-height column on the right-hand side
    # (gs[:, 2]) instead of sharing a cell with a plot panel.  With the previous
    # 3x3 layout the text ran over the panel below it.
    fig = plt.figure(figsize=(16.6, 15.0))
    gs = GridSpec(4, 3, figure=fig, height_ratios=[1.45, 1.0, 1.0, 1.0],
                  hspace=0.40, wspace=0.28,
                  left=0.055, right=0.975, top=0.940, bottom=0.072)

    ax_p = fig.add_subplot(gs[0, 0:2])       # (a) chamber pressure (wide)
    ax_key = fig.add_subplot(gs[:, 2])       # key results (full-height column)
    ax_F = fig.add_subplot(gs[1, 0])         # (b) thrust
    ax_md = fig.add_subplot(gs[2, 0])        # (c) mass-flow balance
    ax_geo = fig.add_subplot(gs[3, 0])       # (d) grain geometry evolution
    ax_zm = fig.add_subplot(gs[1, 1])        # (e) ignition transient
    ax_nz = fig.add_subplot(gs[2, 1])        # (f) nozzle operating point
    ax_r = fig.add_subplot(gs[3, 1])         # (g) regression of the grain

    # ---------------- (a) chamber pressure ----------------
    _shade(ax_p, spans)
    ax_p.fill_between(t[burn], 0.0, h["p0"][burn] / 1e6, color=C_PRESS, alpha=0.13, lw=0)
    ax_p.plot(t[burn], h["p0"][burn] / 1e6, color=C_PRESS, lw=2.1, label="chamber pressure $p_0$")
    if has_blow:
        ax_p.plot(t[~burn], h["p0"][~burn] / 1e6, color=C_PRESS, lw=1.3, ls=":",
                  label="blowdown (post-burnout)")
    if np.isfinite(h["p_eq"]).any():
        mk = np.isfinite(h["p_eq"])
        ax_p.plot(t[mk], h["p_eq"][mk] / 1e6, color=C_GREY, lw=1.1, ls="--",
                  label="quasi-steady equilibrium")
    if s["choke_limit_pressure_Pa"]:
        ax_p.axhline(s["choke_limit_pressure_Pa"] / 1e6, color="#7a7a7a", lw=0.9, ls="-.")
        ax_p.text(0.995, s["choke_limit_pressure_Pa"] / 1e6, "  choking limit",
                  transform=ax_p.get_yaxis_transform(), ha="right", va="bottom",
                  fontsize=7.4, color="#5a5a5a")
    if s["property_extrapolation"]:
        ex = h["extrap"]
        ax_p.plot(t[ex], h["p0"][ex] / 1e6, color="#c0392b", lw=2.6, alpha=0.55,
                  solid_capstyle="butt", label="outside property-table range")
    ax_p.axvline(t_bo, color="#404040", lw=1.0, ls="--")
    ax_p.annotate("burnout", xy=(t_bo, ax_p.get_ylim()[1]),
                  xytext=(-4, -10), textcoords="offset points",
                  ha="right", va="top", fontsize=8, color="#404040", rotation=90)
    ax_p.plot([s["t_p_max_s"]], [s["p_max_Pa"] / 1e6], "o", ms=5.0,
              mfc="white", mec=C_PRESS, mew=1.6, zorder=5)
    ax_p.annotate(f"MEOP {s['p_max_Pa']/1e6:.3f} MPa @ {s['t_p_max_s']:.3f} s",
                  xy=(s["t_p_max_s"], s["p_max_Pa"] / 1e6),
                  xytext=(-14, 14), textcoords="offset points", ha="right",
                  fontsize=8.2, color=C_PRESS,
                  arrowprops=dict(arrowstyle="->", color=C_PRESS, lw=0.9))
    ax_p.set_xlim(0.0, t[-1])
    ax_p.set_ylim(0.0, max(h["p0"]) / 1e6 * 1.22)
    ax_p.set_xlabel("time  $t$  [s]")
    ax_p.set_ylabel("chamber pressure  $p_0$  [MPa]")
    ax_p.set_title("(a)  Chamber pressure history", loc="left")
    ax_p.legend(loc="lower right", ncol=1)

    # ---------------- key results panel ----------------
    ax_key.axis("off")
    ev = s["event_times"]
    lines = [
        "KEY RESULTS",
        "-" * 42,
        f"{'web thickness':<24}{s['web_mm']:>11.3f} mm",
        f"{'propellant density':<24}{s['propellant_density_kg_m3']:>11.2f} kg/m3",
        f"{'propellant mass':<24}{s['propellant_mass_g']:>11.3f} g",
        f"{'web burn time':<24}{s['burn_time_s']:>11.4f} s",
        f"{'total time computed':<24}{s['total_time_s']:>11.4f} s",
        "",
        f"{'max chamber pressure':<24}{s['p_max_Pa']/1e6:>11.4f} MPa",
        f"{'min chamber pressure':<24}{s['p_min_burn_Pa']/1e6:>11.4f} MPa",
        f"{'mean chamber pressure':<24}{s['p_mean_burn_Pa']/1e6:>11.4f} MPa",
        f"{'max thrust':<24}{s['F_max_N']:>11.3f} N",
        f"{'mean thrust (burn)':<24}{s['F_mean_burn_N']:>11.3f} N",
        "",
        f"{'total impulse':<24}{s['impulse_total_Ns']:>11.3f} N s",
        f"{'  of which blowdown':<24}{s['impulse_blowdown_Ns']:>11.3f} N s",
        f"{'effective Isp':<24}{s['Isp_s']:>11.2f} s",
        f"{'effective c*':<24}{s['cstar_eff_m_s']:>11.1f} m/s",
        f"{'effective CF':<24}{s['CF_eff']:>11.4f} -",
        "",
        f"{'Kn  (start / max / end)':<24}",
        f"{'':<10}{s['Kn_initial']:>8.1f} /{s['Kn_max']:>8.1f} /{s['Kn_final']:>8.1f}",
        f"{'nozzle regimes':<24}{', '.join(s['regimes_visited']):>11s}",
        "",
        "THERMOCHEMISTRY",
        f"  backend            {s['thermochemistry_backend']}",
        f"  property range     {s['property_range_min_Pa']/1e5:.3g}"
        f"-{s['property_range_max_Pa']/1e5:.3g} bar",
        "",
        "VERIFICATION",
        f"{'propellant mass balance':<24}{s['propellant_mass_balance_error']:>11.2e}",
        f"{'gas mass consistency':<24}{s['gas_mass_consistency_error']:>11.2e}",
        f"{'property extrapolation':<24}{str(s['property_extrapolation']):>11s}",
    ]
    if cea_diag:
        lines += [
            "",
            "CEA theoretical (comparison only)",
            f"  {'c*':<18}{cea_diag['cea_theoretical_cstar_m_s']:>11.1f} m/s",
            f"  {'CF':<18}{cea_diag['cea_theoretical_cf']:>11.4f} -",
            f"  {'Isp':<18}{cea_diag['cea_theoretical_isp_s']:>11.2f} s",
            f"  {'GRIBS c* (effective)':<18}{s['cstar_eff_m_s']:>11.1f} m/s",
        ]
    if ev.get("unchoked"):
        lines.append(f"{'first unchoking at':<24}{ev['unchoked'][0]*1e3:>11.3f} ms")
    if ev.get("rechoked"):
        lines.append(f"{'re-choking at':<24}{ev['rechoked'][0]*1e3:>11.3f} ms")
    _draw_key_results(fig, ax_key, lines)

    # ---------------- (b) thrust ----------------
    _shade(ax_F, spans)
    ax_F.fill_between(t[burn], 0.0, h["F"][burn], color=C_THRUST, alpha=0.15, lw=0,
                      label=f"impulse = {s['impulse_burn_Ns']:.2f} N s")
    ax_F.plot(t[burn], h["F"][burn], color=C_THRUST, lw=2.0, label="thrust $F$")
    if has_blow:
        ax_F.plot(t[~burn], h["F"][~burn], color=C_THRUST, lw=1.2, ls=":",
                  label="blowdown")
    ax_F.axvline(t_bo, color="#404040", lw=1.0, ls="--")
    ax_F.plot([s["t_F_max_s"]], [s["F_max_N"]], "o", ms=5, mfc="white",
              mec=C_THRUST, mew=1.6, zorder=5)
    ax_F.annotate(f"{s['F_max_N']:.2f} N", xy=(s["t_F_max_s"], s["F_max_N"]),
                  xytext=(-12, 10), textcoords="offset points", ha="right",
                  fontsize=8.2, color=C_THRUST,
                  arrowprops=dict(arrowstyle="->", color=C_THRUST, lw=0.9))
    ax_F.set_xlim(0.0, t[-1])
    ax_F.set_ylim(0.0, max(h["F"]) * 1.22)
    ax_F.set_xlabel("time  $t$  [s]")
    ax_F.set_ylabel("thrust  $F$  [N]")
    ax_F.set_title("(b)  Thrust history", loc="left")
    ax_F.legend(loc="upper left")

    # ---------------- (c) mass flow balance ----------------
    _shade(ax_md, spans)
    pos = h["mdot_gen"] > 0
    ax_md.plot(t[pos], h["mdot_gen"][pos] * 1e3, color=C_GEN, lw=1.8,
               label=r"generated  $\dot{m}_{gen}=\rho_p A_b r$")
    ax_md.plot(t, np.maximum(h["mdot_out"], 1e-12) * 1e3, color=C_OUT, lw=1.6,
               ls="--", label=r"nozzle  $\dot{m}_{out}$")
    acc = h["mdot_gen"] > h["mdot_out"]
    ax_md.fill_between(t, np.maximum(h["mdot_out"], 1e-12) * 1e3,
                       np.maximum(h["mdot_gen"], 1e-12) * 1e3, where=acc,
                       color=C_ACC, alpha=0.25, lw=0, label="chamber filling")
    ax_md.axvline(t_bo, color="#404040", lw=1.0, ls="--")
    ax_md.set_yscale("log")
    ax_md.set_xlim(0.0, t[-1])
    ax_md.set_xlabel("time  $t$  [s]")
    ax_md.set_ylabel(r"mass flow  [g/s]")
    ax_md.set_title("(c)  Mass-flow balance", loc="left")
    ax_md.legend(loc="lower right")

    # ---------------- (d) burning area / Kn / free volume ----------------
    ax_geo.plot(t[burn], h["Ab"][burn] * 1e4, color="#1b7f79", lw=1.9,
                label="burning area $A_b$")
    ax_geo.set_xlabel("time  $t$  [s]")
    ax_geo.set_ylabel("$A_b$  [cm$^2$]", color="#1b7f79")
    ax_geo.tick_params(axis="y", colors="#1b7f79")
    ax_geo.set_xlim(0.0, t[-1])
    ax_geo.axvline(t_bo, color="#404040", lw=1.0, ls="--")
    ax_g2 = ax_geo.twinx()
    ax_g2.grid(False)
    ax_g2.plot(t[burn], h["Kn"][burn], color="#8a5a00", lw=1.5, ls="--",
               label="$K_n=A_b/A_t$")
    ax_g2.plot(t[burn], h["Vg"][burn] * 1e6, color="#7a7a7a", lw=1.2, ls=":",
               label="free volume $V_g$ [cm$^3$]")
    ax_g2.set_ylabel("$K_n$ [-]   /   $V_g$ [cm$^3$]", color="#8a5a00")
    ax_g2.tick_params(axis="y", colors="#8a5a00")
    ax_g2.spines["right"].set_visible(True)
    hl = ax_geo.get_legend_handles_labels()
    h2 = ax_g2.get_legend_handles_labels()
    ax_geo.legend(hl[0] + h2[0], hl[1] + h2[1], loc="center right")
    ax_geo.set_title("(d)  Grain geometry evolution", loc="left")

    # ---------------- (e) start-up transient ----------------
    R0, T00, g00 = gas_props(c.p0_init, c)
    gam = np.sqrt(g00) * (2.0 / (g00 + 1.0)) ** ((g00 + 1.0) / (2.0 * (g00 - 1.0)))
    tau = c.V_g0 * gam / (res["At0"] * math.sqrt(R0 * T00))
    t_zoom = float(min(t[-1], max(25.0 * tau, 3.0e-3)))
    mz = t <= t_zoom + 1e-15
    _shade(ax_zm, spans, 0.0, t_zoom, scale=1e3)
    ax_zm.plot(t[mz] * 1e3, h["p0"][mz] / 1e3, color=C_PRESS, lw=2.0,
               label="$p_0$ [kPa]")
    ax_zm.axhline(c.p_a / 1e3, color=C_GREY, lw=0.9, ls=":", label="ambient $p_a$")
    if s["choke_limit_pressure_Pa"]:
        ax_zm.axhline(s["choke_limit_pressure_Pa"] / 1e3, color="#7a7a7a",
                      lw=0.9, ls="-.", label="choking limit")
    ax_z2 = ax_zm.twinx()
    ax_z2.grid(False)
    ax_z2.plot(t[mz] * 1e3, h["F"][mz], color=C_THRUST, lw=1.4, ls="--",
               label="$F$ [N]")
    ax_z2.set_ylabel("thrust  $F$  [N]", color=C_THRUST)
    ax_z2.tick_params(axis="y", colors=C_THRUST)
    ax_z2.spines["right"].set_visible(True)
    ax_zm.set_xlim(0.0, t_zoom * 1e3)
    ax_zm.set_xlabel("time  $t$  [ms]   (start-up zoom)")
    ax_zm.set_ylabel("$p_0$  [kPa]", color=C_PRESS)
    ax_zm.tick_params(axis="y", colors=C_PRESS)
    ax_zm.set_title(f"(e)  Ignition transient  ($\\tau_{{fill}}\\approx${tau*1e3:.2f} ms)",
                    loc="left")
    hl = ax_zm.get_legend_handles_labels()
    h2 = ax_z2.get_legend_handles_labels()
    ax_zm.legend(hl[0] + h2[0], hl[1] + h2[1], loc="best")

    # ---------------- (f) nozzle diagnostics ----------------
    _shade(ax_nz, spans)
    ax_nz.plot(t, h["Me"], color="#2d6a9f", lw=1.8, label="exit Mach $M_e$")
    ax_nz.axhline(1.0, color="#7a7a7a", lw=0.8, ls="-.")
    ax_nz.set_ylim(0.0, max(1.15, float(np.nanmax(h["Me"])) * 1.15))
    ax_nz.set_xlim(0.0, t[-1])
    ax_nz.set_xlabel("time  $t$  [s]")
    ax_nz.set_ylabel("$M_e$  [-]", color="#2d6a9f")
    ax_nz.tick_params(axis="y", colors="#2d6a9f")
    ax_n2 = ax_nz.twinx()
    ax_n2.grid(False)
    ax_n2.plot(t, h["pe"] / c.p_a, color="#9b59b6", lw=1.3, ls="--",
               label="$p_e/p_a$")
    ax_n2.plot(t, h["CF"], color="#c0392b", lw=1.3, ls=":", label="$C_F$")
    ax_n2.set_ylabel("$p_e/p_a$  ,  $C_F$  [-]", color="#7b4fa8")
    ax_n2.tick_params(axis="y", colors="#7b4fa8")
    ax_n2.spines["right"].set_visible(True)
    hl = ax_nz.get_legend_handles_labels()
    h2 = ax_n2.get_legend_handles_labels()
    ax_nz.legend(hl[0] + h2[0], hl[1] + h2[1], loc="center right")
    ax_nz.set_title("(f)  Nozzle operating point", loc="left")

    # ---------------- (g) burn rate / burn depth ----------------
    ax_r.plot(t[burn], h["x"][burn] * 1e3, color="#d2691e", lw=1.9,
              label="burn depth $x$")
    ax_r.axhline(web_thickness(c) * 1e3, color="#7a7a7a", lw=0.9, ls="-.",
                 label="web")
    ax_r.set_xlim(0.0, t[-1])
    ax_r.set_xlabel("time  $t$  [s]")
    ax_r.set_ylabel("$x$  [mm]", color="#d2691e")
    ax_r.tick_params(axis="y", colors="#d2691e")
    ax_r.axvline(t_bo, color="#404040", lw=1.0, ls="--")
    ax_r2 = ax_r.twinx()
    ax_r2.grid(False)
    ax_r2.plot(t[burn], h["r"][burn] * 1e3, color="#6b4f2a", lw=1.4, ls="--",
               label="burn rate $r$")
    ax_r2.set_ylabel("$r$  [mm/s]", color="#6b4f2a")
    ax_r2.tick_params(axis="y", colors="#6b4f2a")
    ax_r2.spines["right"].set_visible(True)
    hl = ax_r.get_legend_handles_labels()
    h2 = ax_r2.get_legend_handles_labels()
    ax_r.legend(hl[0] + h2[0], hl[1] + h2[1], loc="center right")
    ax_r.set_title("(g)  Regression of the grain", loc="left")

    # ---------------- title, regime legend, footer ----------------
    fig.suptitle("Solid Rocket Motor - Internal Ballistics\n"
                 f"bore grain ($R_i$={c.R_i0*1e3:.2f} mm, $R_p$={c.R_p*1e3:.2f} mm, "
                 f"$L_p$={c.L_p0*1e3:.1f} mm, {c.n_end:g} burning end face(s))   |   "
                 f"nozzle $R_t$={c.R_t0*1e3:.2f} mm, $A_e/A_t$={c.eps_nozzle:g}   |   "
                 f"$r=${c.a_burn*1e3:.3g} mm/s $\\times(p_0/{c.p_ref/1e6:g}\\,$MPa$)^{{{c.n_burn:g}}}$",
                 fontsize=12.5, fontweight="bold", y=0.985)

    used = [k for k in ("choked", "subsonic", "separated", "no-flow")
            if any(sp[2] == k for sp in spans)]
    handles = [Patch(facecolor=REGIME_FACE[k], edgecolor="#b8c2cc",
                     label=REGIME_NAME[k]) for k in used]
    handles.append(Line2D([0], [0], color="#404040", lw=1.0, ls="--", label="burnout"))
    fig.legend(handles=handles, loc="lower center", ncol=len(handles),
               bbox_to_anchor=(0.5, 0.012), fontsize=8.6, frameon=False)

    note = (f"GRIBS v{PROGRAM_VERSION}  |  thermochemistry: {s['thermochemistry_backend']}   |   "
            f"solver: {c.method}, rtol={c.rtol:g}   |   "
            f"unchoked policy: {c.unchoked_policy}   |   "
            f"property policy: {c.property_policy} "
            f"({s['property_range_min_Pa']/1e5:g}-{s['property_range_max_Pa']/1e5:g} bar)   |   "
            f"Cd={c.Cd:g}, eta_F={c.eta_thrust:g}, eta_T0={c.eta_T0:g}   |   "
            f"erosive burning: {'on' if c.ero_alpha > 0 else 'off'}   |   "
            f"throat erosion: {'on' if c.ero_throat_c > 0 else 'off'}")
    fig.text(0.5, 0.037, note, ha="center", va="center", fontsize=7.8,
             color="#4a5560")

    fig.savefig(path, dpi=190)
    plt.close(fig)


# ==============================================================================
# 11. FILE OUTPUT
#     CSV columns, their units and their meanings are unchanged from v0.3.1-alpha.
# ==============================================================================
CSV_COLUMNS = ["t", "phase", "p0", "x", "Ri", "Lp", "Ab", "Vg", "Kn", "Rt", "At",
               "r", "mdot_gen", "mdot_out", "m_gen", "m_out", "m_gas_eos",
               "m_gas_bal", "regime", "Me", "pe", "ve", "F", "impulse",
               "cstar", "CF", "R", "T0", "gamma", "extrap"]
CSV_UNITS = {"t": "s", "p0": "Pa", "x": "m", "Ri": "m", "Lp": "m", "Ab": "m2",
             "Vg": "m3", "Kn": "-", "Rt": "m", "At": "m2", "r": "m/s",
             "mdot_gen": "kg/s", "mdot_out": "kg/s", "m_gen": "kg", "m_out": "kg",
             "m_gas_eos": "kg", "m_gas_bal": "kg", "Me": "-", "pe": "Pa",
             "ve": "m/s", "F": "N", "impulse": "N.s", "cstar": "m/s", "CF": "-",
             "R": "J/(kg.K)", "T0": "K", "gamma": "-"}


def write_csv(h: dict, path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([f"{k}[{CSV_UNITS[k]}]" if k in CSV_UNITS else k
                    for k in CSV_COLUMNS])
        n = len(h["t"])
        for i in range(n):
            row = []
            for k in CSV_COLUMNS:
                v = h[k][i]
                row.append(v if isinstance(v, (str, np.str_, bool, np.bool_))
                           else f"{float(v):.9e}")
            w.writerow(row)


def summary_text(s: dict, c: Config, thermo_md: Optional[dict] = None,
                 cea_diag: Optional[dict] = None,
                 migration_notes: Optional[Sequence[str]] = None) -> str:
    L = []
    A = L.append
    A("=" * 78)
    A("SOLID ROCKET MOTOR - INTERNAL BALLISTICS SUMMARY")
    A(f"GRIBS v{PROGRAM_VERSION}   (schema {SCHEMA_VERSION})")
    A("=" * 78)
    A("")
    if migration_notes:
        A("[ configuration migration ]")
        for note in migration_notes:
            A(f"  WARNING: {note}")
        A("")
    A("[ configuration ]")
    A(f"  propellant      rho_p = {c.rho_p:.2f} kg/m3,  "
      f"r = {c.a_burn*1e3:.4g} mm/s x (p0/{c.p_ref/1e6:g} MPa)^{c.n_burn:g}")
    A(f"                  sigma_p = {c.sigma_p:g} 1/K, T_grain = {c.T_grain:g} K, "
      f"erosive alpha = {c.ero_alpha:g}")
    A(f"  grain           Ri0 = {c.R_i0*1e3:.3f} mm, Rp = {c.R_p*1e3:.3f} mm, "
      f"Lp0 = {c.L_p0*1e3:.3f} mm, end faces = {c.n_end:g}")
    A(f"                  web = {s['web_mm']:.3f} mm, propellant mass = "
      f"{s['propellant_mass_g']:.3f} g, Vg0 = {c.V_g0*1e6:.3f} cm3")
    A(f"  nozzle          Rt0 = {c.R_t0*1e3:.3f} mm, At0 = {math.pi*c.R_t0**2:.4e} m2, "
      f"Ae/At = {c.eps_nozzle:g}, Cd = {c.Cd:g}, eta_F = {c.eta_thrust:g}")
    A(f"  environment     pa = {c.p_a/1e3:.3f} kPa, p0(0) = {c.p0_init/1e3:.3f} kPa")
    A("")
    A("[ thermochemistry ]")
    A(f"  backend                          : {s['thermochemistry_backend']}")
    if thermo_md:
        A(f"  property range                   : "
          f"{s['property_range_min_Pa']/1e5:.6g} - {s['property_range_max_Pa']/1e5:.6g} bar "
          f"({s['property_range_min_Pa']:.6g} - {s['property_range_max_Pa']:.6g} Pa)")
        if thermo_md.get("cea_python_version"):
            A(f"  official CEA package version     : {thermo_md['cea_python_version']} "
              f"(library {thermo_md.get('cea_library_version')})")
            A(f"  official CEA module path         : {thermo_md.get('cea_python_module_path')}")
            A(f"  equilibrium formulation          : "
              f"{thermo_md.get('equilibrium_constraint')}")
        if thermo_md.get("compatibility_backend"):
            A(f"  NOTE                             : compatibility backend "
              f"({thermo_md.get('warning')})")
        A(f"  ideal mixture density            : "
          f"{thermo_md.get('ideal_mixture_density_kg_m3'):.6f} kg/m3")
        A(f"  packing fraction                 : {thermo_md.get('packing_fraction'):.6f}")
        A(f"  bulk propellant density          : "
          f"{thermo_md.get('bulk_propellant_density_kg_m3'):.6f} kg/m3")
        A(f"  cache key                        : {thermo_md.get('cache_key')}")
        A(f"  cache file                       : {thermo_md.get('cache_file')}")
        A(f"  cache reused                     : {thermo_md.get('cache_used')}")
    A(f"  property policy                  : {c.property_policy}")
    A(f"  values outside property range    : {s['property_extrapolation']} "
      f"({s['property_extrapolation_points']} sample points, "
      f"{100.0*s['property_extrapolation_fraction']:.3f} %)")
    A(f"  pressure range visited           : {s['pressure_min_Pa']:.6g} - "
      f"{s['pressure_max_Pa']:.6g} Pa")
    if cea_diag:
        A("  CEA theoretical comparison only (never used by the GRIBS model):")
        A(f"    cea_theoretical_cstar_m_s      : {cea_diag['cea_theoretical_cstar_m_s']:.3f}")
        A(f"    cea_theoretical_cf             : {cea_diag['cea_theoretical_cf']:.5f}")
        A(f"    cea_theoretical_isp_s          : {cea_diag['cea_theoretical_isp_s']:.4f}")
        A(f"    cea_theoretical_isp_vacuum_s   : {cea_diag['cea_theoretical_isp_vacuum_s']:.4f}")
        A(f"    evaluation                     : "
          f"{cea_diag['cea_theoretical_exit_definition']} at "
          f"p_c = {cea_diag['cea_theoretical_chamber_pressure_Pa']/1e6:.6f} MPa")
        A(f"    GRIBS effective c*             : {s['cstar_eff_m_s']:.3f} m/s")
        A(f"    GRIBS effective CF             : {s['CF_eff']:.5f}")
        A(f"    GRIBS effective Isp            : {s['Isp_s']:.4f} s")
        A("    note: both sets are reported side by side for context.  The CEA")
        A("          station values follow the conventions of the official package")
        A("          (throat c* is the theoretical characteristic velocity and the")
        A("          throat Cf equals v_e/c*), while cstar_eff_m_s, CF_eff and Isp_s")
        A("          are GRIBS results that include the pressure-thrust term of the")
        A("          unsteady nozzle model; they are not expected to be equal and are")
        A("          never substituted for one another.")
    A("")
    A("[ initial balance ]")
    A(f"  initial nozzle regime            : {s['initial_regime']}")
    if s["choke_limit_pressure_Pa"]:
        A(f"  choking-limit chamber pressure   : {s['choke_limit_pressure_Pa']/1e3:.3f} kPa")
    A(f"  mdot_gen(0) = {s['initial_mdot_gen_kg_s']:.6e} kg/s")
    A(f"  mdot_out(0) = {s['initial_mdot_out_kg_s']:.6e} kg/s"
      f"   ->  chamber { 'fills (pressure rises)' if s['initial_mdot_gen_kg_s'] > s['initial_mdot_out_kg_s'] else 'empties (pressure falls)'}")
    A("")
    A("[ events ]")
    for name, times in s["event_times"].items():
        if times:
            A(f"  {name:<10s}: " + ", ".join(f"{v:.6f} s" for v in times[:6]))
    A(f"  burning phase terminated by      : {s['burn_stop_reason']}")
    A(f"  blowdown phase terminated by     : {s['blowdown_stop_reason']}")
    A("")
    A("[ performance ]")
    A(f"  web burn time                    : {s['burn_time_s']:.6f} s")
    A(f"  total computed duration          : {s['total_time_s']:.6f} s")
    A(f"  max chamber pressure (MEOP)      : {s['p_max_Pa']/1e6:.6f} MPa "
      f"at t = {s['t_p_max_s']:.6f} s")
    A(f"  min chamber pressure (burning)   : {s['p_min_burn_Pa']/1e6:.6f} MPa "
      f"at t = {s['t_p_min_s']:.6f} s")
    A(f"  mean chamber pressure (burning)  : {s['p_mean_burn_Pa']/1e6:.6f} MPa")
    A(f"  max thrust                       : {s['F_max_N']:.4f} N "
      f"at t = {s['t_F_max_s']:.6f} s")
    A(f"  mean thrust (burning)            : {s['F_mean_burn_N']:.4f} N")
    A(f"  impulse, burning phase           : {s['impulse_burn_Ns']:.4f} N.s")
    A(f"  impulse, blowdown                : {s['impulse_blowdown_Ns']:.4f} N.s")
    A(f"  total impulse                    : {s['impulse_total_Ns']:.4f} N.s")
    A(f"  effective specific impulse       : {s['Isp_s']:.2f} s")
    A(f"  effective c*                     : {s['cstar_eff_m_s']:.2f} m/s")
    A(f"  effective thrust coefficient CF  : {s['CF_eff']:.5f}")
    A(f"  Kn = Ab/At  (start/max/end)      : {s['Kn_initial']:.1f} / "
      f"{s['Kn_max']:.1f} / {s['Kn_final']:.1f}")
    A(f"  final throat radius              : {s['throat_radius_final_mm']:.5f} mm")
    A(f"  nozzle regimes visited           : {', '.join(s['regimes_visited'])}")
    A("")
    A("[ verification ]")
    A(f"  propellant mass balance error    : {s['propellant_mass_balance_error']:.3e}")
    A(f"  gas-mass (EOS vs balance) error  : {s['gas_mass_consistency_error']:.3e}")
    A(f"  property-range extrapolation used: {s['property_extrapolation']}"
      f"   (policy = {c.property_policy})")
    if s["property_extrapolation"] and c.property_policy == "extrapolate":
        A("  WARNING: results rely on extrapolated thermochemical properties.")
    A("=" * 78)
    return "\n".join(L) + "\n"


# ==============================================================================
# 12. JSON CONFIGURATION (v0.4 schema, backend-specific validation)
# ==============================================================================
def _reject_unknown(node: dict, allowed: Sequence[str], path: str) -> None:
    unknown = sorted(k for k in node if k not in allowed)
    if unknown:
        raise ConfigurationError(
            f"{path}: unknown key(s) {unknown}; allowed keys: {sorted(allowed)}")


def _mapping(node: Any, path: str) -> dict:
    if not isinstance(node, dict):
        raise ConfigurationError(
            f"{path}: expected a JSON object, got {type(node).__name__}")
    return node


def _section(node: dict, key: str, path: str, required: bool = True) -> Optional[dict]:
    if key not in node:
        if required:
            raise ConfigurationError(f"{path}.{key}: missing required configuration section")
        return None
    return _mapping(node[key], f"{path}.{key}")


def _get(node: dict, key: str, path: str, required: bool = True) -> Any:
    if key not in node:
        if required:
            raise ConfigurationError(f"{path}.{key}: missing required configuration key")
        return None
    return node[key]


def _number(node: dict, key: str, path: str, *, required: bool = True,
            default: Optional[float] = None, minimum: Optional[float] = None,
            exclusive_minimum: Optional[float] = None,
            maximum: Optional[float] = None) -> Optional[float]:
    if key not in node:
        if required:
            raise ConfigurationError(f"{path}.{key}: missing required configuration key")
        return default
    value = node[key]
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigurationError(
            f"{path}.{key}: expected a number, got {value!r} ({type(value).__name__})")
    v = float(value)
    if not math.isfinite(v):
        raise ConfigurationError(f"{path}.{key}: expected a finite number, got {value!r}")
    if exclusive_minimum is not None and not v > exclusive_minimum:
        raise ConfigurationError(
            f"{path}.{key}: must be > {exclusive_minimum:g}, got {v:g}")
    if minimum is not None and v < minimum:
        raise ConfigurationError(
            f"{path}.{key}: must be >= {minimum:g}, got {v:g}")
    if maximum is not None and v > maximum:
        raise ConfigurationError(
            f"{path}.{key}: must be <= {maximum:g}, got {v:g}")
    return v


def _integer(node: dict, key: str, path: str, *, required: bool = True,
             default: Optional[int] = None, minimum: Optional[int] = None) -> Optional[int]:
    if key not in node:
        if required:
            raise ConfigurationError(f"{path}.{key}: missing required configuration key")
        return default
    value = node[key]
    if isinstance(value, bool) or not isinstance(value, int):
        if isinstance(value, float) and float(value).is_integer():
            value = int(value)
        else:
            raise ConfigurationError(
                f"{path}.{key}: expected an integer, got {value!r}")
    if minimum is not None and value < minimum:
        raise ConfigurationError(
            f"{path}.{key}: must be >= {minimum}, got {value}")
    return int(value)


def _boolean(node: dict, key: str, path: str, *, required: bool = True,
             default: Optional[bool] = None) -> Optional[bool]:
    if key not in node:
        if required:
            raise ConfigurationError(f"{path}.{key}: missing required configuration key")
        return default
    value = node[key]
    if not isinstance(value, bool):
        raise ConfigurationError(
            f"{path}.{key}: expected true/false, got {value!r}")
    return bool(value)


def _string(node: dict, key: str, path: str, *, required: bool = True,
            default: Optional[str] = None, allow_empty: bool = True) -> Optional[str]:
    if key not in node:
        if required:
            raise ConfigurationError(f"{path}.{key}: missing required configuration key")
        return default
    value = node[key]
    if not isinstance(value, str):
        raise ConfigurationError(f"{path}.{key}: expected a string, got {value!r}")
    if not allow_empty and not value.strip():
        raise ConfigurationError(f"{path}.{key}: must not be empty")
    return value


def _choice(node: dict, key: str, path: str, allowed: Sequence[str], *,
            required: bool = True, default: Optional[str] = None) -> Optional[str]:
    value = _string(node, key, path, required=required, default=default)
    if value is None:
        return None
    if value not in allowed:
        raise ConfigurationError(
            f"{path}.{key}: {value!r} is not supported; choose one of {list(allowed)}")
    return value


def _backend_choice(node: dict, path: str = "thermochemistry") -> str:
    """Read the mandatory, explicit ``backend`` selection.

    There is no default: the user must state which thermochemistry backend the
    configuration is written for, and only the v0.4 backends are accepted.
    A missing, empty, removed or unknown value produces an actionable error that
    lists the admissible choices - a backend is never assumed or substituted.
    """
    key = "backend"
    if key not in node:
        raise ConfigurationError(
            f"{path}.{key}: missing required configuration key.\n"
            f"  Set it explicitly, in the configuration file, to one of:\n"
            + "".join(f"    {name}\n" for name in KNOWN_BACKENDS)
            + "  cea_python            - official NASA CEA Python package "
              "(recommended, no external executable)\n"
              "  cea_legacy_executable - external fcea2 executable "
              "(transitional compatibility with v0.3)")
    value = node[key]
    if not isinstance(value, str):
        raise ConfigurationError(
            f"{path}.{key}: must be a string, got {type(value).__name__} "
            f"({value!r}). Choose one of {list(KNOWN_BACKENDS)}.")
    stripped = value.strip()
    if not stripped:
        raise ConfigurationError(
            f"{path}.{key}: not set (empty string).\n"
            "  The thermochemistry backend must be chosen explicitly - there is no "
            "default.\n"
            "  Set it to 'cea_python' (official NASA CEA Python package, "
            "recommended)\n"
            "  or to 'cea_legacy_executable' (external fcea2 executable).")
    if stripped in REMOVED_BACKENDS:
        raise ConfigurationError(f"{path}.{key}: {stripped!r} - "
                                 f"{REMOVED_BACKENDS[stripped]}")
    if stripped in LEGACY_BACKEND_ALIASES:
        raise ConfigurationError(
            f"{path}.{key}: {stripped!r} is the v0.3 name of "
            f"{LEGACY_BACKEND_ALIASES[stripped]!r}. v0.4 requires the new name to be "
            "written explicitly (the v0.3 name is only understood by the automatic "
            f"{LEGACY_SCHEMA_VERSION} -> {SCHEMA_VERSION} migration).")
    if stripped not in KNOWN_BACKENDS:
        raise ConfigurationError(
            f"{path}.{key}: {value!r} is not a supported backend; choose one of "
            f"{list(KNOWN_BACKENDS)}.")
    return stripped


def _string_array(node: dict, key: str, path: str, *, required: bool = True,
                  default: Sequence[str] = ()) -> Tuple[str, ...]:
    if key not in node:
        if required:
            raise ConfigurationError(f"{path}.{key}: missing required configuration key")
        return tuple(default)
    value = node[key]
    if not isinstance(value, list) or any(not isinstance(v, str) for v in value):
        raise ConfigurationError(
            f"{path}.{key}: expected a JSON array of species-name strings")
    return tuple(value)


# ---- v0.3 -> v0.4 migration --------------------------------------------------
DEFAULT_CEA_PYTHON_BLOCK: Dict[str, Any] = {
    "ions": False,
    "transport": False,
    "trace": 1.0e-10,
    "products_from_reactants": True,
    "product_species": [],
    "omit_species": [],
    "insert_species": [],
    "molecular_weight": "gas_phase_M",
    "smooth_truncation": False,
    "truncation_width": -1.0,
    "rocket_diagnostics": True,
    "cache_enabled": True,
    "rebuild_cache": False,
    "cache_directory": "",
}

_MIGRATION_BACKEND_NOTE = (
    "thermochemistry.backend: 'cea2' -> 'cea_legacy_executable'. In v0.3 the "
    "backend name 'cea2' meant the external fcea2 executable; the physical "
    "calculation is unchanged. Select 'cea_python' to use the official NASA CEA "
    "Python package.")


def migrate_v030_document(doc: dict) -> Tuple[dict, List[str]]:
    """Explicitly migrate a v0.3 configuration document to the v0.4 schema.

    Nothing is reinterpreted silently: every structural change is recorded in the
    returned note list, which is printed as a warning and stored in summary.json.
    """
    notes: List[str] = []
    new = copy.deepcopy(_mapping(doc, "root"))
    th = _mapping(new.get("thermochemistry", {}), "thermochemistry")

    old_backend = th.get("backend")
    if old_backend in REMOVED_BACKENDS:
        # Removed backends are never silently replaced by another one.
        raise ConfigurationError(
            f"thermochemistry.backend: {old_backend!r} in a "
            f"{LEGACY_SCHEMA_VERSION} configuration - {REMOVED_BACKENDS[old_backend]}\n"
            "No other backend is substituted automatically: edit the configuration "
            "and choose explicitly.")
    if old_backend in LEGACY_BACKEND_ALIASES:
        th["backend"] = LEGACY_BACKEND_ALIASES[old_backend]
        notes.append(_MIGRATION_BACKEND_NOTE)
    elif old_backend in KNOWN_BACKENDS:
        notes.append(f"thermochemistry.backend: {old_backend!r} is already a v0.4 name.")
    else:
        raise ConfigurationError(
            f"thermochemistry.backend: cannot migrate unknown v0.3 backend "
            f"{old_backend!r}. Set it explicitly to one of: {', '.join(KNOWN_BACKENDS)}.")

    if "legacy_fit" in th:
        th.pop("legacy_fit")
        notes.append(
            "thermochemistry.legacy_fit: dropped. The manual R/T0/gamma correlation "
            "backend was removed in v0.4.0-alpha together with its "
            "propellant_density_kg_m3 setting; the bulk propellant density is now "
            "always rho_p = packing_fraction / sum_i(w_i/rho_i).")

    cea2 = th.pop("cea2", None)
    if cea2 is not None:
        cea2 = _mapping(cea2, "thermochemistry.cea2")
        th["pressure_points"] = int(cea2.get("pressure_points", 81))
        th["cea_legacy_executable"] = {
            "executable": str(cea2.get("executable", "fcea2")),
            "data_directory": str(cea2.get("data_directory", "")),
            "timeout_s": float(cea2.get("timeout_s", 120.0)),
            "batch_size": 8,
            "trace": float(cea2.get("trace", 1.0e-10)),
            "cache_enabled": True,
            "rebuild_cache": bool(cea2.get("rebuild_cache", False)),
            "cache_directory": str(cea2.get("cache_directory", "")),
        }
        notes.append(
            "thermochemistry.cea2.* -> thermochemistry.cea_legacy_executable.* "
            "(executable, data_directory, timeout_s, trace, cache settings); "
            "cea2.pressure_points -> thermochemistry.pressure_points; the v0.3 "
            "eight-pressure batch workaround is retained as batch_size = 8. "
            "New fields (cache_enabled) take their documented defaults.")
    else:
        th.setdefault("cea_legacy_executable", {
            "executable": "fcea2", "data_directory": "", "timeout_s": 120.0,
            "batch_size": 8, "trace": 1.0e-10, "cache_enabled": True,
            "rebuild_cache": False, "cache_directory": "",
        })
        notes.append("thermochemistry.cea_legacy_executable: not present in the v0.3 "
                     "document, defaults added.")
    th.setdefault("pressure_points", 81)

    if "cea_python" not in th:
        th["cea_python"] = copy.deepcopy(DEFAULT_CEA_PYTHON_BLOCK)
        notes.append("thermochemistry.cea_python: not present in the v0.3 document, "
                     "defaults added (the official package backend is available but "
                     "not selected by this migration).")
    new["schema_version"] = SCHEMA_VERSION
    notes.append(f"schema_version: {LEGACY_SCHEMA_VERSION!r} -> {SCHEMA_VERSION!r}.")
    return new, notes


@dataclass
class Configuration:
    config: Config
    output_dir: Path
    output_names: Dict[str, str]
    document: dict
    migration_notes: List[str]


_ROOT_KEYS = ("schema_version", "notes", "propellant", "grain", "nozzle",
              "environment", "igniter", "thermochemistry", "solver", "output")
_THERMO_ROOT_KEYS = ("backend", "temperature_efficiency", "outside_range_policy",
                     "pressure_range", "pressure_points", "cea_python",
                     "cea_legacy_executable")
_PROP_ROOT_KEYS = ("description", "packing_fraction", "reactants", "burn_law")
_BURN_LAW_KEYS = ("coefficient_m_s", "pressure_exponent", "reference_pressure_Pa",
                  "temperature_sensitivity_1_K", "grain_temperature_K",
                  "reference_temperature_K", "erosive_burning")
_ERO_BURN_KEYS = ("enabled", "alpha", "beta")
_REACTANT_KEYS = ("name", "wt_percent", "temperature_K", "density_kg_m3")
_GRAIN_KEYS = ("initial_bore_radius_m", "outer_radius_m", "initial_length_m",
               "initial_free_volume_m3", "burning_end_faces")
_NOZZLE_KEYS = ("initial_throat_radius_m", "expansion_ratio", "discharge_coefficient",
                "thrust_efficiency", "flow_separation", "throat_erosion")
_SEPARATION_KEYS = ("enabled", "pressure_ratio")
_THROAT_EROSION_KEYS = ("enabled", "rate_m_s_at_reference_pressure", "pressure_exponent")
_ENV_KEYS = ("ambient_pressure_Pa", "initial_chamber_pressure_Pa")
_IGNITER_KEYS = ("mass_flow_kg_s", "duration_s")
_PRESSURE_RANGE_KEYS = ("minimum_Pa", "maximum_Pa")
_CEA_PY_KEYS = ("ions", "transport", "trace", "products_from_reactants",
                "product_species", "omit_species", "insert_species",
                "molecular_weight", "smooth_truncation", "truncation_width",
                "rocket_diagnostics", "cache_enabled", "rebuild_cache",
                "cache_directory")
_CEA_LEGACY_KEYS = ("executable", "data_directory", "timeout_s", "batch_size",
                    "trace", "cache_enabled", "rebuild_cache", "cache_directory")
_SOLVER_KEYS = ("method", "relative_tolerance", "pressure_absolute_tolerance_Pa",
                "burn_depth_absolute_tolerance_m", "burning_time_limit_s",
                "maximum_burning_step_s", "unchoked_policy", "blowdown")
_BLOWDOWN_KEYS = ("enabled", "time_limit_s", "maximum_step_s")
_OUTPUT_KEYS = ("directory", "figure_filename", "history_filename",
                "summary_text_filename", "summary_json_filename")


def _resolve_from_script(value: str) -> Path:
    p = Path(value).expanduser()
    return p if p.is_absolute() else Path(__file__).resolve().parent / p


def _load_document(path: Path) -> Tuple[dict, List[str]]:
    path = path.expanduser().resolve()
    if not path.is_file():
        raise ConfigurationError(
            f"Configuration file not found: {path}\n"
            "Place gribs_config.json beside the Python file or use --config.")
    try:
        doc = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ConfigurationError(
            f"{path}: configuration is not valid JSON.\nOriginal error: {exc}") from exc
    if not isinstance(doc, dict):
        raise ConfigurationError(f"{path}: the configuration root must be a JSON object.")
    schema = doc.get("schema_version")
    notes: List[str] = []
    if schema == LEGACY_SCHEMA_VERSION:
        doc, notes = migrate_v030_document(doc)
        notes.insert(0, f"{path.name}: {LEGACY_SCHEMA_VERSION} configuration migrated "
                        f"automatically to {SCHEMA_VERSION}.")
    elif schema != SCHEMA_VERSION:
        raise ConfigurationError(
            f"schema_version: expected {SCHEMA_VERSION!r} (or the migratable legacy "
            f"value {LEGACY_SCHEMA_VERSION!r}), got {schema!r}.\n"
            "Update the configuration file; see docs/MIGRATION_v0.3_to_v0.4.md.")
    return doc, notes


def configuration_from_document(doc: dict) -> Configuration:
    """Validate a v0.4 document and build the runtime configuration."""
    root = _mapping(doc, "root")
    _reject_unknown(root, _ROOT_KEYS, "root")
    # optional free-text notes carried by the configuration (documentation only)
    _string_array(root, "notes", "root", required=False)

    prop = _section(root, "propellant", "root")
    _reject_unknown(prop, _PROP_ROOT_KEYS, "propellant")
    burn = _mapping(_get(prop, "burn_law", "propellant"), "propellant.burn_law")
    _reject_unknown(burn, _BURN_LAW_KEYS, "propellant.burn_law")
    ero = _mapping(_get(burn, "erosive_burning", "propellant.burn_law"),
                   "propellant.burn_law.erosive_burning")
    _reject_unknown(ero, _ERO_BURN_KEYS, "propellant.burn_law.erosive_burning")

    grain = _section(root, "grain", "root")
    _reject_unknown(grain, _GRAIN_KEYS, "grain")
    nozzle = _section(root, "nozzle", "root")
    _reject_unknown(nozzle, _NOZZLE_KEYS, "nozzle")
    sep = _mapping(_get(nozzle, "flow_separation", "nozzle"), "nozzle.flow_separation")
    _reject_unknown(sep, _SEPARATION_KEYS, "nozzle.flow_separation")
    tero = _mapping(_get(nozzle, "throat_erosion", "nozzle"), "nozzle.throat_erosion")
    _reject_unknown(tero, _THROAT_EROSION_KEYS, "nozzle.throat_erosion")
    env = _section(root, "environment", "root")
    _reject_unknown(env, _ENV_KEYS, "environment")
    ign = _section(root, "igniter", "root")
    _reject_unknown(ign, _IGNITER_KEYS, "igniter")
    th = _section(root, "thermochemistry", "root")
    _reject_unknown(th, _THERMO_ROOT_KEYS, "thermochemistry")
    solver = _section(root, "solver", "root")
    _reject_unknown(solver, _SOLVER_KEYS, "solver")
    blow = _mapping(_get(solver, "blowdown", "solver"), "solver.blowdown")
    _reject_unknown(blow, _BLOWDOWN_KEYS, "solver.blowdown")
    out = _section(root, "output", "root")
    _reject_unknown(out, _OUTPUT_KEYS, "output")

    # ---------------- propellant composition (always required) ----------------
    packing = float(_number(prop, "packing_fraction", "propellant",
                            exclusive_minimum=0.0, maximum=1.0))
    reactants_raw = _get(prop, "reactants", "propellant")
    if not isinstance(reactants_raw, list) or not reactants_raw:
        raise ConfigurationError(
            "propellant.reactants: expected a non-empty JSON array of constituents")
    reactants: List[ReactantSpec] = []
    total = 0.0
    for i, entry in enumerate(reactants_raw):
        path = f"propellant.reactants[{i}]"
        item = _mapping(entry, path)
        _reject_unknown(item, _REACTANT_KEYS, path)
        name = _string(item, "name", path, allow_empty=False)
        wt = float(_number(item, "wt_percent", path, exclusive_minimum=0.0))
        temp = float(_number(item, "temperature_K", path, exclusive_minimum=0.0))
        rho = float(_number(item, "density_kg_m3", path, exclusive_minimum=0.0))
        total += wt
        reactants.append(ReactantSpec(name=name, wt_percent=wt,
                                      temperature_K=temp, density_kg_m3=rho))
    if abs(total - 100.0) > 1e-6:
        raise ConfigurationError(
            f"propellant.reactants[*].wt_percent: the mass percentages must sum to "
            f"100, got {total:.12g} (difference {total-100.0:+.3e} percentage points).")

    # ---------------- thermochemistry block ----------------------------------
    backend = _backend_choice(th, "thermochemistry")
    th_keys = _THERMO_ROOT_KEYS
    prange = _mapping(_get(th, "pressure_range", "thermochemistry"),
                      "thermochemistry.pressure_range")
    _reject_unknown(prange, _PRESSURE_RANGE_KEYS, "thermochemistry.pressure_range")
    pmin = float(_number(prange, "minimum_Pa", "thermochemistry.pressure_range",
                         exclusive_minimum=0.0))
    pmax = float(_number(prange, "maximum_Pa", "thermochemistry.pressure_range",
                         exclusive_minimum=0.0))
    if pmax <= pmin:
        raise ConfigurationError(
            "thermochemistry.pressure_range: maximum_Pa must exceed minimum_Pa "
            f"(got minimum_Pa = {pmin:g}, maximum_Pa = {pmax:g})")
    points = _integer(th, "pressure_points", "thermochemistry", required=False, default=81,
                      minimum=12)

    cea_py = _section(th, "cea_python", "thermochemistry",
                      required=(backend == BACKEND_CEA_PYTHON))
    if cea_py is not None:
        _reject_unknown(cea_py, _CEA_PY_KEYS, "thermochemistry.cea_python")
    cea_legacy = _section(th, "cea_legacy_executable", "thermochemistry",
                          required=(backend == BACKEND_CEA_LEGACY_EXECUTABLE))
    if cea_legacy is not None:
        _reject_unknown(cea_legacy, _CEA_LEGACY_KEYS, "thermochemistry.cea_legacy_executable")

    # Backend-specific defaults are only validated when the section is present, so
    # a cea_python configuration is never forced to carry fcea2 settings and vice
    # versa.
    def _py(key, default):
        if cea_py is None:
            return default
        return cea_py.get(key, default)

    def _legacy(key, default):
        if cea_legacy is None:
            return default
        return cea_legacy.get(key, default)

    cea_py_block = _mapping(cea_py or copy.deepcopy(DEFAULT_CEA_PYTHON_BLOCK),
                            "thermochemistry.cea_python")
    py_ions = _boolean(cea_py_block, "ions", "thermochemistry.cea_python",
                       required=False, default=False)
    py_transport = _boolean(cea_py_block, "transport", "thermochemistry.cea_python",
                            required=False, default=False)
    py_trace = float(_number(cea_py_block, "trace", "thermochemistry.cea_python",
                             required=False, default=1.0e-10))
    py_products = _boolean(cea_py_block, "products_from_reactants",
                           "thermochemistry.cea_python", required=False, default=True)
    py_species = _string_array(cea_py_block, "product_species",
                               "thermochemistry.cea_python", required=False)
    py_omit = _string_array(cea_py_block, "omit_species",
                            "thermochemistry.cea_python", required=False)
    py_insert = _string_array(cea_py_block, "insert_species",
                              "thermochemistry.cea_python", required=False)
    py_mw = _choice(cea_py_block, "molecular_weight", "thermochemistry.cea_python",
                    ("gas_phase_M", "total_MW"), required=False, default="gas_phase_M")
    py_smooth = _boolean(cea_py_block, "smooth_truncation",
                         "thermochemistry.cea_python", required=False, default=False)
    py_width = float(_number(cea_py_block, "truncation_width",
                             "thermochemistry.cea_python", required=False, default=-1.0))
    py_rocket = _boolean(cea_py_block, "rocket_diagnostics",
                         "thermochemistry.cea_python", required=False, default=True)
    py_cache = _boolean(cea_py_block, "cache_enabled", "thermochemistry.cea_python",
                        required=False, default=True)
    py_rebuild = _boolean(cea_py_block, "rebuild_cache", "thermochemistry.cea_python",
                          required=False, default=False)
    py_cache_dir = _string(cea_py_block, "cache_directory", "thermochemistry.cea_python",
                           required=False, default="")

    legacy_block = _mapping(cea_legacy or {}, "thermochemistry.cea_legacy_executable")
    lg_exe = _string(legacy_block, "executable", "thermochemistry.cea_legacy_executable",
                     required=False, default="fcea2")
    lg_dir = _string(legacy_block, "data_directory",
                     "thermochemistry.cea_legacy_executable", required=False, default="")
    lg_timeout = float(_number(legacy_block, "timeout_s",
                               "thermochemistry.cea_legacy_executable",
                               required=False, default=120.0, exclusive_minimum=0.0))
    lg_batch = _integer(legacy_block, "batch_size",
                        "thermochemistry.cea_legacy_executable",
                        required=False, default=8, minimum=1)
    lg_trace = float(_number(legacy_block, "trace",
                             "thermochemistry.cea_legacy_executable",
                             required=False, default=1.0e-10))
    lg_cache = _boolean(legacy_block, "cache_enabled",
                        "thermochemistry.cea_legacy_executable",
                        required=False, default=True)
    lg_rebuild = _boolean(legacy_block, "rebuild_cache",
                          "thermochemistry.cea_legacy_executable",
                          required=False, default=False)
    lg_cache_dir = _string(legacy_block, "cache_directory",
                           "thermochemistry.cea_legacy_executable",
                           required=False, default="")

    if backend == BACKEND_CEA_LEGACY_EXECUTABLE:
        # the compatibility backend really needs these two to be usable
        _string(legacy_block, "executable", "thermochemistry.cea_legacy_executable",
                allow_empty=False)
        _string(legacy_block, "data_directory", "thermochemistry.cea_legacy_executable")

    if backend == BACKEND_CEA_PYTHON:
        cache_enabled, rebuild_cache, cache_dir = bool(py_cache), bool(py_rebuild), str(py_cache_dir)
    elif backend == BACKEND_CEA_LEGACY_EXECUTABLE:
        cache_enabled, rebuild_cache, cache_dir = bool(lg_cache), bool(lg_rebuild), str(lg_cache_dir)
    else:
        cache_enabled, rebuild_cache, cache_dir = False, False, ""

    data = dict(
        rho_p=0.0,          # recomputed from the mixture density below
        a_burn=float(_number(burn, "coefficient_m_s", "propellant.burn_law")),
        n_burn=float(_number(burn, "pressure_exponent", "propellant.burn_law")),
        p_ref=float(_number(burn, "reference_pressure_Pa", "propellant.burn_law",
                            exclusive_minimum=0.0)),
        sigma_p=float(_number(burn, "temperature_sensitivity_1_K", "propellant.burn_law")),
        T_grain=float(_number(burn, "grain_temperature_K", "propellant.burn_law",
                              exclusive_minimum=0.0)),
        T_ref=float(_number(burn, "reference_temperature_K", "propellant.burn_law",
                            exclusive_minimum=0.0)),
        ero_alpha=(float(_number(ero, "alpha", "propellant.burn_law.erosive_burning"))
                   if _boolean(ero, "enabled", "propellant.burn_law.erosive_burning") else 0.0),
        ero_beta=float(_number(ero, "beta", "propellant.burn_law.erosive_burning")),
        R_i0=float(_number(grain, "initial_bore_radius_m", "grain", exclusive_minimum=0.0)),
        R_p=float(_number(grain, "outer_radius_m", "grain", exclusive_minimum=0.0)),
        L_p0=float(_number(grain, "initial_length_m", "grain", exclusive_minimum=0.0)),
        V_g0=float(_number(grain, "initial_free_volume_m3", "grain", exclusive_minimum=0.0)),
        n_end=float(_number(grain, "burning_end_faces", "grain", minimum=0.0)),
        R_t0=float(_number(nozzle, "initial_throat_radius_m", "nozzle", exclusive_minimum=0.0)),
        eps_nozzle=float(_number(nozzle, "expansion_ratio", "nozzle", minimum=1.0)),
        Cd=float(_number(nozzle, "discharge_coefficient", "nozzle", exclusive_minimum=0.0)),
        eta_thrust=float(_number(nozzle, "thrust_efficiency", "nozzle", exclusive_minimum=0.0)),
        use_separation=bool(_boolean(sep, "enabled", "nozzle.flow_separation")),
        sep_ratio=float(_number(sep, "pressure_ratio", "nozzle.flow_separation",
                                exclusive_minimum=0.0)),
        ero_throat_c=(float(_number(tero, "rate_m_s_at_reference_pressure",
                                    "nozzle.throat_erosion", minimum=0.0))
                      if _boolean(tero, "enabled", "nozzle.throat_erosion") else 0.0),
        ero_throat_m=float(_number(tero, "pressure_exponent", "nozzle.throat_erosion")),
        p_a=float(_number(env, "ambient_pressure_Pa", "environment", exclusive_minimum=0.0)),
        p0_init=float(_number(env, "initial_chamber_pressure_Pa", "environment",
                              exclusive_minimum=0.0)),
        ign_mdot=float(_number(ign, "mass_flow_kg_s", "igniter", minimum=0.0)),
        ign_time=float(_number(ign, "duration_s", "igniter", minimum=0.0)),
        thermo_backend=backend,
        cea_pressure_points=int(points),
        cea_cache_enabled=cache_enabled,
        cea_rebuild_cache=rebuild_cache,
        cea_cache_dir=cache_dir,
        cea_py_ions=bool(py_ions),
        cea_py_transport=bool(py_transport),
        cea_py_trace=float(py_trace),
        cea_py_products_from_reactants=bool(py_products),
        cea_py_product_species=tuple(py_species),
        cea_py_omit_species=tuple(py_omit),
        cea_py_insert_species=tuple(py_insert),
        cea_py_molecular_weight=str(py_mw),
        cea_py_smooth_truncation=bool(py_smooth),
        cea_py_truncation_width=float(py_width),
        cea_py_rocket_diagnostics=bool(py_rocket),
        cea_legacy_executable_path=str(lg_exe),
        cea_legacy_data_directory=str(lg_dir),
        cea_legacy_timeout_s=float(lg_timeout),
        cea_legacy_batch_size=int(lg_batch),
        cea_legacy_trace=float(lg_trace),
        p_fit_min=float(pmin), p_fit_max=float(pmax),
        eta_T0=float(_number(th, "temperature_efficiency", "thermochemistry",
                             exclusive_minimum=0.0, maximum=1.5)),
        property_policy=str(_choice(th, "outside_range_policy", "thermochemistry",
                                    ("clamp", "extrapolate"))),
        unchoked_policy=str(_choice(solver, "unchoked_policy", "solver",
                                    ("switch", "stop"))),
        method=str(_choice(solver, "method", "solver", ("LSODA", "BDF", "Radau"))),
        rtol=float(_number(solver, "relative_tolerance", "solver", exclusive_minimum=0.0)),
        atol_p=float(_number(solver, "pressure_absolute_tolerance_Pa", "solver",
                             exclusive_minimum=0.0)),
        atol_x=float(_number(solver, "burn_depth_absolute_tolerance_m", "solver",
                             exclusive_minimum=0.0)),
        t_max=float(_number(solver, "burning_time_limit_s", "solver", exclusive_minimum=0.0)),
        max_step_burn=float(_number(solver, "maximum_burning_step_s", "solver",
                                    exclusive_minimum=0.0)),
        blowdown=bool(_boolean(blow, "enabled", "solver.blowdown")),
        blowdown_tmax=float(_number(blow, "time_limit_s", "solver.blowdown",
                                    exclusive_minimum=0.0)),
        max_step_blow=float(_number(blow, "maximum_step_s", "solver.blowdown",
                                    exclusive_minimum=0.0)),
    )
    c = Config(**data)
    # The bulk propellant density is fully determined by the configuration: the
    # additive-volume mixture density of the constituents times the packing
    # fraction.  It never comes from a gas property and there is no manual override.
    c.rho_p = ideal_mixture_density(reactants) * packing
    c.configuration_source = ""                     # type: ignore[attr-defined]
    c.original_configuration = copy.deepcopy(doc)  # type: ignore[attr-defined]
    c.reactants_spec_list = reactants               # type: ignore[attr-defined]
    c.packing_fraction = packing                    # type: ignore[attr-defined]

    output_dir = _resolve_from_script(str(_get(out, "directory", "output"))).resolve()
    names = {k: str(_get(out, k, "output")) for k in
             ("figure_filename", "history_filename", "summary_text_filename",
              "summary_json_filename")}
    return Configuration(config=c, output_dir=output_dir, output_names=names,
                         document=doc, migration_notes=[])


def load_json_configuration(path: Path) -> Configuration:
    doc, notes = _load_document(path)
    loaded = configuration_from_document(doc)
    loaded.migration_notes = list(notes)
    loaded.config.configuration_source = str(Path(path).expanduser().resolve())
    return loaded


# ==============================================================================
# 13. CLI AND MAIN PROGRAM
# ==============================================================================
def environment_info() -> Dict[str, Any]:
    info: Dict[str, Any] = {
        "program": PROGRAM_NAME,
        "program_version": PROGRAM_VERSION,
        "schema_version": SCHEMA_VERSION,
        "python_version": sys.version.split()[0],
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "system": platform.system(),
        "machine": platform.machine(),
        "numpy_version": np.__version__,
        "scipy_version": __import__("scipy").__version__,
        "matplotlib_version": matplotlib.__version__,
    }
    try:
        import cea  # noqa: F401  (only for provenance)
        info["cea_python_version"] = getattr(cea, "__version__", None)
        info["cea_python_module_path"] = getattr(cea, "__file__", None)
        info["cea_library_version"] = str(cea.lib_version())
    except Exception as exc:
        info["cea_python_version"] = None
        info["cea_python_import_error"] = f"{type(exc).__name__}: {exc}"
    return info


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=f"GRIBS v{PROGRAM_VERSION} JSON-configured internal ballistics")
    p.add_argument("--config", type=Path,
                   default=Path(__file__).resolve().parent / "gribs_config.json",
                   help="complete JSON configuration file")
    p.add_argument("--backend", choices=list(KNOWN_BACKENDS), default=None,
                   help="override thermochemistry.backend for this run")
    p.add_argument("--selftest", action="store_true", help="run numerical self-tests only")
    p.add_argument("--validate-config", action="store_true",
                   help="validate the JSON configuration and exit")
    p.add_argument("--dump-thermo", type=Path, default=None, metavar="FILE",
                   help="write the chamber-property table to FILE and exit")
    p.add_argument("--no-cache", action="store_true",
                   help="force a rebuild of the chamber-property table")
    p.add_argument("--version", action="version",
                   version=f"{PROGRAM_NAME} {PROGRAM_VERSION}")
    return p


def _print_backend_banner(cfg: Config, backend: ThermochemistryBackend) -> None:
    print(f"{PROGRAM_NAME} v{PROGRAM_VERSION}  (schema {SCHEMA_VERSION})")
    print(f"thermochemistry backend : {backend.name}")
    if backend.name == BACKEND_CEA_PYTHON:
        print(f"  official CEA package  : {getattr(backend.cea, '__version__', '?')} "
              f"at {getattr(backend.cea, '__file__', '?')}")
        print(f"  CEA library version   : {backend.cea.lib_version()}")
        print(f"  equilibrium           : HP (assigned enthalpy and pressure)")
    elif backend.name == BACKEND_CEA_LEGACY_EXECUTABLE:
        print("  COMPATIBILITY BACKEND (transitional, v0.4-alpha only).")
        print("  Prefer 'cea_python' (official package) for new work.")
    print(f"  property table        : {backend.c.p_fit_min/1e5:g} - "
          f"{backend.c.p_fit_max/1e5:g} bar, {backend.c.cea_pressure_points} points, "
          f"eta_T0 = {backend.c.eta_T0:g}, policy = {backend.c.property_policy}")
    if getattr(backend, "cache_path", None):
        print(f"  cache file            : {backend.cache_path}")
        print(f"  cache reused          : {backend.cache_used}")
    print(f"  ideal mixture density : {backend.ideal_mixture_density:.6f} kg/m3")
    print(f"  packing fraction      : {backend.packing_fraction:.6f}")
    print(f"  bulk rho_p used       : {cfg.rho_p:.6f} kg/m3")


def _cea_theoretical_diagnostics(cfg: Config, backend: ThermochemistryBackend,
                                 summ: dict) -> Optional[dict]:
    """Theoretical CEA rocket comparison values (never fed back into GRIBS)."""
    if not getattr(cfg, "cea_py_rocket_diagnostics", False):
        return None
    if not hasattr(backend, "theoretical_rocket_performance"):
        return None
    p_ref = summ.get("p_mean_burn_Pa")
    if not (isinstance(p_ref, float) and math.isfinite(p_ref) and p_ref > 0.0):
        p_ref = summ.get("p_max_Pa")
    try:
        return backend.theoretical_rocket_performance(float(p_ref), float(cfg.eps_nozzle))
    except Exception as exc:                                   # diagnostic only
        print(f"WARNING: CEA theoretical rocket diagnostics unavailable: "
              f"{type(exc).__name__}: {exc}")
        return {"cea_theoretical_error": f"{type(exc).__name__}: {exc}"}


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        loaded = load_json_configuration(args.config)
    except ConfigurationError as exc:
        print(f"Configuration error:\n{exc}", file=sys.stderr)
        return 2

    cfg = loaded.config
    if args.backend:
        cfg.thermo_backend = args.backend
        loaded.migration_notes.append(
            f"thermochemistry.backend overridden on the command line to "
            f"{args.backend!r}.")
    if args.no_cache:
        cfg.cea_rebuild_cache = True
        loaded.migration_notes.append("Chamber-property cache rebuild forced by --no-cache.")

    for note in loaded.migration_notes:
        print(f"WARNING: {note}", file=sys.stderr)

    try:
        validate(cfg)
    except (ConfigurationError, ValueError) as exc:
        print(f"Configuration error:\n{exc}", file=sys.stderr)
        return 2

    if args.validate_config:
        print(f"Configuration {loaded.config.configuration_source} is valid "
              f"(schema {SCHEMA_VERSION}, backend {cfg.thermo_backend}).")
        return 0

    cache_dir = loaded.output_dir / "thermo_cache"
    try:
        backend = initialize_thermochemistry(cfg, cache_dir)
    except (ThermochemistryError, ConfigurationError) as exc:
        print(f"Thermochemistry backend error:\n{exc}", file=sys.stderr)
        return 3

    _print_backend_banner(cfg, backend)

    if args.dump_thermo is not None:
        payload = {
            "backend": backend.name,
            "provenance": backend.metadata(),
            "pressure_Pa": list(backend.table.get("pressure_Pa", [])),
            "temperature_K": list(backend.table.get("temperature_K", [])),
            "gamma_s": list(backend.table.get("gamma_s", [])),
            "molecular_weight_kg_kmol": list(backend.table.get("molecular_weight_kg_kmol", [])),
            "gas_constant_J_kgK": list(backend.table.get("gas_constant_J_kgK", [])),
        }
        args.dump_thermo.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        print(f"Chamber-property table written to {args.dump_thermo}")
        return 0

    if args.selftest:
        return 0 if self_tests(cfg) else 1

    if not self_tests(cfg):
        print("Self-tests failed; aborting before the production run.")
        return 1

    print("Running the internal-ballistics simulation ...")
    try:
        res = run_model(cfg)
        hist = sample(res, cfg)
    except (GribsError, ValueError, RuntimeError) as exc:
        print(f"Simulation failed: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 4
    summ = summarize(res, hist, cfg)
    thermo_md = backend.metadata()
    cea_diag = _cea_theoretical_diagnostics(cfg, backend, summ)

    outdir = loaded.output_dir
    outdir.mkdir(parents=True, exist_ok=True)
    png = outdir / loaded.output_names["figure_filename"]
    csvp = outdir / loaded.output_names["history_filename"]
    txtp = outdir / loaded.output_names["summary_text_filename"]
    jsnp = outdir / loaded.output_names["summary_json_filename"]

    write_csv(hist, csvp)
    make_figure(res, hist, summ, cfg, png, cea_diag)
    text = summary_text(summ, cfg, thermo_md, cea_diag, loaded.migration_notes)
    txtp.write_text(text, encoding="utf-8")

    derived = {
        "propellant_density_kg_m3": cfg.rho_p,
        "ideal_mixture_density_kg_m3": backend.ideal_mixture_density,
        "packing_fraction": backend.packing_fraction,
        "bulk_propellant_density_kg_m3": backend.bulk_mixture_density,
        "propellant_density_model": (
            "rho_p = packing_fraction / sum_i(w_i/rho_i); independent of any CEA gas density"),
    }
    summary_json = {
        "program": {"name": PROGRAM_NAME, "version": PROGRAM_VERSION,
                    "schema_version": SCHEMA_VERSION},
        "environment": environment_info(),
        "configuration_source": cfg.configuration_source,
        "configuration_migration": ({"applied": True, "notes": list(loaded.migration_notes)}
                                    if loaded.migration_notes else {"applied": False,
                                                                    "notes": []}),
        "original_configuration": loaded.document,
        "resolved_configuration": asdict(cfg),
        "derived_inputs": derived,
        "thermochemistry": thermo_md,
        "cea_theoretical_diagnostics": cea_diag,
        "results": summ,
        "output_files": {"figure": str(png), "history_csv": str(csvp),
                         "summary_text": str(txtp), "summary_json": str(jsnp)},
    }
    jsnp.write_text(json.dumps(summary_json, indent=2, ensure_ascii=False) + "\n",
                    encoding="utf-8")

    print()
    print(text)
    print("Output files")
    for pth in (png, csvp, txtp, jsnp):
        print(f"  {pth}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
