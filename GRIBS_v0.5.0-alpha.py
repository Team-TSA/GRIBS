#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GRIBS v0.5.0-alpha
================================================================================
JSON-CONFIGURED INTERNAL BALLISTICS OF A SOLID ROCKET MOTOR
  Grain    : cylindrical bore + N burning end face(s) (outer surface inhibited)
  Nozzle   : converging (Ae = At) or converging-diverging (Ae/At > 1)
  Chamber  : unsteady mass/state equation, HOMOGENEOUS EQUILIBRIUM TWO-PHASE
             model B: the conserved chamber mass is the TOTAL product mass
             (gas + condensed), the equation of state acts on the gas phase only

--------------------------------------------------------------------------------
GOVERNING MODEL  (two-phase model B, new in v0.5.0-alpha)
--------------------------------------------------------------------------------
  Chamber composition at every tabulated pressure comes from the CEA HP
  equilibrium and is split, with the official CEA phase information, into
      Yg(p)  gas-phase mass fraction      (sum of the gas-species mass fractions)
      Yc(p)  condensed mass fraction      (sum of the condensed-species mass
                                           fractions),  Yg + Yc = 1
  The conserved chamber mass is the TOTAL product mass  mt = mg + mc  with
      mg = Yg(p)*mt ,  mc = Yc(p)*mt .
  The condensed phase volume is neglected against the chamber free volume, so
  the gas phase occupies Vg and the equation of state is

      p0 * Vg = mt * Yg(p) * Rg(p) * T0(p)  =  mt * Psi(p),
      Psi(p) = Yg(p) * Rg(p) * T0(p)          [J/kg]   (gas-phase R and T only)

  Differentiating  p0*Vg = mt*Psi(p0)  exactly in time,
      Vg*dp0/dt + p0*dVg/dt = Psi*dmt/dt + mt*dPsi/dp*dp0/dt ,
  and substituting  mt = p0*Vg/Psi ,  dVg/dt = Ab*r ,
  dmt/dt = mdot_gen_total - mdot_out_total  gives the pressure equation used
  here (valid for any nozzle regime):

      dp0/dt = [ Psi*(mdot_gen_total - mdot_out_total) - p0*Ab*r ]
               / [ Vg*(1 - p0*Psi'/Psi) ]

  During blowdown (r = 0, Ab = 0, mdot_gen_total = 0):
      dp0/dt = - Psi*mdot_out_total / [ Vg*(1 - p0*Psi'/Psi) ]

  The pressure dependence of Yg(p) (equilibrium condensation / re-vaporization)
  enters the ODE through Psi'(p); dPsi/dp is the ANALYTIC derivative of the
  PCHIP interpolator of Psi in ln(p), never a finite difference.

  Two-phase assumptions (explicit, see summary outputs):
    * gas and condensed phase are in local thermochemical equilibrium and share
      one representative temperature (the CEA HP equilibrium temperature);
    * the phases are homogeneously mixed in the chamber (no settling, no wall
      deposition, no slag, no particle-size distribution);
    * the condensed volume is negligible against the chamber free volume;
    * the nozzle flow is the "complete entrainment" first approximation: the
      gas-phase choked reference flow divided by the chamber gas mass fraction,
      both phases leaving with the SAME exit velocity (no particle slip, no
      non-equilibrium nozzle chemistry);
    * no separate energy equation: T0 = eta_T0 * T_CEA(p) as before;
    * results are NOT experimentally validated - independent validation is
      required before any engineering use.

  LEGACY SINGLE-PHASE MODE (regression comparison only):
      two_phase_model.mode = "single_phase_legacy" reproduces the v0.4.3-alpha
      chamber model exactly (Yg := 1, Psi := Theta = R*T0, total flow = gas
      reference flow).  It is retained ONLY for regression comparison; it still
      feeds the full generated product mass into a gas-only equation of state,
      which is inconsistent whenever Yc > 0.

  Geometry (exactly consistent: dVg/dx = Ab is verified numerically)
      Ri = Ri0 + x ,  Lp = Lp0 - N_end*x
      Ab = 2*pi*Ri*Lp + N_end*pi*(Rp^2 - Ri^2)
      Vg = Vg0 + Vp0 - pi*(Rp^2 - Ri^2)*Lp

  Burn rate (Saint-Robert + temperature sensitivity + optional erosive burning,
  Lenoir-Robert form, solved implicitly):
      r0 = a*(p0/p_ref)^n * exp(sigma_p*(T_grain - T_ref))
      r  = r0 + alpha*G^0.8/D_h^0.2 * exp(-beta*rho_p*r/G)

  Nozzle (isentropic GAS-PHASE reference flow, discharge coefficient Cd, thrust
  efficiency eta_F; two-phase mass split by the chamber composition)
      gas reference flow:
        choked      :  mdot_gas_ref = Cd*At*p0*sqrt(g/(R*T0))
                                   *(2/(g+1))^((g+1)/(2(g-1)))
        unchoked    :  pe = pa, subsonic exit Mach from the pressure ratio
                       (this branch is *continuous* with the choked branch)
        separated   :  Summerfield criterion pe < f_sep*pa (only for Ae/At > 1)
      homogeneous-equilibrium mode (standard, complete entrainment):
        mdot_total     = mdot_gas_ref / Yg_chamber
        mdot_gas       = Yg_chamber * mdot_total
        mdot_condensed = Yc_chamber * mdot_total
      single-phase legacy mode:  mdot_total = mdot_gas = mdot_gas_ref (Yg = 1)
      thrust (complete velocity equilibrium, particle slip neglected):
        F = eta_F*mdot_total*ve + (pe - pa)*Ae_effective
      NOTE: with complete entrainment the momentum term can OVERESTIMATE the
      thrust of propellants with a large condensed fraction (real particles lag
      the gas); eta_F may additionally lump two-phase losses - avoid double
      correction.

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
  analytic derivative of Psi = Yg*R*T0 (two-phase mode) and of Theta = R*T0
  (single-phase legacy mode) used by the ODE pressure equation.

  Gas / condensed mass fractions are taken from the OFFICIAL CEA phase
  information: the product species list of the mixture is ordered gas species
  first, then condensed species (``EqSolver.num_gas`` / ``num_condensed``), and
  ``EqSolution.mass_fractions`` is summed over each block.  Yc is NEVER
  estimated from the molecular-weight difference M vs MW.  The legacy fcea2
  backend cannot provide the phase split from its .plt output and is therefore
  restricted to the single-phase legacy mode, with an explicit error otherwise.

--------------------------------------------------------------------------------
NUMERICS
--------------------------------------------------------------------------------
  State y = [p0, x, Rt, m_out_total, Impulse, m_gen_total] integrated with
  solve_ivp (LSODA/BDF/Radau, dense output, tight tolerances).  The integrated
  mass quadratures are TOTAL product masses (gas + condensed) in the
  homogeneous-equilibrium mode.  Quadratures for impulse and integrated masses
  are carried as ODE states, so they inherit the solver error control instead of
  relying on post-hoc trapezoidal sums.
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
      python GRIBS_v0.5.0-alpha.py --config another_config.json
  Run numerical self-tests only with:
      python GRIBS_v0.5.0-alpha.py --selftest

--------------------------------------------------------------------------------
STATUS / LIMITATIONS
--------------------------------------------------------------------------------
  * GRIBS keeps its own unsteady internal-ballistics and nozzle-flow model.
    CEA theoretical rocket performance (c*, Cf, Isp) is reported ONLY as a
    clearly-named diagnostic and never feeds the GRIBS results.
  * CEA does not provide burn-rate coefficients; the Saint-Robert inputs remain
    empirical user data.
  * The two-phase model is a homogeneous-equilibrium approximation: no particle
    slip, no particle-size distribution, no wall deposition/slag, no separate
    condensed-phase energy equation, instantaneous chemical equilibrium, and
    the HP equilibrium temperature is used as a function of pressure only.
  * The initial chamber gas is assumed to be equilibrium products of the main
    propellant; an initial fill of air or igniter gas of a different
    composition is not represented exactly by the single-composition model.
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
PROGRAM_VERSION = "0.5.0-alpha"
SCHEMA_VERSION = "0.5.0-alpha"
#: Schema string of the pre-migration configuration; used only by the explicit
#: v0.3 -> v0.4 migration helper (never for normal operation).
LEGACY_SCHEMA_VERSION = "0.3.0-alpha"
#: Schema string of the v0.4.3 single-phase configuration; used only by the
#: explicit v0.4.3 -> v0.4.4 migration helper (never for normal operation).
PRE_TWO_PHASE_SCHEMA_VERSION = "0.4.3-alpha"

#: Schema string immediately preceding the v0.5.0-alpha schema.  A v0.4.4
#: document already contains the two-phase model introduced in that release.
PRE_V050_SCHEMA_VERSION = "0.4.4-alpha"

#: Chamber two-phase model: conserved mass = total product mass (gas +
#: condensed), EOS acts on the gas phase through Psi = Yg*Rg*T0 (standard).
TWO_PHASE_HOMOGENEOUS = "homogeneous_equilibrium"
#: Chamber model identical to GRIBS v0.4.3-alpha (Yg := 1, Psi := Theta): the
#: full generated product mass feeds a gas-only EOS.  Kept ONLY as an isolated
#: regression-comparison pathway; inconsistent whenever condensed products exist.
TWO_PHASE_SINGLE_LEGACY = "single_phase_legacy"
KNOWN_TWO_PHASE_MODES = (TWO_PHASE_HOMOGENEOUS, TWO_PHASE_SINGLE_LEGACY)
#: The only condensed-volume treatment implemented (declared in the schema).
TWO_PHASE_CONDENSED_VOLUME_CHOICES = ("neglected",)
#: The only nozzle two-phase treatment implemented (declared in the schema):
#: complete entrainment - both phases leave with the same velocity.
TWO_PHASE_NOZZLE_ENTRAINMENT_CHOICES = ("complete",)

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
#: Rounding-level tolerance for phase mass fractions.  Yg or Yc may overshoot
#: [0, 1] by at most this absolute amount (double rounding of CEA values that
#: are exactly 0 or 1); such an overshoot is corrected by projection.  Anything
#: larger is a ThermochemistryError, never a silent wide clip.
_PHASE_FRACTION_ROUNDING_TOL = 1.0e-12
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

    # ---------------- two-phase chamber / nozzle model ----------------
    # 'homogeneous_equilibrium' (standard): the conserved chamber mass is the
    #   TOTAL product mass mt = mg + mc; EOS p*Vg = mt*Yg(p)*Rg(p)*T0(p); the
    #   nozzle total flow is the gas reference flow divided by Yg_chamber
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
    "cea_py_molecular_weight": ("'gas_phase_M' only: the gas-phase equation of state "
                                "must use the gas-phase molecular weight. The former "
                                "'total_MW' choice is rejected (v0.4.4): MW is tabulated "
                                "as a diagnostic and never enters the EOS or nozzle."),
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
                       "EOS on the gas phase via Psi = Yg*Rg*T0) | "
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
#      props(p)                  -> (Rg [J/(kg K)], T0 [K], gamma_s [-])
#                                   GAS-PHASE properties only
#      psi_and_derivative(p)     -> (Psi = Yg*Rg*T0 [J/kg], dPsi/dp [J/(kg Pa)])
#                                   (two-phase mode; Psi := Theta in legacy mode)
#      theta_and_derivative(p)   -> (Theta = Rg*T0 [J/kg], dTheta/dp [J/(kg Pa)])
#                                   (single-phase legacy pathway)
#      phase_fractions(p)        -> (Yg [-], Yc [-]),  Yg + Yc = 1
#                                   (1, 0) in the single-phase legacy mode
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
                        "gas_phase_molecular_weight_kg_kmol", "gas_constant_J_kgK")
#: Extra columns every cache must carry once the homogeneous-equilibrium
#: two-phase model is used (cache schema >= 3).  A cache without them cannot be
#: used for the two-phase model - it is rejected with an actionable message.
_TABLE_TWO_PHASE_KEYS = ("gas_mass_fraction", "condensed_mass_fraction")


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
    #: True when the backend can tabulate the gas/condensed mass fractions from
    #: the official CEA phase information (required by the two-phase model).
    supports_two_phase = False

    #: lower/upper validity limit of the property representation [Pa]
    pmin: float = 1.0e5
    pmax: float = 8.0e6

    def props(self, pressure_pa: float) -> Tuple[float, float, float]:
        raise NotImplementedError

    def theta_and_derivative(self, pressure_pa: float) -> Tuple[float, float]:
        raise NotImplementedError

    def psi_and_derivative(self, pressure_pa: float) -> Tuple[float, float]:
        raise NotImplementedError

    def phase_fractions(self, pressure_pa: float) -> Tuple[float, float]:
        raise NotImplementedError

    def is_extrapolated(self, pressure_pa: float) -> bool:
        raise NotImplementedError

    def metadata(self) -> Dict[str, Any]:
        raise NotImplementedError

    # -- convenience ---------------------------------------------------------
    def theta(self, pressure_pa: float) -> float:
        return self.theta_and_derivative(pressure_pa)[0]

    def psi(self, pressure_pa: float) -> float:
        return self.psi_and_derivative(pressure_pa)[0]

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
    interpolators in z = ln(p).  The derivatives of Psi = Yg*Rg*T0 (two-phase
    mode) and Theta = Rg*T0 (single-phase legacy mode) are taken from the
    analytic derivative of the interpolators and converted with

        dPsi/dp = dPsi/dln(p) / p

    so that the ODE pressure equation never sees finite-difference noise and CEA is
    never called from inside solve_ivp.

    Cache schema 3 (v0.5.0-alpha) adds the equilibrium phase split
    (gas_mass_fraction, condensed_mass_fraction) and the gas-phase molecular
    weight key name; caches from schema 2 (v0.4.3-alpha, no phase split) are
    rejected for the two-phase model and never silently reused.
    """

    cache_schema = 3
    table_stem = "thermo_table"
    #: human-readable identification of the equilibrium formulation (metadata)
    equilibrium_formulation = "unspecified"

    def __init__(self, c: Config, reactants: Sequence[ReactantSpec],
                 packing_fraction: float, cache_root: Path):
        if c.two_phase_mode == TWO_PHASE_HOMOGENEOUS and not self.supports_two_phase:
            raise ConfigurationError(
                f"two_phase_model.mode = '{TWO_PHASE_HOMOGENEOUS}' requires a "
                f"thermochemistry backend that provides the gas/condensed mass "
                f"fractions, but the selected backend '{self.name}' cannot supply "
                f"them (its CEA output does not contain the phase split, and GRIBS "
                f"never estimates it).\n"
                f"Either select backend '{BACKEND_CEA_PYTHON}', or set "
                f"two_phase_model.mode = '{TWO_PHASE_SINGLE_LEGACY}' (regression "
                f"comparison only).")
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
            # factors that change the tabulated two-phase properties:
            "phase_fractions_tabulated": bool(self.supports_two_phase),
            "molecular_weight_definition": "gas_phase_M (gas-phase M; total MW is "
                                           "diagnostic only since v0.5.0-alpha)",
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
                f"Chamber-property cache is missing {missing}: {self.cache_path}\n"
                "This looks like a cache built by an older GRIBS version (schema < 3, "
                "before the gas/condensed phase split). Old caches are never reused "
                "for the two-phase model.\n"
                "Rebuild it: set 'rebuild_cache' to true once, or delete the file, "
                "or run with --no-cache.")
        if self.c.two_phase_mode == TWO_PHASE_HOMOGENEOUS:
            missing2p = [k for k in _TABLE_TWO_PHASE_KEYS if k not in table]
            if missing2p:
                raise ThermochemistryCacheError(
                    f"Chamber-property cache is missing the two-phase columns "
                    f"{missing2p}: {self.cache_path}\n"
                    "The homogeneous-equilibrium two-phase model requires the "
                    "gas/condensed mass fractions in the cache (schema >= 3).\n"
                    "Rebuild it: set 'rebuild_cache' to true once, or delete the "
                    "file, or run with --no-cache.")
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
        M = np.asarray(table["gas_phase_molecular_weight_kg_kmol"], dtype=float)
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

        # ---- equilibrium phase split (two-phase model) -----------------------
        # Physical range: 0 < Yg <= 1 and 0 <= Yc < 1 with Yg + Yc = 1 exactly.
        # Only a rounding-level overshoot (<= _PHASE_FRACTION_ROUNDING_TOL) is
        # corrected, and it is corrected by projection onto Yg in [0, 1] with
        # Yc = 1 - Yg so the identity holds by construction.  Anything larger is
        # a ThermochemistryError, never a silent wide clip.
        two_phase = (self.c.two_phase_mode == TWO_PHASE_HOMOGENEOUS
                     and self.supports_two_phase)
        if two_phase:
            if any(k not in table for k in _TABLE_TWO_PHASE_KEYS):
                raise ThermochemistryCacheError(
                    f"Chamber-property table lacks the two-phase columns "
                    f"{list(_TABLE_TWO_PHASE_KEYS)}: rebuild the cache "
                    "(rebuild_cache = true or --no-cache).")
            Yg = np.asarray(table["gas_mass_fraction"], dtype=float)
            Yc = np.asarray(table["condensed_mass_fraction"], dtype=float)
            if Yg.shape != p.shape or Yc.shape != p.shape:
                raise ThermochemistryCacheError(
                    "Chamber-property table phase-fraction arrays have "
                    "inconsistent shapes.")
            if not (np.all(np.isfinite(Yg)) and np.all(np.isfinite(Yc))):
                raise ThermochemistryError(
                    "The chamber-property table contains non-finite phase "
                    "mass fractions.")
            if not np.allclose(Yg + Yc, 1.0, rtol=0.0, atol=1.0e-12):
                raise ThermochemistryError(
                    "The chamber-property table violates Yg + Yc = 1 beyond "
                    f"1e-12: max deviation = "
                    f"{float(np.max(np.abs(Yg + Yc - 1.0))):.3e}.")
            bad = (Yg <= 0.0) | (Yg > 1.0 + _PHASE_FRACTION_ROUNDING_TOL) | \
                  (Yc < -_PHASE_FRACTION_ROUNDING_TOL) | (Yc >= 1.0)
            if np.any(bad):
                i = int(np.argmax(bad))
                raise ThermochemistryError(
                    f"Non-physical phase mass fractions at p = {p[i]:.6g} Pa: "
                    f"Yg = {Yg[i]:.6g}, Yc = {Yc[i]:.6g}. "
                    "Require 0 < Yg <= 1 and 0 <= Yc < 1; this is not corrected "
                    "silently.")
            # rounding-level projection only (see the comment above)
            Yg = np.clip(Yg, 0.0, 1.0)
            Yc = 1.0 - Yg
            if np.any(Yg <= 0.0):
                raise ThermochemistryError(
                    "The gas-phase mass fraction reached zero: a purely "
                    "condensed equilibrium cannot be represented by this model.")
        else:
            # single-phase legacy pathway: Yg := 1 exactly (v0.4.3 behaviour)
            Yg = np.ones_like(p)
            Yc = np.zeros_like(p)

        self.pmin = float(p[0])
        self.pmax = float(p[-1])
        z = np.log(p)
        self.Ri = PchipInterpolator(z, R, extrapolate=True)
        self.Ti = PchipInterpolator(z, T, extrapolate=True)
        self.gi = PchipInterpolator(z, g, extrapolate=True)
        self.thi = PchipInterpolator(z, R * T, extrapolate=True)
        self.dthi = self.thi.derivative()
        # Psi = Yg * Rg * T_used with T_used = eta_T0 * T_CEA.  Interpolated as
        # its own PCHIP column so that the analytic derivative dPsi/dp used by
        # the pressure ODE is consistent with the interpolated Psi itself
        # (never a finite difference, never a product of separately
        # interpolated factors).
        psi = Yg * R * T
        self.psi_tab = psi
        self.psii = PchipInterpolator(z, psi, extrapolate=True)
        self.dpsii = self.psii.derivative()
        self.Ygi = PchipInterpolator(z, Yg, extrapolate=True)

    def _validate_table(self) -> Dict[str, Any]:
        """Numerical sanity report for the table (recorded in summary.json)."""
        p = np.asarray(self.table["pressure_Pa"], dtype=float)
        T_used = (np.asarray(self.table["temperature_K"], dtype=float) * self.c.eta_T0)
        R = np.asarray(self.table["gas_constant_J_kgK"], dtype=float)
        theta = R * T_used
        two_phase = (self.c.two_phase_mode == TWO_PHASE_HOMOGENEOUS
                     and self.supports_two_phase)
        if two_phase:
            Yg = np.asarray(self.table["gas_mass_fraction"], dtype=float)
            Yc = np.asarray(self.table["condensed_mass_fraction"], dtype=float)
            psi = Yg * theta
        else:
            Yg = np.ones_like(p)
            Yc = np.zeros_like(p)
            psi = theta

        # Pressure factor D(p) = 1 - p*Psi'/Psi using the ANALYTIC PCHIP
        # derivative (the same object the ODE consumes), evaluated on a dense
        # grid so knot-to-knot minima are captured.
        z = np.log(np.geomspace(p[0], p[-1], max(201, 4 * p.size)))
        psi_f = np.asarray(self.psii(z), dtype=float)
        dpsi_f = np.asarray(self.dpsii(z), dtype=float) / np.exp(z)
        factor = 1.0 - np.exp(z) * dpsi_f / psi_f
        # same check for the legacy Theta factor (metadata / regression)
        theta_f = np.asarray(self.thi(z), dtype=float)
        dtheta_f = np.asarray(self.dthi(z), dtype=float) / np.exp(z)
        factor_theta = 1.0 - np.exp(z) * dtheta_f / theta_f
        # Psi consistency with Yg*R*T at the table nodes (self-test 3 preview)
        psi_err = float(np.max(np.abs(self.psi_tab - psi) / np.maximum(psi, 1e-300))) \
            if two_phase else 0.0

        report = {
            "activation_clamped": self.c.property_policy == "clamp",
            "two_phase_model": bool(two_phase),
            "pressure_points": int(p.size),
            "pressure_min_Pa": float(p[0]),
            "pressure_max_Pa": float(p[-1]),
            "temperature_min_K": float(np.min(self.table["temperature_K"])),
            "temperature_max_K": float(np.max(self.table["temperature_K"])),
            "gamma_min": float(np.min(self.table["gamma_s"])),
            "gamma_max": float(np.max(self.table["gamma_s"])),
            "gas_phase_molecular_weight_min_kg_kmol":
                float(np.min(self.table["gas_phase_molecular_weight_kg_kmol"])),
            "gas_phase_molecular_weight_max_kg_kmol":
                float(np.max(self.table["gas_phase_molecular_weight_kg_kmol"])),
            "gas_mass_fraction_min": float(np.min(Yg)),
            "gas_mass_fraction_max": float(np.max(Yg)),
            "condensed_mass_fraction_min": float(np.min(Yc)),
            "condensed_mass_fraction_max": float(np.max(Yc)),
            "theta_min_J_kg": float(np.min(theta)),
            "theta_max_J_kg": float(np.max(theta)),
            "psi_min_J_kg": float(np.min(psi)),
            "psi_max_J_kg": float(np.max(psi)),
            "psi_vs_YgRT_max_rel_error": psi_err,
            "min_pressure_factor_1_minus_pPsi_over_Psi": float(np.min(factor)),
            "min_pressure_factor_1_minus_pTheta_over_Theta": float(np.min(factor_theta)),
            "all_points_finite": bool(np.all(np.isfinite(p))),
            "pressure_strictly_increasing": bool(np.all(np.diff(p) > 0.0)),
        }
        if report["min_pressure_factor_1_minus_pPsi_over_Psi"] <= 0.0:
            raise ThermochemistryError(
                "The chamber-pressure equation is non-physical over the configured "
                "property table: min(1 - p*Psi'/Psi) = "
                f"{report['min_pressure_factor_1_minus_pPsi_over_Psi']:.6g} <= 0.\n"
                "Check the thermochemistry inputs (reactants, temperatures, "
                "pressure range) before continuing.")
        if report["min_pressure_factor_1_minus_pTheta_over_Theta"] <= 0.0:
            raise ThermochemistryError(
                "The single-phase (Theta) pressure equation is non-physical over "
                "the configured property table: min(1 - p*Theta'/Theta) = "
                f"{report['min_pressure_factor_1_minus_pTheta_over_Theta']:.6g} <= 0.")
        if two_phase and psi_err > 1.0e-12:
            raise ThermochemistryError(
                "Cached Psi column is inconsistent with Yg*Rg*T0 (max rel. err "
                f"{psi_err:.3e} > 1e-12); rebuild the chamber-property cache.")
        return report

    # -- solver-facing API ---------------------------------------------------
    def _effective_pressure(self, p: float) -> float:
        if self.c.property_policy == "clamp":
            return min(max(p, self.pmin), self.pmax)
        return max(p, 1.0)

    def _outside_range(self, p: float) -> bool:
        return bool(p < self.pmin or p > self.pmax)

    def _clamp_flags(self, p: float) -> Tuple[float, bool]:
        """Return (effective pressure, clamped?) for a single query.

        All interpolators (R, T0, gamma_s, Yg, Psi) are evaluated at the SAME
        effective pressure, so clamping can never desynchronize the properties
        from Psi and its derivative; when the clamp is active the derivative of
        every clamped quantity is zero by definition.
        """
        pe = self._effective_pressure(p)
        return pe, (self.c.property_policy == "clamp" and p != pe)

    def _clamp_phase_fraction(self, y: float, name: str, p: float) -> float:
        """Enforce the physical range with a rounding-only correction policy."""
        if not math.isfinite(y):
            raise ThermochemistryError(
                f"Non-finite {name} at p = {p:.6g} Pa.")
        if -_PHASE_FRACTION_ROUNDING_TOL <= y < 0.0:
            return 0.0
        if 1.0 < y <= 1.0 + _PHASE_FRACTION_ROUNDING_TOL:
            return 1.0
        if not (0.0 <= y <= 1.0):
            raise ThermochemistryError(
                f"Interpolated {name} = {y:.12g} is outside [0, 1] beyond the "
                f"rounding tolerance {_PHASE_FRACTION_ROUNDING_TOL:g} at "
                f"p = {p:.6g} Pa; not corrected silently.")
        return y

    def props(self, p: float) -> Tuple[float, float, float]:
        """GAS-PHASE properties: Rg [J/(kg K)], T0 [K], gamma_s [-]."""
        pe, _ = self._clamp_flags(p)
        z = math.log(pe)
        R = float(self.Ri(z))
        T = float(self.Ti(z))
        g = float(self.gi(z))
        if not (_finite_positive(R) and _finite_positive(T) and math.isfinite(g) and g > 1.0):
            raise ThermochemistryError(
                f"Non-physical interpolated gas properties at p = {p:.6g} Pa: "
                f"R = {R:.6g}, T0 = {T:.6g}, gamma_s = {g:.6g}.")
        return R, T, g

    def phase_fractions(self, p: float) -> Tuple[float, float]:
        """Equilibrium phase mass fractions (Yg, Yc) at pressure p [Pa].

        Yc is evaluated as 1 - Yg so that Yg + Yc = 1 holds exactly for the
        interpolated state, and both values are range-checked with the
        rounding-only correction policy (no silent wide clipping).
        """
        if self.c.two_phase_mode != TWO_PHASE_HOMOGENEOUS or not self.supports_two_phase:
            return 1.0, 0.0
        pe, _ = self._clamp_flags(p)
        z = math.log(pe)
        Yg = self._clamp_phase_fraction(float(self.Ygi(z)), "Yg", p)
        if Yg <= 0.0:
            raise ThermochemistryError(
                f"Interpolated gas-phase mass fraction is zero at p = {p:.6g} Pa; "
                "a purely condensed equilibrium is outside this model.")
        return Yg, 1.0 - Yg

    def psi_and_derivative(self, p: float) -> Tuple[float, float]:
        """Psi = Yg*Rg*T0 [J/kg] and its ANALYTIC derivative dPsi/dp [J/(kg Pa)].

        In the single-phase legacy mode Psi := Theta = Rg*T0, i.e. the exact
        v0.4.3-alpha quantity.
        """
        if self.c.two_phase_mode != TWO_PHASE_HOMOGENEOUS or not self.supports_two_phase:
            return self.theta_and_derivative(p)
        pe, clamped = self._clamp_flags(p)
        z = math.log(pe)
        psi = float(self.psii(z))
        dpsi = 0.0 if clamped else float(self.dpsii(z)) / pe
        if not (math.isfinite(psi) and psi > 0.0 and math.isfinite(dpsi)):
            raise ThermochemistryError(
                f"Non-physical Psi or dPsi/dp at p = {p:.6g} Pa "
                f"(Psi = {psi:.6g} J/kg, dPsi/dp = {dpsi:.6g}).")
        return psi, dpsi

    def theta_and_derivative(self, p: float) -> Tuple[float, float]:
        pe, clamped = self._clamp_flags(p)
        z = math.log(pe)
        theta = float(self.thi(z))
        dtheta = 0.0 if clamped else float(self.dthi(z)) / pe
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
            "two_phase_model": {
                "mode": self.c.two_phase_mode,
                "condensed_volume": self.c.two_phase_condensed_volume,
                "nozzle_entrainment": self.c.two_phase_nozzle_entrainment,
                "phase_fractions_available": bool(self.supports_two_phase),
            },
            "interpolation": {
                "variable": "z = ln(p)",
                "method": "PCHIP (monotone cubic Hermite, scipy PchipInterpolator)",
                "quantities": ["Rg", "T0", "gamma_s", "Yg", "Psi = Yg*Rg*T0",
                               "Theta = Rg*T0"],
                "derivative": ("dPsi/dp = d/dz[PCHIP(Psi)](z) / p (analytic; "
                               "never a finite difference)"),
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
            "specific_gas_constant_definition": ("Rg(p) = Ru / M_gas(p); the "
                                                 "GAS-PHASE molecular weight only "
                                                 "enters the EOS and the nozzle"),
            "phase_fraction_definition": ("Yg(p) = sum of the gas-species mass "
                                          "fractions, Yc(p) = sum of the "
                                          "condensed-species mass fractions, from "
                                          "the official CEA phase ordering "
                                          "(gas species first, then condensed)"),
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

    TWO-PHASE MODEL: the legacy .plt output (``plot p t gam m``) provides only
    p, T, gamma_s and the gas-phase molecular weight - it does NOT provide the
    gas/condensed mass fractions, and GRIBS never estimates them from the
    molecular-weight difference.  This backend is therefore restricted to
    ``two_phase_model.mode = 'single_phase_legacy'``; selecting the
    homogeneous-equilibrium model with this backend is a configuration error.
    """

    name = BACKEND_CEA_LEGACY_EXECUTABLE
    uses_cea = True
    supports_two_phase = False
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
            "phase_fractions_available": False,
            "phase_fractions_note": ("the legacy .plt columns (p t gam m) do not "
                                     "provide the gas/condensed mass fractions; "
                                     "this backend is single-phase legacy only"),
            "pressure_Pa": pp.tolist(),
            "temperature_K": T.tolist(),
            "gamma_s": gamma.tolist(),
            "gas_phase_molecular_weight_kg_kmol": M.tolist(),
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
    * Molecular weight: ``EqSolution.M`` - the GAS-PHASE molecular weight in
      kg/kmol (= g/mol), i.e. the same quantity the legacy ``plot m`` column
      provided.  The v0.4.3 option ``total_MW`` is REJECTED in v0.4.4: the
      total molecular weight ``EqSolution.MW`` is tabulated as a diagnostic
      only and never enters the equation of state or the nozzle model.
    * gamma: ``EqSolution.gamma_s`` (isentropic exponent of the equilibrium gas).
    * Specific gas constant: Rg = Ru/M with Ru = 8314.51 J/(kmol K) - the value
      of ``cea.R``, verified at runtime.
    * Phase split (new in v0.4.4): ``Mixture.species_names`` lists the product
      species gas species first, then condensed species; the counts are
      ``EqSolver.num_gas`` / ``EqSolver.num_condensed``.  ``EqSolution.
      mass_fractions`` (which sums to 1 over ALL species) is partitioned
      accordingly:
          Yg = sum of the gas-species mass fractions
          Yc = sum of the condensed-species mass fractions = 1 - Yg
      Yc is NEVER estimated from the molecular-weight difference (M vs MW);
      that ratio carries no reliable condensed-mass information.
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
    supports_two_phase = True
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
        if c.cea_py_molecular_weight == "total_MW":
            raise ConfigurationError(
                "thermochemistry.cea_python.molecular_weight = 'total_MW' is no "
                "longer accepted (v0.5.0-alpha).\n"
                "  The chamber equation of state p*Vg = mt*Yg(p)*Rg(p)*T0(p) and "
                "the nozzle model must use the GAS-PHASE molecular weight; using "
                "the total molecular weight there is physically inconsistent "
                "whenever condensed products exist.\n"
                "  Set molecular_weight = 'gas_phase_M'.  The total molecular "
                "weight (EqSolution.MW) is still tabulated and reported as the "
                "diagnostic 'total_molecular_weight_kg_kmol'.")
        if c.cea_py_molecular_weight != "gas_phase_M":
            raise ConfigurationError(
                "thermochemistry.cea_python.molecular_weight must be "
                "'gas_phase_M' (the only accepted value since v0.5.0-alpha).")
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

    # -- phase split ----------------------------------------------------------
    @staticmethod
    def _phase_split(product_mixture, solver, solution) -> Tuple[float, float, Dict[str, float]]:
        """Partition ``EqSolution.mass_fractions`` into gas and condensed blocks.

        Method (official CEA phase information, verified against the installed
        package - no guessing):
          * ``Mixture.species_names`` lists every product species; the official
            ordering is gas species first, then condensed species.
          * ``EqSolver.num_gas`` / ``EqSolver.num_condensed`` give the block
            sizes; ``num_gas + num_condensed == num_products`` is asserted.
          * ``EqSolution.mass_fractions`` maps species name -> mass fraction
            (kg of species per kg of TOTAL products; the values sum to 1).
        Species are assigned to the blocks by SET MEMBERSHIP of the names, so
        the result never depends on dictionary ordering.
        """
        names = list(product_mixture.species_names)
        ng = int(solver.num_gas)
        nc = int(solver.num_condensed)
        if len(names) != ng + nc or len(names) != int(product_mixture.num_species):
            raise ThermochemistryError(
                "Inconsistent CEA species bookkeeping: "
                f"len(species_names) = {len(names)}, num_gas = {ng}, "
                f"num_condensed = {nc}, num_species = "
                f"{int(product_mixture.num_species)}.")
        gas_names = set(names[:ng])
        condensed_names = set(names[ng:])
        if gas_names & condensed_names:
            raise ThermochemistryError(
                "Overlapping gas/condensed species sets from the CEA mixture; "
                "refusing to classify the phases by guesswork.")
        mf = solution.mass_fractions
        if set(mf.keys()) != set(names):
            raise ThermochemistryError(
                "EqSolution.mass_fractions does not cover exactly the product "
                "species list of the mixture; the phase split cannot be formed "
                "reliably.")
        Yg = sum(float(mf[n]) for n in gas_names)
        Yc = sum(float(mf[n]) for n in condensed_names)
        condensed_fractions = {n: float(mf[n]) for n in sorted(condensed_names)
                               if float(mf[n]) > 0.0}
        return Yg, Yc, condensed_fractions

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
        Yg_out: List[float] = []
        Yc_out: List[float] = []
        rho_out: List[float] = []
        cp_out: List[float] = []
        a_out: List[float] = []
        h_out: List[float] = []
        species_top: Dict[str, float] = {}
        condensed_top: Dict[str, float] = {}
        n_rounding_fixes = 0

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
            M = float(solution.M)          # gas-phase molecular weight [kg/kmol]
            MW = float(solution.MW)        # total molecular weight (diagnostic)
            gamma = float(solution.gamma_s)
            if not _finite_positive(T):
                raise ThermochemistryError(
                    f"CEA returned a non-physical chamber temperature "
                    f"T = {T!r} K at p = {p_pa:.8g} Pa.")
            if not _finite_positive(M):
                raise ThermochemistryError(
                    f"CEA returned a non-physical gas-phase molecular weight "
                    f"M = {M!r} kg/kmol at p = {p_pa:.8g} Pa.")
            if not _finite_positive(MW):
                raise ThermochemistryError(
                    f"CEA returned a non-physical total molecular weight "
                    f"MW = {MW!r} kg/kmol at p = {p_pa:.8g} Pa.")
            if not (math.isfinite(gamma) and gamma > 1.0):
                raise ThermochemistryError(
                    f"CEA returned a non-physical isentropic exponent "
                    f"gamma_s = {gamma!r} at p = {p_pa:.8g} Pa.")

            # equilibrium phase split from the official CEA phase information
            Yg, Yc, condensed_fractions = self._phase_split(
                product_mixture, solver, solution)
            if not (math.isfinite(Yg) and math.isfinite(Yc)):
                raise ThermochemistryError(
                    f"Non-finite phase mass fractions at p = {p_pa:.8g} Pa: "
                    f"Yg = {Yg!r}, Yc = {Yc!r}.")
            if abs(Yg + Yc - 1.0) > 1.0e-12:
                raise ThermochemistryError(
                    f"CEA phase mass fractions do not sum to 1 at p = "
                    f"{p_pa:.8g} Pa: Yg + Yc = {Yg + Yc:.15g}.")
            # Rounding-level projection only: summing ~200 double-precision
            # species fractions can overshoot the [0, 1] endpoints by a few
            # ULP; ONLY that (<= _PHASE_FRACTION_ROUNDING_TOL) is corrected,
            # and Yc is then recomputed as 1 - Yg so the identity holds
            # exactly.  Larger violations are an error, never clipped.
            if 1.0 < Yg <= 1.0 + _PHASE_FRACTION_ROUNDING_TOL:
                n_rounding_fixes += 1
                Yg = 1.0
            elif -_PHASE_FRACTION_ROUNDING_TOL <= Yg < 0.0:
                n_rounding_fixes += 1
                Yg = 0.0
            Yc = 1.0 - Yg
            if not (0.0 < Yg <= 1.0) or not (0.0 <= Yc < 1.0):
                raise ThermochemistryError(
                    f"Non-physical phase mass fractions at p = {p_pa:.8g} Pa: "
                    f"Yg = {Yg!r}, Yc = {Yc!r} (require 0 < Yg <= 1, "
                    f"0 <= Yc < 1; the violation exceeds the rounding "
                    f"tolerance {_PHASE_FRACTION_ROUNDING_TOL:g} and is not "
                    f"corrected silently).")

            T_out.append(T)
            M_out.append(M)
            MW_out.append(MW)
            g_out.append(gamma)
            Yg_out.append(Yg)
            Yc_out.append(Yc)
            rho_out.append(float(solution.density))
            cp_out.append(float(solution.cp_eq))
            a_out.append(float(solution.sonic_velocity) if hasattr(solution, "sonic_velocity")
                         else float("nan"))
            h_out.append(float(solution.enthalpy))
            if p_pa == pressures[-1]:
                species_top = {k: float(v) for k, v in
                               sorted(solution.mass_fractions.items(),
                                      key=lambda kv: -kv[1])[:8]}
                condensed_top = condensed_fractions

        R_out = [R_UNIVERSAL / m for m in M_out]
        # Psi = Yg * Rg * (eta_T0 * T_CEA); stored for a cache round-trip
        # consistency check - it is recomputed from Yg, R, T at load time.
        psi_out = [yg * r * c.eta_T0 * t
                   for yg, r, t in zip(Yg_out, R_out, T_out)]
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
            "gas_phase_molecular_weight_kg_kmol": M_out,
            "gas_constant_J_kgK": R_out,
            "gas_mass_fraction": Yg_out,
            "condensed_mass_fraction": Yc_out,
            "psi_J_kg": psi_out,
            "phase_fraction_method": ("EqSolution.mass_fractions partitioned by the "
                                      "official species ordering (gas species first, "
                                      "then condensed; EqSolver.num_gas / "
                                      "num_condensed), summed per block"),
            "phase_fraction_rounding_corrections": int(n_rounding_fixes),
            "reactant_enthalpy_J_kg": h_c,
            "reactant_enthalpy_over_R_K": h_R,
            "reactant_enthalpy_recomputed_kJ_kg": h_out,
            "molecular_weight_definition": "gas_phase_M (EqSolution.M); EqSolution.MW "
                                           "is tabulated as a diagnostic only",
            "product_species_count": int(product_mixture.num_species),
            "num_gas": int(solver.num_gas),
            "num_condensed": int(solver.num_condensed),
            "num_elements": int(solver.num_elements),
            "diagnostics": {
                "density_kg_m3": rho_out,
                "density_note": ("EqSolution.density equals the gas-phase ideal "
                                 "density p/(Rg*T); it is NOT the two-phase "
                                 "chamber density"),
                "cp_eq_kJ_kgK": cp_out,
                "sonic_velocity_m_s": a_out,
                "total_molecular_weight_kg_kmol": MW_out,
                "top_mass_fractions_at_max_pressure": species_top,
                "condensed_species_mass_fractions_at_max_pressure": condensed_top,
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
            "molecular_weight_definition": "EqSolution.M (gas-phase, kg/kmol); "
                                           "EqSolution.MW is diagnostic only",
            "gamma_definition": "EqSolution.gamma_s (equilibrium isentropic exponent)",
            "heat_capacity_treatment": "equilibrium (cp_eq reported as a diagnostic)",
            "phase_fractions": ("Yg = sum of gas-species mass fractions, Yc = sum "
                                "of condensed-species mass fractions, from "
                                "EqSolution.mass_fractions partitioned by the "
                                "official species ordering (num_gas first, then "
                                "num_condensed); never estimated from M vs MW"),
            "condensed_species": ("included in the equilibrium and consumed by the "
                                  "two-phase model through Yg(p), Yc(p); the "
                                  "condensed volume is neglected"),
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
    """Return GAS-PHASE chamber Rg [J/(kg K)], equilibrium T0 [K], gamma_s [-]."""
    return thermochemistry().props(p0)


def psi_and_deriv(p0: float, c: Config) -> Tuple[float, float]:
    """Return Psi = Yg*Rg*T0 [J/kg] and its analytic dPsi/dp [J/(kg Pa)].

    In the single-phase legacy mode this is Theta = Rg*T0 and dTheta/dp.
    """
    return thermochemistry().psi_and_derivative(p0)


def phase_fractions(p0: float, c: Config) -> Tuple[float, float]:
    """Return the equilibrium phase mass fractions (Yg, Yc) at p0 [Pa]."""
    return thermochemistry().phase_fractions(p0)


def theta_and_deriv(p0: float, c: Config) -> Tuple[float, float]:
    """Return Theta = Rg*T0 [J/kg] and dTheta/dp [J/(kg Pa)] (legacy pathway)."""
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
# 5. NOZZLE MODEL
#    Gas-phase isentropic reference flow (unchanged from v0.3.1-alpha) plus the
#    v0.4.4 two-phase mass split.  MASS BASIS: the flow returned to the ODE is
#    mdot_total (gas + condensed), consistent with the pressure equation that
#    conserves the TOTAL product mass.
# ==============================================================================
class NozzleTwoPhaseModel:
    """Two-phase mass split of the nozzle flow (isolated for future upgrades).

    Implemented model: HOMOGENEOUS EQUILIBRIUM, COMPLETE ENTRAINMENT.
      * The gas-phase isentropic reference flow mdot_gas_ref is computed from
        the gas-phase properties Rg, T0, gamma_s (see ``_gas_reference``).
      * Both phases are assumed to leave with the same velocity, so the total
        flow is the gas flow divided by the chamber gas mass fraction:
            mdot_total = mdot_gas_ref / Yg_chamber
            mdot_gas       = Yg_chamber * mdot_total
            mdot_condensed = Yc_chamber * mdot_total
      * This is a FIRST APPROXIMATION: it is not a rigorous two-phase choking
        model (no particle slip, no two-phase sound speed, no non-equilibrium
        nozzle chemistry, no particle-size effects).  For large condensed
        fractions it tends to OVERESTIMATE the delivered mass flow and thrust.
      * In 'single_phase_legacy' mode Yg := 1 and mdot_total = mdot_gas_ref,
        i.e. exactly the v0.4.3-alpha behaviour (regression comparison only).

    A future rigorous model (e.g. frozen/non-equilibrium two-phase choking,
    particle lag) should replace ``split`` without touching the ODE, which only
    consumes (mdot_total, mdot_gas, mdot_condensed).
    """

    def __init__(self, mode: str):
        if mode not in KNOWN_TWO_PHASE_MODES:
            raise ConfigurationError(
                f"Unknown two-phase mode {mode!r}; expected one of "
                f"{list(KNOWN_TWO_PHASE_MODES)}.")
        self.mode = mode

    @property
    def entrainment_assumption(self) -> str:
        return ("complete entrainment: mdot_total = mdot_gas_ref / Yg_chamber, "
                "both phases share the exit velocity (no particle slip)"
                if self.mode == TWO_PHASE_HOMOGENEOUS else
                "single-phase legacy: the total flow equals the gas reference "
                "flow (Yg := 1)")

    def split(self, mdot_gas_reference: float,
              Yg_chamber: float) -> Tuple[float, float, float]:
        """Return (mdot_total, mdot_gas, mdot_condensed) [kg/s]."""
        if self.mode == TWO_PHASE_SINGLE_LEGACY:
            return mdot_gas_reference, mdot_gas_reference, 0.0
        if not (0.0 < Yg_chamber <= 1.0):
            raise ThermochemistryError(
                f"Nozzle two-phase split received Yg_chamber = {Yg_chamber!r}; "
                "require 0 < Yg <= 1.")
        mdot_total = mdot_gas_reference / Yg_chamber
        mdot_gas = Yg_chamber * mdot_total
        return mdot_total, mdot_gas, mdot_total - mdot_gas


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
    """Nozzle mass flow, thrust, exit conditions and flow regime.

    MASS BASIS (v0.5.0-alpha): the returned ``mdot_total`` is the TOTAL product
    mass flow (gas + condensed) and is the quantity that must be fed to the
    pressure ODE, which conserves the total product mass.  ``mdot_gas`` and
    ``mdot_condensed`` decompose it with the chamber equilibrium phase split.
    ``mdot`` is retained as a backward-compatible ALIAS of ``mdot_total`` (it is
    NOT the gas-only flow); prefer the explicit names in new code.

    The gas-phase isentropic reference flow uses the GAS-PHASE properties
    Rg(p0), T0(p0), gamma_s(p0).  In the homogeneous-equilibrium mode the total
    flow is mdot_gas_ref / Yg(p0) (complete entrainment, see
    ``NozzleTwoPhaseModel``); in the single-phase legacy mode it equals the
    gas reference flow.
    """
    split_model = NozzleTwoPhaseModel(c.two_phase_mode)
    Yg, Yc = phase_fractions(p0, c)

    if p0 <= c.p_a:
        return dict(mdot=0.0,                   # alias of mdot_total (documented)
                    mdot_total=0.0,
                    mdot_gas=0.0,
                    mdot_condensed=0.0,
                    mdot_gas_reference=0.0,
                    gas_mass_fraction=Yg,
                    condensed_mass_fraction=Yc,
                    F=0.0, Me=0.0, pe=p0, ve=0.0,
                    regime="no-flow", A_eff=Ae)

    R, T0, g = gas_props(p0, c)
    eps = Ae / At
    sq = math.sqrt(g / (R * T0))
    p_thr = choke_threshold_pressure(p0, At, Ae, c)

    if p0 >= p_thr:                                   # ---- choked ----
        mdot_ref = c.Cd * At * p0 * sq \
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
        # Complete velocity equilibrium: both phases leave at ve_eff.  With
        # particle slip the momentum term would be mdot_gas*vg + mdot_c*vc and
        # the thrust would be lower for large condensed fractions; eta_F may
        # additionally lump two-phase losses (avoid double correction).
        mdot_total, mdot_gas, mdot_cond = split_model.split(mdot_ref, Yg)
        F = c.eta_thrust * mdot_total * ve_eff + (pe_eff - c.p_a) * A_eff
        return dict(mdot=mdot_total,
                    mdot_total=mdot_total,
                    mdot_gas=mdot_gas,
                    mdot_condensed=mdot_cond,
                    mdot_gas_reference=mdot_ref,
                    gas_mass_fraction=Yg,
                    condensed_mass_fraction=Yc,
                    F=F, Me=Me, pe=pe_eff, ve=ve_eff,
                    regime=regime, A_eff=A_eff)

    # ---- subsonic (exit pressure equals ambient); continuous with the above ----
    Me = math.sqrt(2.0 / (g - 1.0) * ((p0 / c.p_a) ** ((g - 1.0) / g) - 1.0))
    mdot_ref = c.Cd * Ae * p0 * sq * Me \
        * (1.0 + 0.5 * (g - 1.0) * Me * Me) ** (-(g + 1.0) / (2.0 * (g - 1.0)))
    Te = T0 / (1.0 + 0.5 * (g - 1.0) * Me * Me)
    ve = Me * math.sqrt(g * R * Te)
    mdot_total, mdot_gas, mdot_cond = split_model.split(mdot_ref, Yg)
    F = c.eta_thrust * mdot_total * ve
    return dict(mdot=mdot_total,
                mdot_total=mdot_total,
                mdot_gas=mdot_gas,
                mdot_condensed=mdot_cond,
                mdot_gas_reference=mdot_ref,
                gas_mass_fraction=Yg,
                condensed_mass_fraction=Yc,
                F=F, Me=Me, pe=c.p_a, ve=ve,
                regime="subsonic", A_eff=Ae)


# ==============================================================================
# 6. QUASI-STEADY EQUILIBRIUM PRESSURE (diagnostic only)
# ==============================================================================
def equilibrium_pressure(Ab: float, Ri: float, At: float, Ae: float,
                         c: Config) -> Optional[float]:
    """Lowest pressure satisfying mdot_gen_total(p) = mdot_out_total(p).

    Both sides use the TOTAL product mass basis, consistent with the pressure
    ODE; None if no crossing is found.  Diagnostic only.
    """
    if not np.isfinite(Ab) or Ab <= 0.0:
        return None

    def bal(p):
        try:
            v = (c.rho_p * burn_rate(p, c, Ab, Ri) * Ab
                 - nozzle_state(p, At, Ae, c)["mdot_total"])
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
#    y = [p0, x, Rt, m_out_total, Impulse, m_gen_total]
#    The mass quadratures are TOTAL product masses (gas + condensed) in the
#    homogeneous-equilibrium mode.  Pressure equation (two-phase model B):
#      dp0/dt = [ Psi*(mdot_gen_total - mdot_out_total) - p0*Ab*r ]
#               / [ Vg*(1 - p0*Psi'/Psi) ]
#    with Psi = Yg*Rg*T0 and the ANALYTIC dPsi/dp.  In the single-phase legacy
#    mode Psi := Theta = Rg*T0 (exact v0.4.3-alpha behaviour).
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
        # TOTAL generated product mass flow (gas + condensed) [kg/s]: every
        # kilogram of burned propellant enters mt, regardless of its phase.
        m_gen_total = c.rho_p * r * Ab + igniter_mdot(t, c)

        # TOTAL nozzle outflow (gas + condensed), consistent mass basis.
        nz = nozzle_state(p0, At, Ae, c)
        psi, dpsi = psi_and_deriv(p0, c)
        denom = Vg * (1.0 - p0 * dpsi / psi)
        if denom <= 0.0:
            raise ValueError("Non-physical factor (1 - p0*Psi'/Psi) <= 0.")

        dp = (psi * (m_gen_total - nz["mdot_total"]) - p0 * Ab * r) / denom
        dRt = c.ero_throat_c * (p0 / c.p_ref) ** c.ero_throat_m
        return np.array([dp, r, dRt, nz["mdot_total"], nz["F"], m_gen_total])

    return rhs


def _event(fn, terminal=True, direction=0.0):
    fn.terminal = terminal
    fn.direction = direction
    return fn


def run_model(c: Config) -> dict:
    """Integrate the burning phase and (optionally) the blowdown phase."""
    validate(c)


    require_implemented_transition_policy(c)

    At0 = math.pi * c.R_t0 ** 2
    Ae = c.eps_nozzle * At0                      # geometric exit area, fixed
    xw = web_thickness(c)

    init_regime = nozzle_state(c.p0_init, At0, Ae, c)["regime"]
    if c.unchoked_policy == "stop" and init_regime != "choked":
        raise ValueError("Initial state is not choked: policy 'stop' has "
                         "nothing to integrate.")

    theta0, _ = theta_and_deriv(c.p0_init, c)
    psi0, _ = psi_and_deriv(c.p0_init, c)
    y0 = np.array([c.p0_init, 0.0, c.R_t0, 0.0, 0.0, 0.0])
    # Initial chamber inventory as TOTAL product mass:
    #   mt0 = p0*Vg0/Psi(p0)   (two-phase mode; Psi = Yg*Rg*T0)
    #   m0  = p0*Vg0/Theta(p0) (single-phase legacy mode)
    # ASSUMPTION: the initial chamber gas is equilibrium combustion products of
    # the main propellant at p0_init.  An initial fill of air or igniter gas of
    # a different composition is not represented exactly by this
    # single-composition model.
    m_total0 = c.p0_init * c.V_g0 / psi0
    m_gas0 = c.p0_init * c.V_g0 / theta0      # diagnostic (legacy-basis value)

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
                Ae=Ae, At0=At0, m_total0=m_total0, m_gas0=m_gas0,
                init_regime=init_regime,
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
    if c.transition_policy not in ("jump", "shock"):
        raise ValueError("transition_policy must be 'jump' or 'shock'.")


def require_implemented_transition_policy(c: Config) -> None:
    """Reject valid-but-not-yet-implemented nozzle transition models."""
    if c.transition_policy == "shock":
        raise NotImplementedError(
            "nozzle.transition_policy = 'shock' is valid configuration syntax "
            "but the internal-shock nozzle transition model is not implemented. "
            "Use transition_policy = 'jump' for the current direct-switch model."
        )

# ==============================================================================
# 8. POST-PROCESSING
#    v0.5.0-alpha: history columns use the TOTAL product mass basis
#    (m_total_eos, m_total_balance, mdot_out_total/gas/condensed) plus the
#    equilibrium phase fractions and Psi.
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
    """Dense post-processing of the ODE solution.

    Column meanings (v0.5.0-alpha two-phase model B):
      mdot_gen_total        generated TOTAL product mass flow [kg/s]
      mdot_out_total/gas/condensed  nozzle outflow decomposition [kg/s]
      m_gen, m_out          integrated TOTAL product masses [kg] (ODE states)
      m_total_eos           TOTAL product mass from the EOS: p0*Vg/Psi(p0) [kg]
                            (replaces the v0.4.3 column m_gas_eos, whose meaning
                            changed - the old name is deliberately NOT kept)
      m_total_balance       mt0 + m_gen - m_out [kg] (replaces m_gas_bal)
      m_gas_equilibrium     Yg(p0)*m_total_eos [kg] (gas part of m_total_eos)
      m_condensed           Yc(p0)*m_total_eos [kg] (condensed part)
      Y_gas, Y_condensed    equilibrium phase mass fractions [-]
      Psi                   Yg*Rg*T0 [J/kg]
    """
    Ae = res["Ae"]
    xw = web_thickness(c)
    segs = [(res["sol_burn"], "burn", 4000, True)]
    if res["sol_blow"] is not None:
        segs.append((res["sol_blow"], "blowdown", 1600, True))

    rec = {k: [] for k in
           ("t", "phase", "p0", "x", "Rt", "At", "Kn", "Ri", "Lp", "Ab", "Vg",
            "r", "mdot_gen_total", "mdot_out_total", "mdot_out_gas",
            "mdot_out_condensed", "F", "Me", "pe", "ve", "regime",
            "R", "T0", "gamma", "Y_gas", "Y_condensed", "Psi",
            "m_out", "impulse", "m_gen", "m_total_eos", "m_total_balance",
            "m_gas_equilibrium", "m_condensed", "extrap", "cstar", "CF")}

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
            Yg, Yc = phase_fractions(p0, c)
            psi, _ = psi_and_deriv(p0, c)
            nz = nozzle_state(p0, At, Ae, c)
            if phase == "burn":
                r = burn_rate(p0, c, Ab, Ri)
                m_gen_rate = c.rho_p * r * Ab + igniter_mdot(t, c)
            else:
                r, Ab, m_gen_rate = 0.0, 0.0, 0.0
            m_total_eos = p0 * Vg / psi
            m_total_bal = res["m_total0"] + float(ys[IMG, k]) - float(ys[IMO, k])

            rec["t"].append(float(t));          rec["phase"].append(phase)
            rec["p0"].append(p0);               rec["x"].append(x)
            rec["Rt"].append(Rt);               rec["At"].append(At)
            rec["Kn"].append(Ab / At);          rec["Ri"].append(Ri)
            rec["Lp"].append(Lp);               rec["Ab"].append(Ab)
            rec["Vg"].append(Vg);               rec["r"].append(r)
            rec["mdot_gen_total"].append(m_gen_rate)
            rec["mdot_out_total"].append(nz["mdot_total"])
            rec["mdot_out_gas"].append(nz["mdot_gas"])
            rec["mdot_out_condensed"].append(nz["mdot_condensed"])
            rec["F"].append(nz["F"]);           rec["Me"].append(nz["Me"])
            rec["pe"].append(nz["pe"]);         rec["ve"].append(nz["ve"])
            rec["regime"].append(nz["regime"]); rec["R"].append(R)
            rec["T0"].append(T0);               rec["gamma"].append(g)
            rec["Y_gas"].append(Yg);            rec["Y_condensed"].append(Yc)
            rec["Psi"].append(psi)
            rec["m_out"].append(float(ys[IMO, k]))
            rec["impulse"].append(float(ys[IIMP, k]))
            rec["m_gen"].append(float(ys[IMG, k]))
            rec["m_total_eos"].append(m_total_eos)
            rec["m_total_balance"].append(m_total_bal)
            rec["m_gas_equilibrium"].append(Yg * m_total_eos)
            rec["m_condensed"].append(Yc * m_total_eos)
            rec["extrap"].append(is_extrapolated(p0, c))
            rec["cstar"].append(p0 * At / nz["mdot_total"]
                                if nz["mdot_total"] > 0 else np.nan)
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


def _masked_trapezoid(
    y: np.ndarray,
    x: np.ndarray,
    mask: np.ndarray,
) -> float:
    """Integrate only intervals whose two endpoints are inside the mask."""
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    mask = np.asarray(mask, dtype=bool)

    if y.ndim != 1 or x.ndim != 1 or mask.ndim != 1:
        raise ValueError("y, x, and mask must be one-dimensional.")
    if not (y.size == x.size == mask.size):
        raise ValueError("y, x, and mask must have the same length.")
    if x.size < 2:
        return 0.0

    valid_intervals = mask[:-1] & mask[1:]
    dx = np.diff(x)
    interval_areas = 0.5 * (y[:-1] + y[1:]) * dx

    return float(np.sum(interval_areas[valid_intervals]))


def nozzle_transition_diagnostics(
    h: dict,
    c: Config,
    Ae: float,
    *,
    lower_pressure_Pa: Optional[float] = None,
    upper_pressure_Pa: Optional[float] = None,
) -> dict:
    """Quantify time and impulse spent in the jump-model transition band.

    Explicit pressure bounds are primarily intended for isolated tests. During
    normal calculations, the lower boundary is the choking threshold and the
    upper boundary is the onset of the configured flow-separation criterion.
    """
    result = {
        "policy": c.transition_policy,
        "applicable": False,
        "detected": False,
        "warning_code": None,
        "interval_count": 0,
        "duration_s": 0.0,
        "duration_fraction": 0.0,
        "impulse_Ns": 0.0,
        "impulse_fraction": 0.0,
    }

    if c.transition_policy != "jump":
        return result
    if c.eps_nozzle <= 1.0 + 1.0e-12:
        return result
    if not c.use_separation:
        return result

    t = np.asarray(h["t"], dtype=float)
    p0 = np.asarray(h["p0"], dtype=float)
    At = np.asarray(h["At"], dtype=float)
    thrust = np.asarray(h["F"], dtype=float)

    if not (t.size == p0.size == At.size == thrust.size):
        raise ValueError(
            "Transition diagnostic histories t, p0, At, and F "
            "must have the same length."
        )

    result["applicable"] = True

    if t.size < 2:
        return result

    if lower_pressure_Pa is not None or upper_pressure_Pa is not None:
        if lower_pressure_Pa is None or upper_pressure_Pa is None:
            raise ValueError(
                "lower_pressure_Pa and upper_pressure_Pa must be supplied together."
            )
        if not lower_pressure_Pa < upper_pressure_Pa:
            raise ValueError(
                "Require lower_pressure_Pa < upper_pressure_Pa."
            )

        transition_mask = (
            (p0 > float(lower_pressure_Pa))
            & (p0 < float(upper_pressure_Pa))
        )
    else:
        transition_mask = np.zeros(t.size, dtype=bool)

        for i in range(t.size):
            if not np.isfinite(p0[i]) or not np.isfinite(At[i]):
                continue
            if p0[i] <= c.p_a or At[i] <= 0.0:
                continue

            try:
                if choke_margin(p0[i], At[i], Ae, c) < 0.0:
                    continue

                _, _, gamma = gas_props(p0[i], c)
                eps = Ae / At[i]
                if eps <= 1.0 + 1.0e-12:
                    continue

                exit_mach = mach_from_area(eps, gamma, supersonic=True)
                exit_pressure = p0[i] * (
                    1.0 + 0.5 * (gamma - 1.0) * exit_mach * exit_mach
                ) ** (-gamma / (gamma - 1.0))

                transition_mask[i] = (
                    exit_pressure < c.sep_ratio * c.p_a
                )
            except (ValueError, RuntimeError, ThermochemistryError):
                continue

    starts = transition_mask & np.concatenate(
        (np.array([True]), ~transition_mask[:-1])
    )
    interval_count = int(np.count_nonzero(starts))

    duration = _masked_trapezoid(
        np.ones_like(t),
        t,
        transition_mask,
    )
    transition_impulse = _masked_trapezoid(
        thrust,
        t,
        transition_mask,
    )

    total_duration = float(t[-1] - t[0])
    total_impulse = float(_TRAPZ(thrust, t))

    detected = bool(interval_count > 0 and duration > 0.0)

    result.update(
        detected=detected,
        warning_code="W_NOZZLE_TRANSITION" if detected else None,
        interval_count=interval_count,
        duration_s=duration,
        duration_fraction=(
            duration / total_duration if total_duration > 0.0 else 0.0
        ),
        impulse_Ns=transition_impulse,
        impulse_fraction=(
            transition_impulse / total_impulse
            if abs(total_impulse) > 0.0
            else 0.0
        ),
    )
    return result


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

    # EOS mass vs balance mass (TOTAL product mass basis in the two-phase mode)
    denom = np.maximum(np.abs(h["m_total_eos"]), 1e-15)
    mass_err = float(np.max(np.abs(h["m_total_balance"] - h["m_total_eos"]) / denom))

    # two-phase statistics over the computed history
    Yg_h = np.asarray(h["Y_gas"], dtype=float)
    Yc_h = np.asarray(h["Y_condensed"], dtype=float)

    p_star = choke_limit_pressure(c, At0, Ae)
    transition = nozzle_transition_diagnostics(h, c, Ae)

    mg0 = c.rho_p * base_burn_rate(c.p0_init, c) * geometry(0.0, c)[2]
    mo0 = nozzle_state(c.p0_init, At0, Ae, c)["mdot_total"]

    backend = thermochemistry()
    pv_table = getattr(backend, "property_validation", {})
    two_phase_active = (c.two_phase_mode == TWO_PHASE_HOMOGENEOUS
                        and bool(getattr(backend, "supports_two_phase", False)))
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
        total_mass_consistency_error=mass_err,
        # ---------------- two-phase model ----------------
        two_phase_model={
            "mode": c.two_phase_mode,
            "condensed_volume": c.two_phase_condensed_volume,
            "nozzle_entrainment": c.two_phase_nozzle_entrainment,
            "active": bool(two_phase_active),
        },
        gas_mass_fraction_min=float(np.min(Yg_h)),
        gas_mass_fraction_max=float(np.max(Yg_h)),
        condensed_mass_fraction_min=float(np.min(Yc_h)),
        condensed_mass_fraction_max=float(np.max(Yc_h)),
        maximum_condensed_mass_fraction_during_run=float(np.max(Yc_h)),
        mean_condensed_mass_fraction_during_burn=float(np.mean(Yc_h[burn]))
        if np.any(burn) else float("nan"),
        condensed_mass_fraction_at_max_pressure=float(Yc_h[int(np.argmax(h["p0"]))]),
        gas_mass_fraction_table_min=float(pv_table.get("gas_mass_fraction_min", float("nan"))),
        gas_mass_fraction_table_max=float(pv_table.get("gas_mass_fraction_max", float("nan"))),
        condensed_mass_fraction_table_min=float(pv_table.get("condensed_mass_fraction_min", float("nan"))),
        condensed_mass_fraction_table_max=float(pv_table.get("condensed_mass_fraction_max", float("nan"))),
        min_pressure_factor_1_minus_pPsi_over_Psi=float(
            pv_table.get("min_pressure_factor_1_minus_pPsi_over_Psi", float("nan"))),
        nozzle_entrainment_assumption=(
            "complete entrainment: mdot_total = mdot_gas_ref / Yg_chamber; both "
            "phases leave with the same exit velocity (particle slip NOT "
            "modelled; momentum thrust may be overestimated for large "
            "condensed fractions)" if two_phase_active else
            "single-phase legacy: total flow = gas reference flow (Yg := 1)"),
        condensed_volume_assumption=(
            "condensed-phase volume neglected against the chamber free volume; "
            "the gas phase occupies Vg" if two_phase_active else
            "not applicable (single-phase legacy mode)"),
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
        nozzle_transition=transition,
        warnings=(
            [
                {
                    "code": "W_NOZZLE_TRANSITION",
                    "severity": "warning",
                    "message": (
                        "The direct-switch nozzle model traversed the "
                        "C-D nozzle transition band. Thrust and impulse in "
                        "this band are model-dependent because an internal-"
                        "shock solution is not implemented."
                    ),
                    "context": {
                        "policy": c.transition_policy,
                        "duration_s": transition["duration_s"],
                        "duration_fraction": transition["duration_fraction"],
                        "impulse_Ns": transition["impulse_Ns"],
                        "impulse_fraction": transition["impulse_fraction"],
                        "interval_count": transition["interval_count"],
                    },
                }
            ]
            if transition["warning_code"] == "W_NOZZLE_TRANSITION"
            else []
        ),
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


class _SyntheticTestBackend(PressureTableBackend):
    """Synthetic analytic property table for numerical self-tests ONLY.

    Not selectable from any configuration; used to verify:
      * the single-phase limit (Yg = 1 everywhere) reproduces the legacy
        Theta model bit-for-bit (self-test: single-phase limit), and
      * the homogeneous-equilibrium run matches the legacy run when the table
        contains no condensed phase (compatibility test).
    Constant T0, gamma_s, M and a configurable constant condensed fraction are
    sufficient because these tests compare MODELS, not chemistry.
    """

    name = "synthetic_selftest"
    uses_cea = False
    supports_two_phase = True
    table_stem = "thermo_synthetic_selftest"
    equilibrium_formulation = ("synthetic constant-property analytic table "
                               "(numerical self-tests only, never a production "
                               "backend)")

    def __init__(self, c: Config, condensed_fraction: float = 0.0):
        if not (0.0 <= condensed_fraction < 1.0):
            raise ValueError("condensed_fraction must lie in [0, 1).")
        self._ycond = float(condensed_fraction)
        # cache must never interfere with the self-tests
        c2 = Config(**{**asdict(c), "cea_cache_enabled": False,
                       "cea_rebuild_cache": False})
        super().__init__(c2, c.reactants_spec_list, c.packing_fraction,
                         Path(tempfile.gettempdir()) / "gribs_selftest_cache")

    def _build_table(self) -> Dict[str, Any]:
        c = self.c
        pmin = max(c.p_fit_min / 50.0, 1.0e3)
        pmax = c.p_fit_max * 10.0
        pressures = np.geomspace(pmin, pmax, 41)
        T0 = 3000.0
        gamma = 1.2
        M = 25.0
        R = R_UNIVERSAL / M
        Yg = 1.0 - self._ycond
        n = pressures.size
        return {
            "schema": self.cache_schema,
            "backend": self.name,
            "cache_key": self.key,
            "source": "synthetic analytic table (self-tests only)",
            "pressure_Pa": pressures.tolist(),
            "temperature_K": [T0] * n,
            "gamma_s": [gamma] * n,
            "gas_phase_molecular_weight_kg_kmol": [M] * n,
            "gas_constant_J_kgK": [R] * n,
            "gas_mass_fraction": [Yg] * n,
            "condensed_mass_fraction": [self._ycond] * n,
            "psi_J_kg": [Yg * R * T0 * c.eta_T0] * n,
        }

    def metadata(self) -> Dict[str, Any]:
        return {"backend": self.name, "synthetic_selftest_backend": True}


#: Tolerance of the single-phase-limit regression test.  With Yg = 1 the two
#: code paths evaluate bit-identical arithmetic, so any difference beyond
#: double-rounding level indicates a structural inconsistency.
SINGLE_PHASE_LIMIT_TOLERANCE = 1.0e-12


def self_tests(c: Config, verbose=True) -> bool:
    """Numerical self-tests executed explicitly with ``--selftest``.

    All v0.4.3-alpha tests are kept unchanged (T1-T7); the two-phase model
    adds T8-T15 (see the test names).  Post-simulation conservation tests run
    in :func:`self_tests_post_run` after the ODE integration.
    """
    global _THERMO
    out, ok_all = [], True

    def chk(name, ok, detail=""):
        nonlocal ok_all
        ok_all &= bool(ok)
        out.append(f"  [{'PASS' if ok else 'FAIL'}] {name:<46s} {detail}")

    backend = thermochemistry()
    two_phase_active = (c.two_phase_mode == TWO_PHASE_HOMOGENEOUS
                        and bool(getattr(backend, "supports_two_phase", False)))

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

    # T4 choking-boundary continuity and transition classification.
    #
    # For a converging nozzle (Ae/At = 1), mass flow and thrust must both
    # remain continuous across the choking boundary.
    #
    # For a C-D nozzle (Ae/At > 1), the current jump model switches directly
    # between supersonic choked and subsonic solutions. Mass flow must remain
    # continuous. The thrust jump remains a diagnostic until an internal-shock
    # transition model is implemented.
    continuity_tolerance = 1.0e-6
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
        chk("nozzle choking-boundary continuity", False,
            "choking-limit pressure could not be determined")
    else:
        lo = nozzle_state(pstar * (1.0 - 1e-9), At0, Ae, cc)
        hi = nozzle_state(pstar * (1.0 + 1e-9), At0, Ae, cc)
        worst_mdot = abs(lo["mdot"] / hi["mdot"] - 1.0)
        worst_thrust = abs((lo["F"] + 1e-12) / (hi["F"] + 1e-12) - 1.0)

        if math.isclose(eps, 1.0, rel_tol=0.0, abs_tol=1.0e-12):
            worst = max(worst_mdot, worst_thrust)
            chk("converging-nozzle continuity (mdot, F)",
                worst < continuity_tolerance,
                f"Ae/At = {eps:.6g}, mdot jump = {worst_mdot:.2e}, "
                f"F jump = {worst_thrust:.2e}")
        else:
            chk("C-D nozzle transition mass-flow continuity",
                worst_mdot < continuity_tolerance,
                f"Ae/At = {eps:.6g}, mdot jump = {worst_mdot:.2e}, "
                f"F jump = {worst_thrust:.2e} "
                "(diagnostic; known jump-model limitation)")

        # T5 closed-form check of the converging-nozzle choked thrust.
        # Two-phase momentum basis: only the momentum term scales with the
        # entrained condensed mass (mdot_total = mdot_gas/Yg), the pressure
        # thrust does not:
        #   F = mdot_total*ve + (pe - pa)*At
        #     = At*(crit*p*(g/Yg + 1) - pa),  crit = (2/(g+1))^(g/(g-1))
        # With Yg = 1 this reduces to the classical single-phase identity
        # F = At*((g+1)*crit*p - pa).
        cc5 = Config(**{**asdict(c), "eps_nozzle": 1.0, "Cd": 1.0,
                        "eta_thrust": 1.0, "use_separation": False})
        At0 = math.pi * cc5.R_t0 ** 2
        e = 0.0
        for p in np.geomspace(cc5.p_a / 0.5, 50e6, 20):
            nz = nozzle_state(p, At0, At0, cc5)
            _, _, g = gas_props(p, cc5)
            Yg_p, _ = phase_fractions(p, cc5)
            F_ref = At0 * (critical_ratio(g) * p * (g / Yg_p + 1.0) - cc5.p_a)
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
            float(pv.get("min_pressure_factor_1_minus_pTheta_over_Theta", -1.0)) > 0.0 and \
            float(pv.get("min_pressure_factor_1_minus_pPsi_over_Psi", -1.0)) > 0.0
        detail = (f"{pv.get('pressure_points', 'n/a')} points, "
                  f"T0 = {pv.get('temperature_min_K', float('nan')):.1f}-"
                  f"{pv.get('temperature_max_K', float('nan')):.1f} K, "
                  f"gamma = {pv.get('gamma_min', float('nan')):.4f}-"
                  f"{pv.get('gamma_max', float('nan')):.4f}")
        chk("chamber-property table physical", ok, detail)
    except Exception as exc:                                  # pragma: no cover
        chk("chamber-property table physical", False, str(exc))

    # ---- two-phase tests (new in v0.5.0-alpha) ------------------------------
    table = getattr(backend, "table", {})
    if "gas_mass_fraction" in table:
        Yg_tab = np.asarray(table["gas_mass_fraction"], dtype=float)
        Yc_tab = np.asarray(table["condensed_mass_fraction"], dtype=float)
    else:
        npts = len(table.get("pressure_Pa", [0] * 4))
        Yg_tab = np.ones(npts)
        Yc_tab = np.zeros(npts)

    # T8 phase-fraction completeness: Yg + Yc = 1 at every table point
    e = float(np.max(np.abs(Yg_tab + Yc_tab - 1.0)))
    chk("phase fractions sum to 1 (table)", e <= 1.0e-12,
        f"max |Yg+Yc-1| = {e:.2e} (tol 1e-12)")

    # T9 physical ranges: 0 < Yg <= 1, 0 <= Yc < 1
    ok = bool(np.all(Yg_tab > 0.0) and np.all(Yg_tab <= 1.0)
              and np.all(Yc_tab >= 0.0) and np.all(Yc_tab < 1.0))
    chk("phase fractions physical (table)", ok,
        f"Yg in [{Yg_tab.min():.6g}, {Yg_tab.max():.6g}], "
        f"Yc in [{Yc_tab.min():.6g}, {Yc_tab.max():.6g}]")

    # T10 Psi consistency: Psi = Yg*Rg*T_used at every table point.
    # The interpolator check is mode-consistent: in the two-phase mode the
    # solver-facing Psi is Yg*Rg*T_used; in the single-phase legacy mode it is
    # Theta = Rg*T_used by definition (Yg := 1), while the stored table column
    # keeps the physical Yg*Rg*T_used product in both cases.
    R_tab = np.asarray(table["gas_constant_J_kgK"], dtype=float)
    T_tab = np.asarray(table["temperature_K"], dtype=float) * c.eta_T0
    psi_ref = Yg_tab * R_tab * T_tab
    psi_solver_ref = psi_ref if two_phase_active else R_tab * T_tab
    if "psi_J_kg" in table:
        psi_tab = np.asarray(table["psi_J_kg"], dtype=float)
    else:
        psi_tab = psi_ref
    e = float(np.max(np.abs(psi_tab - psi_ref) / np.maximum(psi_ref, 1e-300)))
    # also check the interpolated Psi against the mode-consistent product
    e_node = 0.0
    for p, ref in zip(table["pressure_Pa"], psi_solver_ref):
        e_node = max(e_node, abs(psi_and_deriv(float(p), c)[0] / ref - 1.0))
    chk("Psi = Yg*Rg*T0 (table + interpolator)",
        e <= 1.0e-12 and e_node <= 1.0e-10,
        f"stored-vs-product rel. err = {e:.2e}, node interp err = {e_node:.2e}")

    # T11 analytic dPsi/dp against a centred difference
    e = 0.0
    for p in np.geomspace(max(c.p_fit_min * 1.05, 1.1e5), c.p_fit_max * 0.95, 25):
        hstep = p * 1e-5
        fd = (psi_and_deriv(p + hstep, c)[0] - psi_and_deriv(p - hstep, c)[0]) / (2 * hstep)
        dp = psi_and_deriv(p, c)[1]
        if dp == 0.0:
            e = max(e, abs(fd) / max(psi_and_deriv(p, c)[0], 1e-300) * p)
        else:
            e = max(e, abs(fd / dp - 1.0))
    chk("analytic vs numerical dPsi/dp", e < tol_deriv,
        f"max rel. err = {e:.2e} (tol {tol_deriv:.0e})")

    # T12 pressure-equation denominator D(p) = 1 - p*Psi'/Psi > 0 everywhere
    try:
        dmin = float(backend.property_validation[
            "min_pressure_factor_1_minus_pPsi_over_Psi"])
        chk("pressure factor 1 - p*Psi'/Psi > 0", dmin > 0.0,
            f"min D(p) = {dmin:.6f}")
    except Exception as exc:                                  # pragma: no cover
        chk("pressure factor 1 - p*Psi'/Psi > 0", False, str(exc))

    # T13 single-phase limit: with Yg = 1 everywhere the two-phase model B must
    # reduce EXACTLY to the legacy Theta model (regression test).
    # T14 nozzle decomposition identities (with a condensed phase present).
    # T15 compatibility: no condensed phase -> homogeneous and legacy runs agree.
    saved_thermo = _THERMO
    try:
        synth = _SyntheticTestBackend(c, condensed_fraction=0.0)
        _THERMO = synth
        e_psi, e_dpsi, e_rhs = 0.0, 0.0, 0.0
        c_hom = Config(**{**asdict(c), "two_phase_mode": TWO_PHASE_HOMOGENEOUS})
        c_leg = Config(**{**asdict(c), "two_phase_mode": TWO_PHASE_SINGLE_LEGACY})
        for p in np.geomspace(synth.pmin * 1.01, synth.pmax * 0.99, 30):
            psi_h, dpsi_h = psi_and_deriv(float(p), c_hom)
            th_l, dth_l = theta_and_deriv(float(p), c_leg)
            e_psi = max(e_psi, abs(psi_h / th_l - 1.0))
            if dth_l != 0.0:
                e_dpsi = max(e_dpsi, abs(dpsi_h / dth_l - 1.0))
            else:
                e_dpsi = max(e_dpsi, abs(dpsi_h - dth_l))
            y0 = np.array([float(p), 0.3 * web_thickness(c_hom), c.R_t0,
                           1e-5, 1e-4, 2e-5])
            dy_h = make_rhs(c_hom, Ae, True)(0.0, y0)
            dy_l = make_rhs(c_leg, Ae, True)(0.0, y0)
            e_rhs = max(e_rhs, float(np.max(np.abs(dy_h - dy_l)
                                            / np.maximum(np.abs(dy_l), 1e-300))))
        chk("single-phase limit: model B == Theta model",
            max(e_psi, e_dpsi, e_rhs) < SINGLE_PHASE_LIMIT_TOLERANCE,
            f"Psi err = {e_psi:.2e}, dPsi err = {e_dpsi:.2e}, "
            f"rhs err = {e_rhs:.2e} (tol {SINGLE_PHASE_LIMIT_TOLERANCE:.0e})")

        # T14 nozzle flow decomposition: mdot_total = mdot_gas + mdot_condensed,
        # mdot_gas = Yg*mdot_total, mdot_condensed = Yc*mdot_total.
        synth_c = _SyntheticTestBackend(c, condensed_fraction=0.3)
        _THERMO = synth_c
        e_dec = 0.0
        c_hom_c = Config(**{**asdict(c), "two_phase_mode": TWO_PHASE_HOMOGENEOUS})
        for p in np.geomspace(c.p_a * 1.05, c.p_fit_max, 30):
            nz = nozzle_state(float(p), At0, Ae, c_hom_c)
            if nz["mdot_total"] <= 0.0:
                continue
            Yg_p, Yc_p = phase_fractions(float(p), c_hom_c)
            e_dec = max(e_dec,
                        abs(nz["mdot_gas"] + nz["mdot_condensed"]
                            - nz["mdot_total"]) / nz["mdot_total"],
                        abs(nz["mdot_gas"] / nz["mdot_total"] - Yg_p),
                        abs(nz["mdot_condensed"] / nz["mdot_total"] - Yc_p),
                        abs(nz["mdot_total"] * Yg_p
                            - nz["mdot_gas_reference"]) / nz["mdot_gas_reference"])
        chk("nozzle decomposition mdot_total = gas + condensed",
            e_dec < 1.0e-12, f"max rel. err = {e_dec:.2e}")

        # T15 compatibility regression: with a condensed-FREE table the
        # homogeneous-equilibrium run must match the legacy run (short motor).
        _THERMO = synth
        fast = {**asdict(c), "a_burn": 10.0 * c.a_burn, "blowdown": False}
        c_hom_f = Config(**{**fast, "two_phase_mode": TWO_PHASE_HOMOGENEOUS})
        c_leg_f = Config(**{**fast, "two_phase_mode": TWO_PHASE_SINGLE_LEGACY})
        res_h = run_model(c_hom_f)
        res_l = run_model(c_leg_f)
        hist_h = sample(res_h, c_hom_f)
        hist_l = sample(res_l, c_leg_f)
        s_h = summarize(res_h, hist_h, c_hom_f)
        s_l = summarize(res_l, hist_l, c_leg_f)
        e_run = max(abs(s_h["p_max_Pa"] / s_l["p_max_Pa"] - 1.0),
                    abs(s_h["burn_time_s"] / s_l["burn_time_s"] - 1.0),
                    abs(s_h["impulse_total_Ns"] / s_l["impulse_total_Ns"] - 1.0),
                    abs(s_h["Isp_s"] / s_l["Isp_s"] - 1.0))
        t_cmp = np.linspace(0.0, 0.98 * min(s_h["burn_time_s"], s_l["burn_time_s"]), 200)
        p_h = res_h["sol_burn"].sol(t_cmp)[IP]
        p_l = res_l["sol_burn"].sol(t_cmp)[IP]
        e_traj = float(np.max(np.abs(p_h / p_l - 1.0)))
        e_run = max(e_run, e_traj)
        chk("no-condensed compatibility (hom. vs legacy run)", e_run < 1.0e-9,
            f"summary/trajectory rel. err = {e_run:.2e} (tol 1e-9)")
    except Exception as exc:                                  # pragma: no cover
        chk("two-phase regression tests", False, f"{type(exc).__name__}: {exc}")
    finally:
        _THERMO = saved_thermo

    if verbose:
        print(f"Self-tests  (thermochemistry backend: {backend.name}, "
              f"two-phase mode: {c.two_phase_mode}"
              f"{'' if two_phase_active else ' (phase split inactive)'})")
        print("\n".join(out))
        print(f"  -> {'all checks passed' if ok_all else 'FAILURES DETECTED'}\n")
    return ok_all


#: Tolerance of the post-simulation mass-conservation tests.  The integrated
#: mass quadratures inherit the ODE error control (rtol ~ 1e-9), so a relative
#: consistency of 1e-6 is far above the numerical noise floor while still
#: catching structural mass-basis errors (which are O(1) or O(Yc)).
POST_RUN_MASS_TOLERANCE = 1.0e-6


def self_tests_post_run(res: dict, h: dict, c: Config, verbose=True) -> bool:
    """Conservation tests that require the completed simulation (tests 7-8)."""
    out, ok_all = [], True

    def chk(name, ok, detail=""):
        nonlocal ok_all
        ok_all &= bool(ok)
        out.append(f"  [{'PASS' if ok else 'FAIL'}] {name:<46s} {detail}")

    m_total0 = float(res["m_total0"])
    m_gen_final = float(h["m_gen"][-1])
    m_out_final = float(h["m_out"][-1])
    m_eos_final = float(h["m_total_eos"][-1])

    # Test 7: total mass conservation  mt0 + m_gen = m_out + mt_remaining
    lhs = m_total0 + m_gen_final
    rhs = m_out_final + m_eos_final
    e7 = abs(lhs - rhs) / max(abs(lhs), 1e-300)
    chk("total mass conservation (final)", e7 < POST_RUN_MASS_TOLERANCE,
        f"mt0 + m_gen = {lhs:.9e} kg vs m_out + mt = {rhs:.9e} kg "
        f"(rel. err = {e7:.2e})")

    # Test 8: EOS mass vs balance mass over the whole history
    denom = np.maximum(np.abs(np.asarray(h["m_total_eos"], dtype=float)), 1e-15)
    e8 = float(np.max(np.abs(np.asarray(h["m_total_balance"], dtype=float)
                             - np.asarray(h["m_total_eos"], dtype=float)) / denom))
    chk("EOS mass vs balance mass (history)", e8 < POST_RUN_MASS_TOLERANCE,
        f"max rel. err = {e8:.2e} (tol {POST_RUN_MASS_TOLERANCE:.0e})")

    if verbose:
        print("Post-run conservation tests")
        print("\n".join(out))
        print(f"  -> {'all checks passed' if ok_all else 'FAILURES DETECTED'}\n")
    return ok_all


# ==============================================================================
# 10. PLOTTING (one single, carefully designed figure)
#     Panel layout unchanged from v0.3.1-alpha; v0.4.4 adds the condensed
#     mass fraction as a second axis of panel (a) and the condensed outflow
#     component to panel (c), and the key-results block reports the two-phase
#     model statistics.
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
    # condensed-phase mass fraction on a second axis (two-phase model B);
    # drawn only when a condensed phase is actually present, so the panel is
    # unchanged for condensed-free propellants
    yc_hist = np.asarray(h["Y_condensed"], dtype=float)
    if float(np.max(yc_hist)) > 0.0:
        ax_yc = ax_p.twinx()
        ax_yc.grid(False)
        ax_yc.plot(t, yc_hist, color="#8d6e63", lw=1.4, ls="-.",
                   label=r"condensed fraction $Y_c$")
        ax_yc.set_ylim(0.0, max(float(np.max(yc_hist)) * 1.35, 1e-6))
        ax_yc.set_ylabel(r"condensed mass fraction  $Y_c$  [-]", color="#8d6e63")
        ax_yc.tick_params(axis="y", colors="#8d6e63")
        ax_yc.spines["right"].set_visible(True)
        hl = ax_p.get_legend_handles_labels()
        h2 = ax_yc.get_legend_handles_labels()
        ax_p.legend(hl[0] + h2[0], hl[1] + h2[1], loc="lower right", ncol=1)
    else:
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
        f"  two-phase model    {s['two_phase_model']['mode']}",
        f"{'max Yc (run)':<24}{s['maximum_condensed_mass_fraction_during_run']:>11.4f} -",
        f"{'mean Yc (burn)':<24}{s['mean_condensed_mass_fraction_during_burn']:>11.4f} -",
        "",
        "VERIFICATION",
        f"{'propellant mass balance':<24}{s['propellant_mass_balance_error']:>11.2e}",
        f"{'total mass consistency':<24}{s['total_mass_consistency_error']:>11.2e}",
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
    pos = h["mdot_gen_total"] > 0
    ax_md.plot(t[pos], h["mdot_gen_total"][pos] * 1e3, color=C_GEN, lw=1.8,
               label=r"generated (total)  $\dot{m}_{gen}=\rho_p A_b r$")
    ax_md.plot(t, np.maximum(h["mdot_out_total"], 1e-12) * 1e3, color=C_OUT, lw=1.6,
               ls="--", label=r"nozzle (total)  $\dot{m}_{out}$")
    if np.any(np.asarray(h["mdot_out_condensed"]) > 0.0):
        ax_md.plot(t, np.maximum(h["mdot_out_condensed"], 1e-12) * 1e3,
                   color="#8d6e63", lw=1.2, ls=":",
                   label=r"nozzle condensed  $\dot{m}_{out,c}$")
    acc = h["mdot_gen_total"] > h["mdot_out_total"]
    ax_md.fill_between(t, np.maximum(h["mdot_out_total"], 1e-12) * 1e3,
                       np.maximum(h["mdot_gen_total"], 1e-12) * 1e3, where=acc,
                       color=C_ACC, alpha=0.25, lw=0, label="chamber filling")
    ax_md.axvline(t_bo, color="#404040", lw=1.0, ls="--")
    ax_md.set_yscale("log")
    ax_md.set_xlim(0.0, t[-1])
    ax_md.set_xlabel("time  $t$  [s]")
    ax_md.set_ylabel(r"mass flow (total)  [g/s]")
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
            f"two-phase: {s['two_phase_model']['mode']}"
            f"{' (complete entrainment, no particle slip)' if s['two_phase_model']['active'] else ''}   |   "
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
#     v0.5.0-alpha: the mass columns are renamed to their true meaning (total
#     product mass basis) - see the CSV_COLUMNS note below.
# ==============================================================================
#: v0.5.0-alpha column set.  MEANING CHANGES vs v0.4.3-alpha are RENAMED, not
#: silently reused:  mdot_gen -> mdot_gen_total, mdot_out -> mdot_out_total,
#: m_gas_eos -> m_total_eos, m_gas_bal -> m_total_balance (TOTAL product mass
#: basis: gas + condensed).  m_gen / m_out are the integrated TOTAL masses.
CSV_COLUMNS = ["t", "phase", "p0", "x", "Ri", "Lp", "Ab", "Vg", "Kn", "Rt", "At",
               "r", "mdot_gen_total", "mdot_out_total", "mdot_out_gas",
               "mdot_out_condensed", "m_gen", "m_out", "m_total_eos",
               "m_total_balance", "m_gas_equilibrium", "m_condensed",
               "Y_gas", "Y_condensed", "Psi",
               "regime", "Me", "pe", "ve", "F", "impulse",
               "cstar", "CF", "R", "T0", "gamma", "extrap"]
CSV_UNITS = {"t": "s", "p0": "Pa", "x": "m", "Ri": "m", "Lp": "m", "Ab": "m2",
             "Vg": "m3", "Kn": "-", "Rt": "m", "At": "m2", "r": "m/s",
             "mdot_gen_total": "kg/s", "mdot_out_total": "kg/s",
             "mdot_out_gas": "kg/s", "mdot_out_condensed": "kg/s",
             "m_gen": "kg", "m_out": "kg",
             "m_total_eos": "kg", "m_total_balance": "kg",
             "m_gas_equilibrium": "kg", "m_condensed": "kg",
             "Y_gas": "-", "Y_condensed": "-", "Psi": "J/kg",
             "Me": "-", "pe": "Pa", "ve": "m/s", "F": "N", "impulse": "N.s",
             "cstar": "m/s", "CF": "-",
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


def nozzle_transition_summary_lines(
    transition: dict,
    warnings: Sequence[dict],
) -> List[str]:
    """Format the nozzle-transition diagnostic for the text summary."""
    lines = [
        "[ nozzle transition diagnostic ]",
        f"  policy                            : {transition['policy']}",
        f"  applicable                        : {transition['applicable']}",
        f"  transition band detected          : {transition['detected']}",
        f"  contiguous intervals              : {transition['interval_count']}",
        (
            "  residence time                    : "
            f"{1.0e3 * transition['duration_s']:.4f} ms"
        ),
        (
            "  fraction of computed duration     : "
            f"{100.0 * transition['duration_fraction']:.6f} %"
        ),
        (
            "  model impulse in transition band  : "
            f"{transition['impulse_Ns']:.6f} N.s"
        ),
        (
            "  fraction of total model impulse   : "
            f"{100.0 * transition['impulse_fraction']:.6f} %"
        ),
    ]

    for warning in warnings:
        if warning.get("code") == "W_NOZZLE_TRANSITION":
            lines.append(f"  WARNING [{warning['code']}]")
            lines.append(f"    {warning['message']}")

    return lines


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
    A("[ two-phase model (homogeneous equilibrium, model B) ]")
    tpm = s["two_phase_model"]
    A(f"  model                            : {tpm['mode']}"
      f"{'  (ACTIVE)' if tpm['active'] else '  (phase split inactive: Yg := 1)'}")
    A(f"  chamber EOS                      : p0*Vg = mt*Yg(p)*Rg(p)*T0(p) = mt*Psi(p)")
    A(f"  conserved chamber mass           : TOTAL product mass mt = m_gas + m_condensed")
    A(f"  max condensed fraction (run)     : "
      f"{s['maximum_condensed_mass_fraction_during_run']:.6f} [-]")
    A(f"  mean condensed fraction (burn)   : "
      f"{s['mean_condensed_mass_fraction_during_burn']:.6f} [-]")
    A(f"  condensed fraction at MEOP       : "
      f"{s['condensed_mass_fraction_at_max_pressure']:.6f} [-]")
    A(f"  Yc over the property table       : "
      f"{s['condensed_mass_fraction_table_min']:.6f} - "
      f"{s['condensed_mass_fraction_table_max']:.6f} [-]")
    A(f"  min (1 - p*Psi'/Psi) on table    : "
      f"{s['min_pressure_factor_1_minus_pPsi_over_Psi']:.6f} [-]  (must be > 0)")
    A(f"  nozzle entrainment assumption    : {s['nozzle_entrainment_assumption']}")
    A(f"  condensed volume assumption      : {s['condensed_volume_assumption']}")
    A("  NOTE: particle slip is NOT modelled; the complete-entrainment momentum")
    A("        term can overestimate the thrust for propellants with a large")
    A("        condensed fraction.  eta_F may additionally lump two-phase losses.")
    A("        These results are NOT experimentally validated - independent")
    A("        validation is required before any engineering use.")
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
    for line in nozzle_transition_summary_lines(
        s["nozzle_transition"],
        s["warnings"],
    ):
        A(line)
    A("")
    A("[ verification ]")
    A(f"  propellant mass balance error    : {s['propellant_mass_balance_error']:.3e}")
    A(f"  total-mass (EOS vs balance) error: {s['total_mass_consistency_error']:.3e}")
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
    new["schema_version"] = PRE_TWO_PHASE_SCHEMA_VERSION
    notes.append(f"schema_version: {LEGACY_SCHEMA_VERSION!r} -> "
                 f"{PRE_TWO_PHASE_SCHEMA_VERSION!r}.")
    return new, notes


#: Default two-phase block inserted by the v0.4.3 -> v0.4.4 migration.
DEFAULT_TWO_PHASE_BLOCK: Dict[str, Any] = {
    "mode": TWO_PHASE_HOMOGENEOUS,
    "condensed_volume": "neglected",
    "nozzle_entrainment": "complete",
}


def migrate_v043_document(doc: dict) -> Tuple[dict, List[str]]:
    """Explicitly migrate a v0.4.3-alpha configuration to the v0.4.4 schema.

    The physical model CHANGES with this migration (single-phase gas-only EOS
    -> homogeneous-equilibrium two-phase model with the total product mass
    conserved), so the inserted block and its consequences are recorded as
    notes, printed as warnings and stored in summary.json.  Nothing is
    reinterpreted silently.
    """
    notes: List[str] = []
    new = copy.deepcopy(_mapping(doc, "root"))
    th = _mapping(new.get("thermochemistry", {}), "thermochemistry")

    # total_MW in the gas-phase EOS is physically inconsistent (v0.4.4 policy):
    # the migration refuses to rewrite it silently.
    cea_py = th.get("cea_python")
    if isinstance(cea_py, dict) and cea_py.get("molecular_weight") == "total_MW":
        raise ConfigurationError(
            "thermochemistry.cea_python.molecular_weight = 'total_MW' cannot be "
            "migrated automatically: since v0.5.0-alpha the gas-phase equation "
            "of state and the nozzle model use the GAS-PHASE molecular weight "
            "only, and a silent substitution is forbidden.\n"
            "Edit the configuration: set molecular_weight = 'gas_phase_M' "
            "(the total molecular weight remains available as a diagnostic).")

    if "two_phase_model" not in new:
        backend = _mapping(new.get("thermochemistry", {}),
                           "thermochemistry").get("backend")
        if backend == BACKEND_CEA_LEGACY_EXECUTABLE:
            # the fcea2 backend cannot provide the phase split; the only
            # migration that preserves the v0.4.3 behaviour is the legacy mode
            mode = TWO_PHASE_SINGLE_LEGACY
            notes.append(
                "two_phase_model: not present in the v0.4.3 document; the "
                "defaults were added with mode = 'single_phase_legacy' because "
                "the selected backend 'cea_legacy_executable' cannot provide "
                "the gas/condensed mass fractions (its .plt output does not "
                "contain them).  This exactly preserves the v0.4.3 physical "
                "model.  Select backend 'cea_python' and mode "
                "'homogeneous_equilibrium' to enable the two-phase model.")
        else:
            mode = TWO_PHASE_HOMOGENEOUS
            notes.append(
                "two_phase_model: not present in the v0.4.3 document; defaults "
                "added with mode = 'homogeneous_equilibrium' (the standard "
                "two-phase model B: the conserved chamber mass is the TOTAL "
                "product mass, EOS p*Vg = mt*Yg(p)*Rg(p)*T0(p), nozzle flow with "
                "complete entrainment).  THIS CHANGES THE PHYSICAL MODEL vs "
                "v0.4.3 whenever condensed products exist.  To reproduce exact "
                "v0.4.3 results for regression comparison, set "
                "two_phase_model.mode = 'single_phase_legacy'.")
        new["two_phase_model"] = {
            "mode": mode,
            "condensed_volume": "neglected",
            "nozzle_entrainment": "complete",
        }
    new["schema_version"] = PRE_V050_SCHEMA_VERSION
    notes.append(f"schema_version: {PRE_TWO_PHASE_SCHEMA_VERSION!r} -> "
                 f"{PRE_V050_SCHEMA_VERSION!r}.")
    return new, notes


def migrate_v044_document(doc: dict) -> Tuple[dict, List[str]]:
    """Migrate a v0.4.4-alpha document to the v0.5.0-alpha schema.

    The migration preserves every physical and numerical input.  Only the
    declared schema version changes because A-0 aligns the program identity,
    file name, and configuration schema without modifying the model.
    """
    new = copy.deepcopy(_mapping(doc, "root"))
    new["schema_version"] = SCHEMA_VERSION
    notes = [
        f"schema_version: {PRE_V050_SCHEMA_VERSION!r} -> "
        f"{SCHEMA_VERSION!r}; no physical or numerical setting changed."
    ]
    return new, notes


@dataclass
class Configuration:
    config: Config
    output_dir: Path
    output_names: Dict[str, str]
    document: dict
    migration_notes: List[str]


_ROOT_KEYS = ("schema_version", "notes", "propellant", "grain", "nozzle",
              "environment", "igniter", "thermochemistry", "two_phase_model",
              "solver", "output")
_TWO_PHASE_KEYS = ("mode", "condensed_volume", "nozzle_entrainment")
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
_NOZZLE_KEYS = (
    "initial_throat_radius_m",
    "expansion_ratio",
    "discharge_coefficient",
    "thrust_efficiency",
    "transition_policy",
    "flow_separation",
    "throat_erosion",
)
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
                        f"automatically to {PRE_TWO_PHASE_SCHEMA_VERSION}.")
        schema = doc.get("schema_version")
    if schema == PRE_TWO_PHASE_SCHEMA_VERSION:
        doc, notes44 = migrate_v043_document(doc)
        notes.extend(notes44)
        migrated_mode = _mapping(doc.get("two_phase_model", {}),
                                 "two_phase_model").get("mode")
        model_note = ("the chamber model changes to the homogeneous-equilibrium "
                      "two-phase model; see the notes below"
                      if migrated_mode == TWO_PHASE_HOMOGENEOUS else
                      "the chamber model stays on the v0.4.3 single-phase "
                      "pathway (two_phase_model.mode = 'single_phase_legacy'); "
                      "see the notes below")
        notes.insert(0, f"{path.name}: {PRE_TWO_PHASE_SCHEMA_VERSION} configuration "
                        f"migrated automatically to {PRE_V050_SCHEMA_VERSION} "
                        f"({model_note}).")
        schema = doc.get("schema_version")
    if schema == PRE_V050_SCHEMA_VERSION:
        doc, notes50 = migrate_v044_document(doc)
        notes.extend(notes50)
        notes.append(
            f"{path.name}: {PRE_V050_SCHEMA_VERSION} configuration migrated "
            f"automatically to {SCHEMA_VERSION}; physical and numerical "
            "settings were preserved."
        )
        schema = doc.get("schema_version")
    if schema != SCHEMA_VERSION:
        raise ConfigurationError(
            f"schema_version: expected {SCHEMA_VERSION!r} (or a migratable legacy "
            f"value: {PRE_V050_SCHEMA_VERSION!r}, "
            f"{PRE_TWO_PHASE_SCHEMA_VERSION!r}, or {LEGACY_SCHEMA_VERSION!r}), "
            f"got {schema!r}.\n"
            "Update the configuration file; see the migration documentation "
            "and CHANGELOG.md.")
    return doc, notes


def configuration_from_document(doc: dict) -> Configuration:
    """Validate a v0.5 document and build the runtime configuration."""
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
    tpm = _section(root, "two_phase_model", "root", required=False)
    if tpm is None:
        raise ConfigurationError(
            "two_phase_model: missing required configuration section (v0.4.4).\n"
            "  Add, for example:\n"
            '    "two_phase_model": {\n'
            '      "mode": "homogeneous_equilibrium",\n'
            '      "condensed_volume": "neglected",\n'
            '      "nozzle_entrainment": "complete"\n'
            "    }\n"
            "  ('homogeneous_equilibrium' is the standard model; "
            "'single_phase_legacy' reproduces the v0.4.3 behaviour for "
            "regression comparison only.)")
    _reject_unknown(tpm, _TWO_PHASE_KEYS, "two_phase_model")
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
    two_phase_mode = str(_choice(tpm, "mode", "two_phase_model",
                                 KNOWN_TWO_PHASE_MODES))
    two_phase_cv = str(_choice(tpm, "condensed_volume", "two_phase_model",
                               TWO_PHASE_CONDENSED_VOLUME_CHOICES))
    two_phase_ne = str(_choice(tpm, "nozzle_entrainment", "two_phase_model",
                               TWO_PHASE_NOZZLE_ENTRAINMENT_CHOICES))
    if backend == BACKEND_CEA_LEGACY_EXECUTABLE and \
            two_phase_mode == TWO_PHASE_HOMOGENEOUS:
        raise ConfigurationError(
            "two_phase_model.mode = 'homogeneous_equilibrium' cannot be used "
            "with backend 'cea_legacy_executable': the legacy fcea2 .plt output "
            "(p t gam m) does not provide the gas/condensed mass fractions, and "
            "GRIBS never estimates them.\n"
            "  Use backend 'cea_python' for the two-phase model, or set "
            "two_phase_model.mode = 'single_phase_legacy' (regression "
            "comparison only).")
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
    # molecular_weight: since v0.5.0-alpha the gas-phase EOS and the nozzle use
    # the GAS-PHASE molecular weight exclusively.  A leftover "total_MW" is
    # rejected with an explicit migration message - it is NEVER substituted
    # silently.
    py_mw_raw = cea_py_block.get("molecular_weight", "gas_phase_M")
    if py_mw_raw == "total_MW":
        raise ConfigurationError(
            "thermochemistry.cea_python.molecular_weight = 'total_MW' is no "
            "longer accepted (v0.5.0-alpha).\n"
            "  The two-phase chamber EOS p*Vg = mt*Yg(p)*Rg(p)*T0(p) and the "
            "nozzle model require the GAS-PHASE molecular weight; feeding the "
            "total molecular weight (which includes condensed species) into a "
            "gas-phase EOS is physically inconsistent.\n"
            "  Set molecular_weight = 'gas_phase_M'.  The total molecular "
            "weight remains available as the diagnostic column "
            "'total_molecular_weight_kg_kmol' in the property table and "
            "summary.json.")
    py_mw = _choice(cea_py_block, "molecular_weight", "thermochemistry.cea_python",
                    ("gas_phase_M",), required=False, default="gas_phase_M")
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
        transition_policy=str(_choice(
            nozzle,
            "transition_policy",
            "nozzle",
            ("jump", "shock"),
            required=False,
            default="jump",
        )),
        use_separation=bool(
            _boolean(sep, "enabled", "nozzle.flow_separation")
        ),
        sep_ratio=float(
            _number(
                sep,
                "pressure_ratio",
                "nozzle.flow_separation",
                exclusive_minimum=0.0,
            )
        ),
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
        two_phase_mode=two_phase_mode,
        two_phase_condensed_volume=two_phase_cv,
        two_phase_nozzle_entrainment=two_phase_ne,
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
    two_phase_active = (cfg.two_phase_mode == TWO_PHASE_HOMOGENEOUS
                        and bool(getattr(backend, "supports_two_phase", False)))
    print(f"  two-phase model       : {cfg.two_phase_mode}"
          + ("  (ACTIVE: total product mass conserved, EOS on the gas phase)"
             if two_phase_active else
             ("  (phase split inactive: Yg := 1)"
              if cfg.two_phase_mode == TWO_PHASE_HOMOGENEOUS else
              "  (v0.4.3 regression pathway)")))
    if two_phase_active:
        print(f"  nozzle entrainment    : {cfg.two_phase_nozzle_entrainment} "
              f"(no particle slip)")
        print(f"  condensed volume      : {cfg.two_phase_condensed_volume}")
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
            "gas_phase_molecular_weight_kg_kmol":
                list(backend.table.get("gas_phase_molecular_weight_kg_kmol", [])),
            "gas_constant_J_kgK": list(backend.table.get("gas_constant_J_kgK", [])),
            "gas_mass_fraction": list(backend.table.get("gas_mass_fraction", [])),
            "condensed_mass_fraction":
                list(backend.table.get("condensed_mass_fraction", [])),
            "psi_J_kg": list(backend.table.get("psi_J_kg", [])),
        }
        args.dump_thermo.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        print(f"Chamber-property table written to {args.dump_thermo}")
        return 0

    if args.selftest:
        return 0 if self_tests(cfg) else 1

    print("Running the internal-ballistics simulation ...")
    try:
        res = run_model(cfg)
        hist = sample(res, cfg)
    except (GribsError, ValueError, RuntimeError) as exc:
        print(f"Simulation failed: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 4
    if not self_tests_post_run(res, hist, cfg):
        print("Post-run conservation tests failed; aborting before writing output.",
              file=sys.stderr)
        return 1
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
        "two_phase_model": summ["two_phase_model"],
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
