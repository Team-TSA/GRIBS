#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GRIBS v0.6.2-alpha
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

  Nozzle (isentropic reference flow on the corrected mass basis, discharge
  coefficient Cd, thrust efficiency eta_F; two-phase split from the chamber
  composition)
      reference flow (R = R_eff = Ru/M_CEA = Yg*R_gas):
        choked      :  mdot_total_ref = Cd*At*p0*sqrt(g/(R*T0))
                                   *(2/(g+1))^((g+1)/(2(g-1)))
        unchoked    :  pe = pa, subsonic exit Mach from the pressure ratio
                       (this branch is *continuous* with the choked branch)
        separated   :  Summerfield criterion pe < f_sep*pa (only for Ae/At > 1)
      homogeneous-equilibrium mode (standard, complete entrainment):
        mdot_total     = mdot_total_ref   (R_eff already carries Yg)
        mdot_gas       = Yg_chamber * mdot_total
        mdot_condensed = Yc_chamber * mdot_total
      single-phase legacy mode:  mdot_total = mdot_gas = mdot_total_ref (Yg = 1)
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
  analytic derivative of Psi = R_eff*T0 (two-phase mode) and of Theta = R*T0
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
  State y = [p0, x, Rt, m_out_total, Impulse, m_gen_total,
  m_igniter_total] integrated with
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
      python GRIBS_v0.6.2-alpha.py --config another_config.json
  Run numerical self-tests only with:
      python GRIBS_v0.6.2-alpha.py --selftest

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
import json
import math
import platform
import sys
import tempfile
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# A-6 (package split, stage 1; ADR_A6_PACKAGE_SPLIT): the package lives in src/gribs next to this entry file.
# Names moved into the package are imported back below so that existing callers keep working. This entry
# is never a patch target for a moved name (tests/test_a6_patch_targets.py).
_GRIBS_PACKAGE_ROOT = Path(__file__).resolve().parent / "src"
if _GRIBS_PACKAGE_ROOT.is_dir() and str(_GRIBS_PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(_GRIBS_PACKAGE_ROOT))
# --- re-exports of moved names (A-6) ---
from gribs.constants import (  # noqa: E402, F401
    PROGRAM_NAME,
    PROGRAM_VERSION,
    SCHEMA_VERSION,
    LEGACY_SCHEMA_VERSION,
    PRE_TWO_PHASE_SCHEMA_VERSION,
    PRE_V050_SCHEMA_VERSION,
    TWO_PHASE_HOMOGENEOUS,
    TWO_PHASE_SINGLE_LEGACY,
    KNOWN_TWO_PHASE_MODES,
    TWO_PHASE_CONDENSED_VOLUME_CHOICES,
    TWO_PHASE_NOZZLE_ENTRAINMENT_CHOICES,
    BACKEND_CEA_PYTHON,
    BACKEND_CEA_LEGACY_EXECUTABLE,
    BACKEND_REFERENCE_CONSTANT,
    KNOWN_BACKENDS,
    LEGACY_BACKEND_ALIASES,
    REMOVED_BACKENDS,
    CEA_PROJECT_URL,
    CEA_INSTALL_COMMAND,
    CEA_MIN_API_VERSION,
    CEA_TESTED_API_VERSION,
    R_UNIVERSAL,
    _PHASE_FRACTION_ROUNDING_TOL,
    _MASS_BASIS_INVARIANT_REL_TOL,
    _MASS_BASIS_PHASE_SUM_ATOL,
    MODEL_STATUS_VALIDATED_YG_THRESHOLD,
    MODEL_STATUS_OUT_OF_SCOPE_YG_THRESHOLD,
    HIGH_CONDENSATE_G4B_UNVALIDATED_YG_THRESHOLD,
    G0,
    _TRAPZ,
    _PHASE_TAGS,
    POST_RUN_MASS_TOLERANCE,
)
from gribs.errors import (  # noqa: E402, F401
    GribsError,
    ConfigurationError,
    ThermochemistryError,
    ThermochemistryCacheError,
    ThermochemistryInvariantError,
    GammaBelowUnityError,
    BackendUnavailableError,
)
from gribs.util import (  # noqa: E402, F401
    _report_value_text,
)
from gribs.config import (  # noqa: E402, F401
    Config,
    PARAM_DOC,
    web_thickness,
    initial_propellant_volume,
    propellant_mass,
    geometry,
    _is_nonfinite_number,
    require_finite_config,
    validate,
)
from gribs.thermo import (  # noqa: E402, F401
    _PROGRAM_DIR,
    MISSING_CEA_PACKAGE_MESSAGE,
    _TABLE_REQUIRED_KEYS,
    _TABLE_TWO_PHASE_KEYS,
    ReactantSpec,
    ideal_mixture_density,
    _finite_positive,
    ThermochemistryBackend,
    PressureTableBackend,
    CEALegacyExecutableBackend,
    import_official_cea_module,
    _cea_package_metadata,
    _cea_database_identity,
    CEAPythonBackend,
    ReferenceConstantBackend,
    _thermo_state,
    _set_thermo,
    _ThermoSwap,
    build_backend,
    _CacheWriteError,
    _CACHE_CLEANUP_TEXT,
    _cache_write_failed,
    _cache_rename_noreplace,
    _cache_write_checkpoint,
    _cache_atomic_write_text,
    initialize_thermochemistry,
    thermochemistry,
    gas_props,
    psi_and_deriv,
    phase_fractions,
    theta_and_deriv,
    is_extrapolated,
    CEA_DENSITY_NOTE,
    N_GAS_PER_KG_TOTAL_DEFINITION,
    _sonic_column_status,
)
from gribs.propellant.burning import (  # noqa: E402, F401
    base_burn_rate,
    burn_rate,
)
from gribs.nozzle import (  # noqa: E402, F401
    NozzleTwoPhaseModel,
    area_mach,
    mach_from_area,
    critical_ratio,
    choke_threshold_pressure,
    choke_margin,
    choke_limit_pressure,
    nozzle_state,
)
from gribs.chamber.equilibrium import (  # noqa: E402, F401
    equilibrium_pressure,
)
from gribs.simulation.integrator import (  # noqa: E402, F401
    IP,
    IX,
    IRT,
    IMO,
    IIMP,
    IMG,
    IMIGN,
    igniter_mdot,
    igniter_injected_mass,
    make_rhs,
    _event,
    INTEGRATION_MAX_RHS_CALLS,
    INTEGRATION_MAX_WALL_S,
    INTEGRATION_MAX_TABLE_RESTARTS,
    _INTEGRATION_RESTART_FIRST_STEP_S,
    IntegrationBudgetError,
    _IntegrationBudget,
    _SegmentedSolution,
    _segmented_solve,
    run_model,
    require_implemented_transition_policy,
)
from gribs.simulation.postprocess import (  # noqa: E402, F401
    _phase_grid,
    sample,
    _masked_trapezoid,
    nozzle_transition_diagnostics,
    post_run_consistency_report,
    model_status_from_yg,
    _finite_number_or_none,
    _validation_range_fields,
    model_status_basis,
    MODEL_STATUS_PREVIOUS_TOKENS,
    SEMANTIC_STATUS_PREVIOUS_TOKENS,
    semantic_status_for_backend,
    out_of_validated_range_error,
    PHASE_FRACTION_VALIDATION_ERROR_CODE,
    _PHASE_FRACTION_KEYS,
    PhaseFractionValidationError,
    _two_phase_is_active,
    phase_fraction_validation_error,
    model_applicability_error,
    append_model_applicability_warnings,
    append_thermal_efficiency_warning,
    summarize,
    append_run_verdict_warnings,
    run_verdict_passed,
    EXIT_CODE_SUCCESS,
    EXIT_CODE_VERIFICATION_FAILED,
    EXIT_CODE_MASS_LEDGER_FAILED,
    run_verification_exit_code,
    MASS_LEDGER_TOTAL_RELATIVE_TOLERANCE,
    MASS_LEDGER_PROPELLANT_RELATIVE_TOLERANCE,
    MASS_LEDGER_PROPELLANT_TOLERANCE_FLOOR,
    MASS_LEDGER_PROPELLANT_TOLERANCE_RTL_FACTOR,
    mass_energy_ledger,
)
# owner modules of moved names that tests patch: the entry reads those names as <module>.<name> (T4)
import gribs.config  # noqa: E402
import gribs.thermo  # noqa: E402
import gribs.propellant.burning  # noqa: E402
import gribs.nozzle  # noqa: E402
import gribs.chamber.equilibrium  # noqa: E402
import gribs.simulation.integrator  # noqa: E402
import gribs.simulation.postprocess  # noqa: E402

# ==============================================================================
# 0. PROGRAM IDENTITY AND CONSTANTS
# ==============================================================================
# moved to src/gribs/constants.py (A-6 unit 1); re-exported above


# ==============================================================================
# 0b. EXCEPTION HIERARCHY
# ==============================================================================
# moved to src/gribs/errors.py (A-6 unit 1); re-exported above


# ==============================================================================
# 1. INTERNAL RUNTIME CONFIGURATION
# Values are populated from the complete JSON configuration. Field defaults are
# implementation fallbacks required by the dataclass constructor, not normal user input.
# ==============================================================================
# moved to src/gribs/config.py (A-6 unit 2); re-exported above


# ==============================================================================
# 2. DERIVED GEOMETRY HELPERS
# ==============================================================================
# moved to src/gribs/config.py (A-6 unit 2); re-exported above


# ==============================================================================
# 3. THERMOCHEMISTRY BACKENDS
# ------------------------------------------------------------------------------
#  The internal-ballistics solver only ever asks a backend for:
#      props(p)                  -> (Rg [J/(kg K)], T0 [K], gamma_s [-])
#                                   GAS-PHASE properties only
#      psi_and_derivative(p)     -> (Psi = R_eff*T0 [J/kg], dPsi/dp [J/(kg Pa)])
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

# moved to src/gribs/thermo.py (A-6 unit 3b); re-exported above. Readers of the names that tests patch call
# gribs.thermo.<name> (T4); the backend state is reached only through _thermo_state() / _set_thermo() (T10).


# ==============================================================================
# 4. BURN RATE (Saint-Robert + temperature sensitivity + erosive burning)
#    Unchanged from v0.3.1-alpha.
# ==============================================================================
# moved to src/gribs/propellant/burning.py (A-6 unit 4); re-exported above. The entry reads burn_rate as
# gribs.propellant.burning.burn_rate (T4).


# ==============================================================================
# 5. NOZZLE MODEL
#    Gas-phase isentropic reference flow (unchanged from v0.3.1-alpha) plus the
#    v0.4.4 two-phase mass split.  MASS BASIS: the flow returned to the ODE is
#    mdot_total (gas + condensed), consistent with the pressure equation that
#    conserves the TOTAL product mass.
# ==============================================================================
# moved to src/gribs/nozzle.py (A-6 unit 4); re-exported above. The entry reads nozzle_state and
# choke_limit_pressure as gribs.nozzle.<name> (T4).


# ==============================================================================
# 6. QUASI-STEADY EQUILIBRIUM PRESSURE (diagnostic only)
# ==============================================================================
# moved to src/gribs/chamber/equilibrium.py (A-6 unit 4); re-exported above. The entry reads it as
# gribs.chamber.equilibrium.equilibrium_pressure (T4).


# ==============================================================================
# 7. ODE SYSTEM
#    y = [p0, x, Rt, m_out_total, Impulse, m_gen_total,
#         m_igniter_total]
#    The mass quadratures are TOTAL product masses (gas + condensed) in the
#    homogeneous-equilibrium mode.  Pressure equation (two-phase model B):
#      dp0/dt = [ Psi*(mdot_gen_total - mdot_out_total) - p0*Ab*r ]
#               / [ Vg*(1 - p0*Psi'/Psi) ]
#    with Psi = R_eff*T0 and the ANALYTIC dPsi/dp.  In the single-phase legacy
#    mode Psi := Theta = Rg*T0 (exact v0.4.3-alpha behaviour).
# ==============================================================================
# moved to src/gribs/simulation/integrator.py (A-6 unit 5); re-exported above. The entry reads run_model as
# gribs.simulation.integrator.run_model (T4); the state indices and igniter_mdot are bare re-exported reads.


# _is_nonfinite_number, require_finite_config, validate: moved to src/gribs/config.py (A-6 unit 2);
# re-exported above. Readers of validate call gribs.config.validate (tests patch it there).


# ==============================================================================
# 8. POST-PROCESSING
#    v0.5.0-alpha: history columns use the TOTAL product mass basis
#    (m_total_eos, m_total_balance, mdot_out_total/gas/condensed) plus the
#    equilibrium phase fractions and Psi.
# ==============================================================================
# moved to src/gribs/simulation/postprocess.py (A-6 unit 6); re-exported above. The entry reads _phase_grid,
# sample, summarize, mass_energy_ledger, model_applicability_error and nozzle_transition_diagnostics as
# gribs.simulation.postprocess.<name> (T4); the other names are bare re-exported reads.


# CEA_DENSITY_NOTE, N_GAS_PER_KG_TOTAL_DEFINITION, _sonic_column_status: moved to src/gribs/thermo.py
# (A-6 unit 3b, report 113 placement); re-exported above.


# _report_value_text: moved to src/gribs/util.py (A-6 unit 1); re-exported above


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
            "equilibrium_molecular_weight_kg_kmol": [M] * n,
            "effective_gas_constant_J_kgK": [R] * n,
            "gas_mass_fraction": [Yg] * n,
            "condensed_mass_fraction": [self._ycond] * n,
            "psi_J_kg": [R * T0 * c.eta_T0] * n,
        }

    def metadata(self) -> Dict[str, Any]:
        return {"backend": self.name, "synthetic_selftest_backend": True}


#: Tolerance of the single-phase-limit regression test.  With Yg = 1 the two
#: code paths evaluate bit-identical arithmetic, so any difference beyond
#: double-rounding level indicates a structural inconsistency.
SINGLE_PHASE_LIMIT_TOLERANCE = 1.0e-12


SELF_TEST_COMMON_NAMES = (
    "analytic vs numerical dTheta/dp",
    "geometry consistency dVg/dx = Ab",
    "area-Mach inversion round trip",
    "closed-form choked thrust identity",
    "finite initial derivatives",
    "chamber-property table physical",
    "phase fractions sum to 1 (table)",
    "phase fractions physical (table)",
    "Psi = R_eff*T0 (table + interpolator)",
    "analytic vs numerical dPsi/dp",
    "pressure factor 1 - p*Psi'/Psi > 0",
    "single-phase limit: model B == Theta model",
    "nozzle decomposition mdot_total = gas + condensed",
    "no-condensed compatibility (hom. vs legacy run)",
)


def expected_self_test_names(c: Config) -> Tuple[str, ...]:
    """Return the exact ordered self-test names for a configuration."""
    if math.isclose(
        c.eps_nozzle,
        1.0,
        rel_tol=0.0,
        abs_tol=1.0e-12,
    ):
        nozzle_name = "converging-nozzle continuity (mdot, F)"
    else:
        nozzle_name = "C-D nozzle transition mass-flow continuity"

    return (
        SELF_TEST_COMMON_NAMES[0],
        SELF_TEST_COMMON_NAMES[1],
        SELF_TEST_COMMON_NAMES[2],
        nozzle_name,
        *SELF_TEST_COMMON_NAMES[3:],
    )


def self_test_name_contract(
    c: Config,
    actual_names: Sequence[str],
) -> dict:
    """Compare reported self-test names with the exact expected contract."""
    expected = expected_self_test_names(c)
    actual = tuple(actual_names)

    missing = [
        name
        for name in expected
        if name not in actual
    ]
    unexpected = [
        name
        for name in actual
        if name not in expected
    ]
    duplicates = sorted({
        name
        for name in actual
        if actual.count(name) > 1
    })

    return {
        "expected_names": list(expected),
        "actual_names": list(actual),
        "expected_count": len(expected),
        "actual_count": len(actual),
        "missing_names": missing,
        "unexpected_names": unexpected,
        "duplicate_names": duplicates,
        "passed": bool(
            actual == expected
            and not missing
            and not unexpected
            and not duplicates
        ),
    }


def self_tests(c: Config, verbose=True) -> bool:
    """Numerical self-tests executed explicitly with ``--selftest``.

    All v0.4.3-alpha tests are kept unchanged (T1-T7); the two-phase model
    adds T8-T15 (see the test names).  Post-simulation conservation tests run
    in :func:`self_tests_post_run` after the ODE integration.
    """
    out, ok_all = [], True
    result_names = []

    def chk(name, ok, detail=""):
        nonlocal ok_all
        ok_all &= bool(ok)
        result_names.append(name)
        out.append(f"  [{'PASS' if ok else 'FAIL'}] {name:<46s} {detail}")

    backend = gribs.thermo.thermochemistry()
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
    pstar = gribs.nozzle.choke_limit_pressure(cc, At0, Ae)

    if pstar is None:
        chk("nozzle choking-boundary continuity", False,
            "choking-limit pressure could not be determined")
    else:
        lo = gribs.nozzle.nozzle_state(pstar * (1.0 - 1e-9), At0, Ae, cc)
        hi = gribs.nozzle.nozzle_state(pstar * (1.0 + 1e-9), At0, Ae, cc)
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
        # Mass basis: mdot_total comes from the isentropic formula evaluated
        # with R_eff and ve = sqrt(g*R_eff*Te), so the momentum term carries
        # the phase split through R_eff alone.  The pressure thrust does not:
        #   F = mdot_total*ve + (pe - pa)*At
        #     = At*(crit*p*(g + 1) - pa),  crit = (2/(g+1))^(g/(g-1))
        # This is the classical single-phase identity; no /Yg factor appears
        # because R_eff already carries Yg.
        cc5 = Config(**{**asdict(c), "eps_nozzle": 1.0, "Cd": 1.0,
                        "eta_thrust": 1.0, "use_separation": False})
        At0 = math.pi * cc5.R_t0 ** 2
        e = 0.0
        for p in np.geomspace(cc5.p_a / 0.5, 50e6, 20):
            nz = gribs.nozzle.nozzle_state(p, At0, At0, cc5)
            _, _, g = gribs.thermo.gas_props(p, cc5)
            F_ref = At0 * (critical_ratio(g) * p * (g + 1.0) - cc5.p_a)
            e = max(e, abs(nz["F"] / F_ref - 1.0))
        chk("closed-form choked thrust identity", e < 1e-12,
            f"max rel. err = {e:.2e}")

        # T6 finite and physical right-hand side at t = 0
        try:
            dy = make_rhs(c, c.eps_nozzle * math.pi * c.R_t0 ** 2, True)(
                0.0, np.array([
                    c.p0_init,
                    0.0,
                    c.R_t0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ]))
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

    # T10 Psi consistency: Psi = R_eff*T_used at every table point.
    # R_eff already carries Yg, so the solver-facing Psi and the stored
    # table column are the same quantity in both the two-phase mode and
    # the single-phase legacy mode (where Yg := 1 by definition).
    R_tab = np.asarray(table["effective_gas_constant_J_kgK"], dtype=float)
    T_tab = np.asarray(table["temperature_K"], dtype=float) * c.eta_T0
    psi_ref = R_tab * T_tab
    psi_solver_ref = psi_ref
    if "psi_J_kg" in table:
        psi_tab = np.asarray(table["psi_J_kg"], dtype=float)
    else:
        psi_tab = psi_ref
    e = float(np.max(np.abs(psi_tab - psi_ref) / np.maximum(psi_ref, 1e-300)))
    # also check the interpolated Psi against the mode-consistent product
    e_node = 0.0
    for p, ref in zip(table["pressure_Pa"], psi_solver_ref):
        e_node = max(e_node, abs(gribs.thermo.psi_and_deriv(float(p), c)[0] / ref - 1.0))
    chk("Psi = R_eff*T0 (table + interpolator)",
        e <= 1.0e-12 and e_node <= 1.0e-10,
        f"stored-vs-product rel. err = {e:.2e}, node interp err = {e_node:.2e}")

    # T11 analytic dPsi/dp against a centred difference
    e = 0.0
    for p in np.geomspace(max(c.p_fit_min * 1.05, 1.1e5), c.p_fit_max * 0.95, 25):
        hstep = p * 1e-5
        fd = (gribs.thermo.psi_and_deriv(p + hstep, c)[0] - gribs.thermo.psi_and_deriv(p - hstep, c)[0]) / (2 * hstep)
        dp = gribs.thermo.psi_and_deriv(p, c)[1]
        if dp == 0.0:
            e = max(e, abs(fd) / max(gribs.thermo.psi_and_deriv(p, c)[0], 1e-300) * p)
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
    # Cleanup P4 v2 (A-1): every comparison builds its backends with an
    # explicit two-phase mode, independent of the mode of the configuration
    # the self-tests were started with.
    def _cfg_mode(mode):
        cc = copy.copy(c)                      # keeps reactants_spec_list
        cc.two_phase_mode = mode
        return cc
    synth = _SyntheticTestBackend(_cfg_mode(TWO_PHASE_HOMOGENEOUS),
                                  condensed_fraction=0.0)
    try:
        with _ThermoSwap(synth):
            e_psi, e_dpsi, e_rhs = 0.0, 0.0, 0.0
            c_hom = Config(**{**asdict(c), "two_phase_mode": TWO_PHASE_HOMOGENEOUS})
            c_leg = Config(**{**asdict(c), "two_phase_mode": TWO_PHASE_SINGLE_LEGACY})
            # Cleanup P4 (A-1): the legacy side is evaluated with a backend
            # built in legacy mode, so the backend's own legacy branch is
            # exercised (the CEA/legacy recursion of A-6 3c24fbc went
            # undetected because both sides used the homogeneous backend).
            synth_leg = _SyntheticTestBackend(_cfg_mode(TWO_PHASE_SINGLE_LEGACY),
                                              condensed_fraction=0.0)
            for p in np.geomspace(synth.pmin * 1.01, synth.pmax * 0.99, 30):
                psi_h, dpsi_h = gribs.thermo.psi_and_deriv(float(p), c_hom)
                with _ThermoSwap(synth_leg):
                    th_l, dth_l = theta_and_deriv(float(p), c_leg)
                e_psi = max(e_psi, abs(psi_h / th_l - 1.0))
                if dth_l != 0.0:
                    e_dpsi = max(e_dpsi, abs(dpsi_h / dth_l - 1.0))
                else:
                    e_dpsi = max(e_dpsi, abs(dpsi_h - dth_l))
                y0 = np.array([
                    float(p),
                    0.3 * web_thickness(c_hom),
                    c.R_t0,
                    1e-5,
                    1e-4,
                    2e-5,
                    0.0,
                ])
                dy_h = make_rhs(c_hom, Ae, True)(0.0, y0)
                with _ThermoSwap(synth_leg):
                    dy_l = make_rhs(c_leg, Ae, True)(0.0, y0)
                e_rhs = max(e_rhs, float(np.max(np.abs(dy_h - dy_l)
                                                / np.maximum(np.abs(dy_l), 1e-300))))
            chk("single-phase limit: model B == Theta model",
                max(e_psi, e_dpsi, e_rhs) < SINGLE_PHASE_LIMIT_TOLERANCE,
                f"Psi err = {e_psi:.2e}, dPsi err = {e_dpsi:.2e}, "
                f"rhs err = {e_rhs:.2e} (tol {SINGLE_PHASE_LIMIT_TOLERANCE:.0e})")

            # T14 nozzle flow decomposition: mdot_total = mdot_gas + mdot_condensed,
            # mdot_gas = Yg*mdot_total, mdot_condensed = Yc*mdot_total.
            synth_c = _SyntheticTestBackend(_cfg_mode(TWO_PHASE_HOMOGENEOUS),
                                            condensed_fraction=0.3)
            with _ThermoSwap(synth_c):
                e_dec = 0.0
                c_hom_c = Config(**{**asdict(c), "two_phase_mode": TWO_PHASE_HOMOGENEOUS})
                for p in np.geomspace(c.p_a * 1.05, c.p_fit_max, 30):
                    nz = gribs.nozzle.nozzle_state(float(p), At0, Ae, c_hom_c)
                    if nz["mdot_total"] <= 0.0:
                        continue
                    Yg_p, Yc_p = gribs.thermo.phase_fractions(float(p), c_hom_c)
                    e_dec = max(e_dec,
                                abs(nz["mdot_gas"] + nz["mdot_condensed"]
                                    - nz["mdot_total"]) / nz["mdot_total"],
                                abs(nz["mdot_gas"] / nz["mdot_total"] - Yg_p),
                                abs(nz["mdot_condensed"] / nz["mdot_total"] - Yc_p),
                                abs(nz["mdot_total"] - nz["mdot_total_reference"])
                                / nz["mdot_total_reference"])
                chk("nozzle decomposition mdot_total = gas + condensed",
                    e_dec < 1.0e-12, f"max rel. err = {e_dec:.2e}")

            # T15 compatibility regression: with a condensed-FREE table the
            # homogeneous-equilibrium run must match the legacy run (short motor).
            fast = {**asdict(c), "a_burn": 10.0 * c.a_burn, "blowdown": False}
            c_hom_f = Config(**{**fast, "two_phase_mode": TWO_PHASE_HOMOGENEOUS})
            c_leg_f = Config(**{**fast, "two_phase_mode": TWO_PHASE_SINGLE_LEGACY})
            res_h = gribs.simulation.integrator.run_model(c_hom_f)
            hist_h = gribs.simulation.postprocess.sample(res_h, c_hom_f)
            s_h = gribs.simulation.postprocess.summarize(res_h, hist_h, c_hom_f)
            with _ThermoSwap(synth_leg):
                res_l = gribs.simulation.integrator.run_model(c_leg_f)
                hist_l = gribs.simulation.postprocess.sample(res_l, c_leg_f)
                s_l = gribs.simulation.postprocess.summarize(res_l, hist_l, c_leg_f)
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

    name_contract = self_test_name_contract(c, result_names)
    ok_all &= bool(name_contract["passed"])

    if verbose:
        print(f"Self-tests  (thermochemistry backend: {backend.name}, "
              f"two-phase mode: {c.two_phase_mode}"
              f"{'' if two_phase_active else ' (phase split inactive)'})")
        print("\n".join(out))
        if not name_contract["passed"]:
            print("SELF-TEST CONTRACT ERROR")
            print(
                "  expected count: "
                f"{name_contract['expected_count']}, "
                "actual count: "
                f"{name_contract['actual_count']}"
            )
            if name_contract["missing_names"]:
                print(
                    "  missing names: "
                    + ", ".join(name_contract["missing_names"])
                )
            if name_contract["unexpected_names"]:
                print(
                    "  unexpected names: "
                    + ", ".join(name_contract["unexpected_names"])
                )
            if name_contract["duplicate_names"]:
                print(
                    "  duplicate names: "
                    + ", ".join(name_contract["duplicate_names"])
                )
        print(f"  -> {'all checks passed' if ok_all else 'FAILURES DETECTED'}\n")
    return ok_all


# POST_RUN_MASS_TOLERANCE: moved to src/gribs/constants.py (A-6 unit 1); re-exported above


def self_tests_post_run(
    res: dict,
    h: dict,
    c: Config,
    verbose=True,
    report: Optional[dict] = None,
) -> bool:
    """Display and return the structured post-run consistency verdict."""
    report = (
        post_run_consistency_report(res, h)
        if report is None
        else report
    )

    final = report["final_total_mass"]
    eos = report["eos_history"]

    out = [
        (
            "  [{status}] {name:<46s} "
            "mt0 + m_gen = {lhs:.9e} kg vs "
            "m_out + mt = {rhs:.9e} kg "
            "(rel. err = {error:.2e})"
        ).format(
            status="PASS" if final["passed"] else "FAIL",
            name="total mass conservation (final)",
            lhs=final["lhs_initial_plus_generated_kg"],
            rhs=final["rhs_discharged_plus_remaining_kg"],
            error=final["relative_error"],
        ),
        (
            "  [{status}] {name:<46s} "
            "max rel. err = {error:.2e} "
            "(fixed reference = {reference:.9e} kg, "
            "tol {tolerance:.0e})"
        ).format(
            status="PASS" if eos["passed"] else "FAIL",
            name="EOS mass vs balance mass (history)",
            error=eos["invariant_relative_error"],
            reference=eos["reference_mass_kg"],
            tolerance=report["tolerance"],
        ),
        (
            "  [INFO] {name:<46s} "
            "max rel. err = {error:.2e}"
        ).format(
            name="legacy shrinking-mass metric (deprecated)",
            error=eos[
                "legacy_instantaneous_mass_relative_error"
            ],
        ),
    ]

    if verbose:
        print("Post-run conservation tests")
        print("\n".join(out))
        print(
            f"  -> {'all checks passed' if report['passed'] else 'FAILURES DETECTED'}\n"
        )

    return bool(report["passed"])


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
    cea_status, _cea_missing, cea_invalid = cea_diagnostics_status(cea_diag)
    if cea_status == "ok":
        lines += [
            "",
            "CEA theoretical (comparison only)",
            f"  {'c*':<18}{cea_diag['cea_theoretical_cstar_m_s']:>11.1f} m/s",
            f"  {'CF':<18}{cea_diag['cea_theoretical_cf']:>11.4f} -",
            f"  {'Isp':<18}{cea_diag['cea_theoretical_isp_s']:>11.2f} s",
            f"  {'GRIBS c* (effective)':<18}{s['cstar_eff_m_s']:>11.1f} m/s",
        ]
    elif cea_status != "not_requested":
        # Show the diagnostic state; never draw a substituted CEA value.
        def _fig_value(key: str, fmt: str, unit: str) -> str:
            if not isinstance(cea_diag, dict) or key not in cea_diag:
                return f"{'n/a':>11}"
            if key in cea_invalid:
                return f"{'invalid':>11}"
            return f"{format(cea_diag[key], fmt):>11} {unit}"
        lines += [
            "",
            f"CEA theoretical (comparison only): {cea_status.upper()}",
        ]
        if cea_status != "failed":
            lines += [
                f"  {'c*':<18}{_fig_value('cea_theoretical_cstar_m_s', '.1f', 'm/s')}",
                f"  {'CF':<18}{_fig_value('cea_theoretical_cf', '.4f', '-')}",
                f"  {'Isp':<18}{_fig_value('cea_theoretical_isp_s', '.2f', 's')}",
            ]
        lines.append(f"  {'GRIBS c* (effective)':<18}{s['cstar_eff_m_s']:>11.1f} m/s")
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
    R0, T00, g00 = gribs.thermo.gas_props(c.p0_init, c)
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
               "mdot_out_condensed", "m_gen", "m_igniter", "m_out",
               "m_total_eos",
               "m_total_balance", "m_gas_equilibrium", "m_condensed",
               "Y_gas", "Y_condensed", "Psi",
               "regime", "Me", "pe", "ve", "F", "impulse",
               "cstar", "CF", "R", "T0", "gamma", "extrap"]
CSV_UNITS = {"t": "s", "p0": "Pa", "x": "m", "Ri": "m", "Lp": "m", "Ab": "m2",
             "Vg": "m3", "Kn": "-", "Rt": "m", "At": "m2", "r": "m/s",
             "mdot_gen_total": "kg/s", "mdot_out_total": "kg/s",
             "mdot_out_gas": "kg/s", "mdot_out_condensed": "kg/s",
             "m_gen": "kg", "m_igniter": "kg", "m_out": "kg",
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
    cea_status, cea_missing, cea_invalid = cea_diagnostics_status(cea_diag)
    if cea_status == "ok":
        A("  CEA theoretical comparison only (never used by the GRIBS model):")
        A(f"    cea_theoretical_cstar_m_s      : {cea_diag['cea_theoretical_cstar_m_s']:.3f}")
        A(f"    cea_theoretical_cf             : {cea_diag['cea_theoretical_cf']:.5f}")
        A(f"    cea_theoretical_isp_s          : {cea_diag['cea_theoretical_isp_s']:.4f}")
        A(f"    cea_theoretical_isp_vacuum_s   : {cea_diag['cea_theoretical_isp_vacuum_s']:.4f}")
        A(f"    evaluation                     : "
          f"{cea_diag['cea_theoretical_exit_definition']} at "
          f"p_c = {cea_diag['cea_theoretical_chamber_pressure_Pa']/1e6:.6f} MPa")
    elif cea_status != "not_requested":
        A("  CEA theoretical comparison only (never used by the GRIBS model):")
        A(f"    status                         : {cea_status} "
          "(no CEA comparison value is substituted)")
        if cea_status == "failed":
            A(f"    error                          : {cea_diag.get('cea_theoretical_error')}")
        else:
            if cea_missing:
                A(f"    missing fields                 : {', '.join(cea_missing)}")
            if cea_invalid:
                A(f"    invalid fields                 : {', '.join(cea_invalid)}")
            if isinstance(cea_diag, dict):
                for _key, _fmt in (("cea_theoretical_cstar_m_s", ".3f"),
                                   ("cea_theoretical_cf", ".5f"),
                                   ("cea_theoretical_isp_s", ".4f"),
                                   ("cea_theoretical_isp_vacuum_s", ".4f")):
                    A(f"    {_key:<31}: "
                      f"{_cea_diag_display(cea_diag, _key, _fmt, cea_invalid)}")
    if cea_status != "not_requested":
        A(f"    GRIBS effective c*             : {s['cstar_eff_m_s']:.3f} m/s")
        A(f"    GRIBS effective CF             : {s['CF_eff']:.5f}")
        A(f"    GRIBS effective Isp            : {s['Isp_s']:.4f} s")
    if cea_status == "ok":
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
    run_warnings = [
        warning
        for warning in s["warnings"]
        if warning.get("code") in {
            "W_POST_RUN_CONSISTENCY",
            "W_MASS_LEDGER_FAILED",
            "W_HIGH_CONDENSATE_MODEL_LIMITATION",
            "W_HIGH_CONDENSATE_G4B_UNVALIDATED",
            "W_TWO_PHASE_NOZZLE_APPROXIMATION",
            "W_THERMAL_EFFICIENCY_PARTIAL",
            "E_OUT_OF_VALIDATED_RANGE",
        }
    ]
    if run_warnings:
        A("[ verification warnings ]")
        for warning in run_warnings:
            A(
                f"  WARNING [{warning['code']}] "
                f"severity={warning['severity']}"
            )
            A(f"    {warning['message']}")
        A("")

    A("[ mass ledger ]")
    ledger = s["ledger"]
    mass = ledger["mass"]
    total = mass["total_system"]
    propellant = mass["propellant"]
    phase = mass["phase_discharge_diagnostic"]
    A(f"  basis                            : {mass['basis']}")
    A(f"  initial chamber inventory        : "
      f"{mass['initial_chamber_inventory_kg']:.9e} kg")
    A(f"  generated products, total        : "
      f"{mass['generated_mass_total_kg']:.9e} kg")
    A(f"    propellant-derived             : "
      f"{mass['generated_mass_propellant_kg']:.9e} kg")
    A(f"    igniter-derived, injected      : "
      f"{mass['igniter_mass_injected_kg']:.9e} kg")
    A(f"    igniter charge, commanded      : "
      f"{mass['igniter_mass_commanded_kg']:.9e} kg")
    A(f"  discharged products, ODE state   : "
      f"{mass['discharged_mass_total_kg']:.9e} kg")
    A(f"  final chamber inventory, EOS     : "
      f"{mass['final_chamber_inventory_kg']:.9e} kg")
    A(f"  initial solid propellant         : "
      f"{mass['initial_propellant_mass_kg']:.9e} kg")
    A(f"  final unburned propellant        : "
      f"{mass['final_unburned_propellant_kg']:.9e} kg")
    A(f"  burned propellant, geometry      : "
      f"{mass['burned_propellant_geometry_kg']:.9e} kg")
    A(f"  total-system signed residual     : "
      f"{total['signed_residual_kg']:+.9e} kg")
    A(f"  total-system relative residual   : "
      f"{total['relative_residual']:.3e} "
      f"(limit {total['tolerance']:.1e}, "
      f"{'PASS' if total['passed'] else 'FAIL'})")
    A(f"  propellant signed residual       : "
      f"{propellant['signed_residual_kg']:+.9e} kg")
    A(f"  propellant relative residual     : "
      f"{propellant['relative_residual']:.3e} "
      f"(limit {propellant['tolerance']:.1e}, "
      f"{'PASS' if propellant['passed'] else 'FAIL'})")
    A(f"  phase-discharge diagnostic       : {phase['method']} "
      f"(grid dependent)")
    A(f"    gas / condensed                : "
      f"{phase['gas_mass_kg']:.9e} / "
      f"{phase['condensed_mass_kg']:.9e} kg")
    A(f"    quadrature mismatch            : "
      f"{phase['relative_quadrature_mismatch']:.3e}")
    A(f"  mass-ledger result               : "
      f"{'PASS' if mass['passed'] else 'FAIL'}")
    A("")
    A("[ energy ledger ]")
    energy = ledger["energy"]
    A(f"  closure residual available       : "
      f"{energy['closure_residual_available']}")
    A(f"  reason                           : {energy['reason']}")
    A(f"  equilibrium constraint           : "
      f"{energy['model']['equilibrium_constraint']}")
    A(f"  temperature treatment            : "
      f"{energy['model']['temperature_treatment']}")
    A(f"  state relation                   : "
      f"{energy['model']['state_relation']}")
    A("  NOTE: the state relation is not a complete energy-conservation law.")
    A("")
    A("[ verification ]")
    A(f"  legacy propellant balance metric : "
      f"{s['propellant_mass_balance_error']:.3e} "
      f"(includes igniter in the legacy numerator)")
    A(f"  EOS-history invariant error      : "
      f"{s['total_mass_consistency_error_invariant']:.3e}")
    A(f"  invariant reference mass         : "
      f"{s['total_mass_consistency_reference_mass_kg']:.9e} kg")
    A(f"  invariant metric version         : "
      f"{s['total_mass_consistency_metric_version']}")
    A(f"  legacy shrinking-mass error      : "
      f"{s['total_mass_consistency_error']:.3e} (deprecated)")
    A(f"  property-range extrapolation used: {s['property_extrapolation']}"
      f"   (policy = {c.property_policy})")
    if s["property_extrapolation"] and c.property_policy == "extrapolate":
        A("  WARNING: results rely on extrapolated thermochemical properties.")
    # --- A-5.1-5: validation status and c* indicators ---
    val = s.get("validation", {})
    if val:
        A("")
        A("[ validation status ]")
        A(f"  semantic status                : {val.get('semantic_status', 'n/a')}")
        A(f"  numerical status               : {val.get('numerical_status', 'n/a')}")
        A(f"  model status                   : {val.get('model_status', 'n/a')}")
        A(f"  table pressure range [not validated] : {val.get('table_pressure_range_Pa', 'n/a')}")
        A(f"  table Yg range [not validated]       : {val.get('table_Yg_range', 'n/a')}")
        A(f"  mass basis semantics version   : {val.get('mass_basis_semantics_version', 'n/a')}")
    lmc = s.get("local_model_cstar_m_s")
    if lmc is not None and not (isinstance(lmc, float) and math.isnan(lmc)):
        A(f"  local model c* (burn median)   : {lmc:.2f} m/s")
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
            "migrated automatically: since v0.5.0-alpha the chamber equation "
            "of state and the nozzle model use M_CEA = EqSolution.M = "
            "1/EqSolution.n [kg-total/kmol-gas] only (EqSolution.MW also "
            "counts condensed moles), and a silent substitution is forbidden.\n"
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
                     "cea_legacy_executable", "reference_constant")
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
_REFERENCE_CONSTANT_KEYS = (
    "temperature_K",
    "gamma",
    "molecular_weight_kg_kmol",
)
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
    def _reject_nonstandard_constant(token):   # NaN / Infinity / -Infinity are not JSON
        raise ConfigurationError(
            f"{path}: non-standard JSON constant {token!r} (NaN / Infinity are not valid JSON, RFC 8259).")
    try:
        doc = json.loads(path.read_text(encoding="utf-8"), parse_constant=_reject_nonstandard_constant)
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


def _reject_nonfinite_values(node, path="root"):
    """Reject any non-finite number in a configuration document (all sections)."""
    if isinstance(node, dict):
        for k, v in node.items():
            _reject_nonfinite_values(v, f"{path}.{k}")
    elif isinstance(node, list):
        for i, v in enumerate(node):
            _reject_nonfinite_values(v, f"{path}[{i}]")
    elif _is_nonfinite_number(node):
        raise ConfigurationError(f"{path}: non-finite number {node!r} is not allowed in a configuration.")


def configuration_from_document(doc: dict) -> Configuration:
    """Validate a v0.5 document and build the runtime configuration."""
    _reject_nonfinite_values(doc)
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
    reference_constant = _section(
        th,
        "reference_constant",
        "thermochemistry",
        required=(backend == BACKEND_REFERENCE_CONSTANT),
    )
    if reference_constant is not None:
        _reject_unknown(
            reference_constant,
            _REFERENCE_CONSTANT_KEYS,
            "thermochemistry.reference_constant",
        )

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
    # molecular_weight: since v0.5.0-alpha the chamber EOS and the nozzle use
    # M_CEA = EqSolution.M = 1/n [kg-total/kmol-gas] exclusively (token name
    # 'gas_phase_M' is a frozen compatibility name).  A leftover "total_MW" is
    # rejected with an explicit migration message - it is NEVER substituted
    # silently.
    py_mw_raw = cea_py_block.get("molecular_weight", "gas_phase_M")
    if py_mw_raw == "total_MW":
        raise ConfigurationError(
            "thermochemistry.cea_python.molecular_weight = 'total_MW' is no "
            "longer accepted (v0.5.0-alpha).\n"
            "  The two-phase chamber EOS and the nozzle model use "
            "M_CEA = EqSolution.M = 1/EqSolution.n [kg-total/kmol-gas] "
            "(R_eff = Ru/M_CEA = Yg*R_gas); feeding EqSolution.MW (whose mole "
            "count includes condensed species) into the EOS is physically "
            "inconsistent.\n"
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

    reference_block = _mapping(
        reference_constant or {},
        "thermochemistry.reference_constant",
    )
    reference_temperature = float(_number(
        reference_block,
        "temperature_K",
        "thermochemistry.reference_constant",
        required=(backend == BACKEND_REFERENCE_CONSTANT),
        default=3000.0,
        exclusive_minimum=0.0,
    ))
    reference_gamma = float(_number(
        reference_block,
        "gamma",
        "thermochemistry.reference_constant",
        required=(backend == BACKEND_REFERENCE_CONSTANT),
        default=1.2,
        exclusive_minimum=1.0,
    ))
    reference_molecular_weight = float(_number(
        reference_block,
        "molecular_weight_kg_kmol",
        "thermochemistry.reference_constant",
        required=(backend == BACKEND_REFERENCE_CONSTANT),
        default=25.0,
        exclusive_minimum=0.0,
    ))

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
        reference_temperature_K=reference_temperature,
        reference_gamma=reference_gamma,
        reference_molecular_weight_kg_kmol=reference_molecular_weight,
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


def load_json_configuration(path: Path,
                            backend_override: Optional[str] = None) -> Configuration:
    """Load and validate a configuration file.

    With ``backend_override`` (the CLI ``--backend`` value) the configuration
    is validated twice:

    1. the loaded document exactly as before (after any schema migration), so
       an override never makes an invalid file acceptable - in particular the
       file must name a valid backend of its own (``_backend_choice``) and
       satisfy that backend's own rules;
    2. a copy with ``thermochemistry.backend`` replaced by the override, so
       every rule that depends on the backend (the legacy/homogeneous
       exclusion, the backend-conditional sections, the cache settings) is
       applied to the backend that will actually run.

    The recorded ``document`` / ``original_configuration`` is the loaded
    document (after migration, before the override).
    """
    doc, notes = _load_document(path)
    effective = doc
    if backend_override is not None:
        configuration_from_document(doc)          # step 1: the file on its own
        effective = copy.deepcopy(doc)
        thermo = effective["thermochemistry"]    # a mapping after step 1
        thermo["backend"] = backend_override
        if backend_override == BACKEND_REFERENCE_CONSTANT:
            # Phase A-5 contract (verification/a5-reference/A5_DESIGN.md):
            # reference values the file does not state come from the
            # documented Config defaults (absent section, or missing keys when
            # the file was written for another backend), exactly as before
            # the override was moved ahead of validation.  A file written for
            # reference_constant already passed step 1 with a complete section.
            fields = Config.__dataclass_fields__
            defaults = {
                "temperature_K": fields["reference_temperature_K"].default,
                "gamma": fields["reference_gamma"].default,
                "molecular_weight_kg_kmol":
                    fields["reference_molecular_weight_kg_kmol"].default,
            }
            if "reference_constant" not in thermo:
                thermo["reference_constant"] = dict(defaults)
            elif isinstance(thermo["reference_constant"], dict):
                for key, value in defaults.items():
                    thermo["reference_constant"].setdefault(key, value)
    loaded = configuration_from_document(effective)
    if backend_override is not None:
        loaded.document = doc
        loaded.config.original_configuration = copy.deepcopy(doc)  # type: ignore[attr-defined]
        notes = list(notes) + [
            f"thermochemistry.backend overridden on the command line to "
            f"{backend_override!r}."]
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
    p.add_argument(
        "--allow-verification-failure",
        action="store_true",
        help=(
            "return exit 0 only for a mass-ledger-only verification failure "
            "(diagnostics remain marked failed)"
        ),
    )
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
    elif backend.name == BACKEND_REFERENCE_CONSTANT:
        print("  property model        : constant analytic properties")
        print(
            "  reference properties : "
            f"T0 = {backend._temperature:g} K, "
            f"gamma = {backend._gamma:g}, "
            f"M = {backend._molecular_weight:g} kg/kmol"
        )
        print("  CEA used              : no")
        print("  property cache        : disabled")
        print(
            f"  validity range        : {backend.pmin/1e5:g} - "
            f"{backend.pmax/1e5:g} bar, "
            f"policy = {backend.outside_range_policy}"
        )
    else:
        print(f"  property table        : {backend.c.p_fit_min/1e5:g} - "
              f"{backend.c.p_fit_max/1e5:g} bar, "
              f"{backend.c.cea_pressure_points} points, "
              f"eta_T0 = {backend.c.eta_T0:g}, "
              f"policy = {backend.c.property_policy}")
    if getattr(backend, "cache_path", None):
        print(f"  cache file            : {backend.cache_path}")
        print(f"  cache reused          : {backend.cache_used}")
    print(f"  ideal mixture density : {backend.ideal_mixture_density:.6f} kg/m3")
    print(f"  packing fraction      : {backend.packing_fraction:.6f}")
    print(f"  bulk rho_p used       : {cfg.rho_p:.6f} kg/m3")


# Fields of the CEA theoretical diagnostics that the text summary and the
# figure display.  Numeric fields must be finite real numbers; the exit
# definition must be a string.
_CEA_DIAG_DISPLAY_NUMERIC_FIELDS = (
    "cea_theoretical_cstar_m_s",
    "cea_theoretical_cf",
    "cea_theoretical_isp_s",
    "cea_theoretical_isp_vacuum_s",
    "cea_theoretical_chamber_pressure_Pa",
)
_CEA_DIAG_DISPLAY_TEXT_FIELDS = ("cea_theoretical_exit_definition",)


def _cea_diag_value_ok(key: str, value: Any) -> bool:
    if key in _CEA_DIAG_DISPLAY_TEXT_FIELDS:
        return isinstance(value, str)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return False
    try:
        return math.isfinite(float(value))
    except (OverflowError, ValueError, TypeError):
        return False


def cea_diagnostics_status(cea_diag: Optional[dict]) -> Tuple[str, List[str], List[str]]:
    """Classify CEA theoretical diagnostics for the output writers.

    Returns ``(status, missing_fields, invalid_fields)`` with status
      'not_requested' : None (diagnostics disabled or not applicable),
      'failed'        : the diagnostic computation raised (error recorded),
      'ok'            : every displayed field is present and valid, and no
                        other cea_theoretical_* number is non-finite,
      'incomplete'    : displayed fields are missing (none invalid),
      'invalid'       : a field is non-finite / of the wrong type.
    Never used by the GRIBS model or the run verdict.
    """
    if cea_diag is None:
        return "not_requested", [], []
    if not isinstance(cea_diag, dict):
        return "invalid", [], ["<not a mapping>"]
    if "cea_theoretical_error" in cea_diag:
        return "failed", [], []
    display = _CEA_DIAG_DISPLAY_NUMERIC_FIELDS + _CEA_DIAG_DISPLAY_TEXT_FIELDS
    missing = [k for k in display if k not in cea_diag]
    invalid = [k for k in display if k in cea_diag and not _cea_diag_value_ok(k, cea_diag[k])]
    for k, v in cea_diag.items():
        if (k not in display and isinstance(v, float) and not math.isfinite(v)):
            invalid.append(k)
    if invalid:
        return "invalid", missing, invalid
    if missing:
        return "incomplete", missing, invalid
    return "ok", [], []


def _cea_diag_display(cea_diag: dict, key: str, fmt: str, invalid: Sequence[str]) -> str:
    """Formatted value, or an explicit marker; never a substituted number."""
    if key not in cea_diag:
        return "unavailable (missing)"
    if key in invalid:
        return f"invalid ({cea_diag[key]!r})"
    return format(cea_diag[key], fmt)


def _cea_theoretical_diagnostics(cfg: Config, backend: ThermochemistryBackend,
                                 summ: dict) -> Optional[dict]:
    """Theoretical CEA rocket comparison values (never fed back into GRIBS).

    A failed, incomplete or invalid diagnostic is returned with an explicit
    ``cea_theoretical_status`` (plus the error or the missing / invalid field
    names) so that it is preserved in summary.json and shown in the text
    summary and the figure; a successful result is returned unchanged.
    """
    if not getattr(cfg, "cea_py_rocket_diagnostics", False):
        return None
    if not hasattr(backend, "theoretical_rocket_performance"):
        return None
    p_ref = summ.get("p_mean_burn_Pa")
    if not (isinstance(p_ref, float) and math.isfinite(p_ref) and p_ref > 0.0):
        p_ref = summ.get("p_max_Pa")
    try:
        result = backend.theoretical_rocket_performance(float(p_ref), float(cfg.eps_nozzle))
    except Exception as exc:                                   # diagnostic only
        print(f"WARNING: CEA theoretical rocket diagnostics unavailable: "
              f"{type(exc).__name__}: {exc}")
        return {"cea_theoretical_status": "failed",
                "cea_theoretical_error": f"{type(exc).__name__}: {exc}"}
    status, missing, invalid = cea_diagnostics_status(result)
    if status == "ok":
        return result
    annotated = dict(result) if isinstance(result, dict) else {"cea_theoretical_raw": repr(result)[:200]}
    annotated["cea_theoretical_status"] = status
    if missing:
        annotated["cea_theoretical_missing_fields"] = list(missing)
    if invalid:
        annotated["cea_theoretical_invalid_fields"] = list(invalid)
    print(f"WARNING: CEA theoretical rocket diagnostics {status}: "
          f"missing={list(missing)} invalid={list(invalid)}")
    return annotated


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        if args.backend:
            loaded = load_json_configuration(args.config,
                                             backend_override=args.backend)
        else:
            loaded = load_json_configuration(args.config)
    except ConfigurationError as exc:
        print(f"Configuration error:\n{exc}", file=sys.stderr)
        return 2

    cfg = loaded.config
    if args.no_cache:
        cfg.cea_rebuild_cache = True
        loaded.migration_notes.append("Chamber-property cache rebuild forced by --no-cache.")

    for note in loaded.migration_notes:
        print(f"WARNING: {note}", file=sys.stderr)

    try:
        gribs.config.validate(cfg)
    except (ConfigurationError, ValueError) as exc:
        print(f"Configuration error:\n{exc}", file=sys.stderr)
        return 2

    if args.validate_config:
        print(f"Configuration {loaded.config.configuration_source} is valid "
              f"(schema {SCHEMA_VERSION}, backend {cfg.thermo_backend}).")
        return 0

    cache_dir = loaded.output_dir / "thermo_cache"
    try:
        backend = gribs.thermo.initialize_thermochemistry(cfg, cache_dir)
    except (ThermochemistryError, ConfigurationError) as exc:
        error_code = getattr(exc, "code", None)
        if error_code:
            print(
                f"Thermochemistry backend error [{error_code}]:\n{exc}",
                file=sys.stderr,
            )
            error_record = {
                "code": error_code,
                "severity": "error",
                "message": str(exc),
                "context": getattr(exc, "context", {}),
            }
            print(
                json.dumps(error_record, ensure_ascii=False, sort_keys=True),
                file=sys.stderr,
            )
        else:
            print(f"Thermochemistry backend error:\n{exc}", file=sys.stderr)
        return 3

    _print_backend_banner(cfg, backend)

    if args.dump_thermo is not None:
        payload = {
            "backend": backend.name,
            "provenance": backend.metadata(),
            **backend.thermo_dump_table(),
        }
        args.dump_thermo.write_text(
            json.dumps(payload, indent=2) + "\n",
            encoding="utf-8",
        )
        print(f"Chamber-property table written to {args.dump_thermo}")
        return 0

    if args.selftest:
        return 0 if self_tests(cfg) else 1

    applicability_error = gribs.simulation.postprocess.model_applicability_error(cfg, backend)
    if applicability_error is not None:
        print(
            f"Model applicability error [{applicability_error['code']}]: "
            f"{applicability_error['message']}",
            file=sys.stderr,
        )
        error_record = dict(applicability_error)
        error_record["context"] = {
            **applicability_error["context"],
            "simulation_started": False,
        }
        print(
            json.dumps(error_record, ensure_ascii=False, sort_keys=True),
            file=sys.stderr,
        )
        return 1

    print("Running the internal-ballistics simulation ...")
    try:
        res = gribs.simulation.integrator.run_model(cfg)
        hist = gribs.simulation.postprocess.sample(res, cfg)
    except (GribsError, ValueError, RuntimeError) as exc:
        print(f"Simulation failed: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 4
    _ib = (res.get("integration_budget") if isinstance(res, dict) else None) or {}
    if _ib.get("table_restarts_burn") or _ib.get("table_restarts_blowdown"):
        # Diagnostic record only (stderr, not JSON-first so that structured error parsers ignore it); the
        # summary / CSV schema is unchanged.
        print("INTEGRATION_TABLE_RESTARTS " + json.dumps(_ib, sort_keys=True), file=sys.stderr)
    summ = gribs.simulation.postprocess.summarize(res, hist, cfg, backend=backend)
    # If migration from v0.4.3 changed the physical model, record it as a
    # machine-readable warning (previously this only appeared in stderr and
    # configuration_migration.notes, making it invisible to automated checks).
    if loaded.migration_notes:
        _model_change_keywords = ("THIS CHANGES THE PHYSICAL MODEL",)
        for _note in loaded.migration_notes:
            if any(kw in _note for kw in _model_change_keywords):
                summ.setdefault("warnings", []).append({
                    "code": "W_MIGRATION_MODEL_CHANGE",
                    "severity": "warning",
                    "message": (
                        "Automatic configuration migration changed the physical "
                        "model.  See configuration_migration.notes in the summary "
                        "for details."
                    ),
                    "context": {"migration_note": _note},
                })
                break
    _cache_dur = getattr(backend, "cache_write_durability", None)      # cache durability warning
    if _cache_dur:
        summ.setdefault("warnings", []).append({
            "code": "W_THERMO_CACHE_DURABILITY_NOT_GUARANTEED", "severity": "warning",
            "message": ("The chamber-property cache was written, but fsync of its directory is not supported on "
                        "this file system; the cache file's durability after a crash or power loss is not "
                        "guaranteed."),
            "context": dict(_cache_dur)})
    append_run_verdict_warnings(summ)
    run_passed = run_verdict_passed(summ)
    verification_exit_code = run_verification_exit_code(
        summ,
        allow_verification_failure=args.allow_verification_failure,
    )

    thermo_md = backend.metadata()
    cea_diag = _cea_theoretical_diagnostics(cfg, backend, summ)

    outdir = loaded.output_dir
    outdir.mkdir(parents=True, exist_ok=True)
    png = outdir / loaded.output_names["figure_filename"]
    csvp = outdir / loaded.output_names["history_filename"]
    txtp = outdir / loaded.output_names["summary_text_filename"]
    jsnp = outdir / loaded.output_names["summary_json_filename"]

    text = summary_text(
        summ,
        cfg,
        thermo_md,
        cea_diag,
        loaded.migration_notes,
    )

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

    # Preserve machine-readable diagnostics before any verdict-based exit
    # and before figure generation.
    write_csv(hist, csvp)
    txtp.write_text(text, encoding="utf-8")
    jsnp.write_text(
        json.dumps(summary_json, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    make_figure(res, hist, summ, cfg, png, cea_diag)

    self_tests_post_run(
        res,
        hist,
        cfg,
        report=summ["post_run_consistency"],
    )

    print()
    print(text)
    print("Output files")
    for pth in (png, csvp, txtp, jsnp):
        print(f"  {pth}")

    if not run_passed:
        failed_codes = [
            warning["code"]
            for warning in summ["warnings"]
            if warning.get("code") in {
                "W_POST_RUN_CONSISTENCY",
                "W_MASS_LEDGER_FAILED",
                "E_OUT_OF_VALIDATED_RANGE",
            }
        ]
        if verification_exit_code == EXIT_CODE_SUCCESS:
            print(
                "Run verification failure allowed by "
                "--allow-verification-failure after diagnostic output was "
                "written: " + ", ".join(failed_codes),
                file=sys.stderr,
            )
            return EXIT_CODE_SUCCESS
        print(
            "Run verification failed after diagnostic output was written: "
            + ", ".join(failed_codes),
            file=sys.stderr,
        )
        return verification_exit_code

    return EXIT_CODE_SUCCESS


if __name__ == "__main__":
    sys.exit(main())
