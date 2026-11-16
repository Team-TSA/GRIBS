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

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

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
    CEADatabaseError,
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
    E_IONS_NOT_SUPPORTED,
    reject_unsupported_cea_ions,
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
    W_EXTERNAL_CEA_DATABASE_IGNORED,
    E_CEA_DATABASE_ALREADY_INITIALIZED,
    import_cea_without_auto_init,
    _sha256_file,
    _external_cea_database_files,
    pin_bundled_cea_database,
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
    exit_area_throat_radius,
    throat_exceeds_exit_area,
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
    append_external_cea_database_warning,
    summarize,
    W_THROAT_EXCEEDS_EXIT_AREA,
    THROAT_EXCEEDS_EXIT_EVENT,
    throat_exceeds_exit_area_time,
    append_throat_exceeds_exit_area_warning,
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
from gribs.plotting import (  # noqa: E402, F401
    C_PRESS,
    C_THRUST,
    C_GEN,
    C_OUT,
    C_ACC,
    C_GREY,
    REGIME_FACE,
    REGIME_NAME,
    _style,
    _regime_spans,
    _shade,
    KEY_RESULTS_MAX_FONTSIZE,
    KEY_RESULTS_MIN_FONTSIZE,
    KEY_RESULTS_LINESPACING,
    KEY_RESULTS_BOX_PAD,
    KEY_RESULTS_CHAR_WIDTH,
    key_results_fontsize,
    _draw_key_results,
    make_figure,
)
from gribs.reporting import (  # noqa: E402, F401
    CSV_COLUMNS,
    CSV_UNITS,
    write_csv,
    nozzle_transition_summary_lines,
    summary_text,
    _CEA_DIAG_DISPLAY_NUMERIC_FIELDS,
    _CEA_DIAG_DISPLAY_TEXT_FIELDS,
    _cea_diag_value_ok,
    cea_diagnostics_status,
    _cea_diag_display,
)
from gribs.selftest import (  # noqa: E402, F401
    DERIVATIVE_TOLERANCE_TABLE,
    _SyntheticTestBackend,
    SINGLE_PHASE_LIMIT_TOLERANCE,
    SELF_TEST_COMMON_NAMES,
    expected_self_test_names,
    self_test_name_contract,
    self_tests,
    self_tests_post_run,
)
from gribs.config.loader import (  # noqa: E402, F401
    _reject_unknown,
    _mapping,
    _section,
    _get,
    _number,
    _integer,
    _boolean,
    _string,
    _choice,
    _backend_choice,
    _string_array,
    DEFAULT_CEA_PYTHON_BLOCK,
    _MIGRATION_BACKEND_NOTE,
    migrate_v030_document,
    DEFAULT_TWO_PHASE_BLOCK,
    migrate_v043_document,
    migrate_v044_document,
    Configuration,
    _ROOT_KEYS,
    _TWO_PHASE_KEYS,
    _THERMO_ROOT_KEYS,
    _PROP_ROOT_KEYS,
    _BURN_LAW_KEYS,
    _ERO_BURN_KEYS,
    _REACTANT_KEYS,
    _GRAIN_KEYS,
    _NOZZLE_KEYS,
    _SEPARATION_KEYS,
    _THROAT_EROSION_KEYS,
    _ENV_KEYS,
    _IGNITER_KEYS,
    _PRESSURE_RANGE_KEYS,
    _CEA_PY_KEYS,
    _CEA_LEGACY_KEYS,
    _REFERENCE_CONSTANT_KEYS,
    _SOLVER_KEYS,
    _BLOWDOWN_KEYS,
    _OUTPUT_KEYS,
    _resolve_from_script,
    _load_document,
    _reject_nonfinite_values,
    configuration_from_document,
    load_json_configuration,
)
from gribs.cli import (  # noqa: E402, F401
    environment_info,
    build_parser,
    _print_backend_banner,
    _cea_theoretical_diagnostics,
    main,
)

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
# moved to src/gribs/config/__init__.py (A-6 unit 2; a package since unit 9); re-exported above


# ==============================================================================
# 2. DERIVED GEOMETRY HELPERS
# ==============================================================================
# moved to src/gribs/config/__init__.py (A-6 unit 2; a package since unit 9); re-exported above


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


# _is_nonfinite_number, require_finite_config, validate: moved to src/gribs/config/__init__.py (A-6
# unit 2; a package since unit 9); re-exported above. Readers of validate call gribs.config.validate
# (tests patch it there).


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
# moved to src/gribs/selftest.py (A-6 unit 8); re-exported above. The entry reads self_tests_post_run as
# gribs.selftest.self_tests_post_run (T4); the other names are bare re-exported reads.


# POST_RUN_MASS_TOLERANCE: moved to src/gribs/constants.py (A-6 unit 1); re-exported above


# ==============================================================================
# 10. PLOTTING (one single, carefully designed figure)
#     Panel layout unchanged from v0.3.1-alpha; v0.4.4 adds the condensed
#     mass fraction as a second axis of panel (a) and the condensed outflow
#     component to panel (c), and the key-results block reports the two-phase
#     model statistics.
# ==============================================================================
# moved to src/gribs/plotting.py (A-6 unit 7); re-exported above. The entry reads make_figure as
# gribs.plotting.make_figure (T4); the other names are bare re-exported reads.


# ==============================================================================
# 11. FILE OUTPUT
#     v0.5.0-alpha: the mass columns are renamed to their true meaning (total
#     product mass basis) - see the CSV_COLUMNS note below.
# ==============================================================================
# moved to src/gribs/reporting.py (A-6 unit 7); re-exported above. The entry reads summary_text and
# write_csv as gribs.reporting.<name> (T4); the other names are bare re-exported reads.


# ==============================================================================
# 12. JSON CONFIGURATION (v0.4 schema, backend-specific validation)
# ==============================================================================
# moved to src/gribs/config/loader.py (A-6 unit 9); re-exported above. Tests patch load_json_configuration
# there ("gribs.config.loader.load_json_configuration").


# ==============================================================================
# 13. CLI AND MAIN PROGRAM
# ==============================================================================
# moved to src/gribs/cli.py (A-6 unit 9); re-exported above (main is called below). Tests patch
# _print_backend_banner, environment_info and _cea_theoretical_diagnostics there ("gribs.cli.<name>").


# _CEA_DIAG_DISPLAY_NUMERIC_FIELDS, _CEA_DIAG_DISPLAY_TEXT_FIELDS, _cea_diag_value_ok, cea_diagnostics_status,
# _cea_diag_display: moved to src/gribs/reporting.py (A-6 unit 7, U-35 placement); re-exported above.


if __name__ == "__main__":
    sys.exit(main())
