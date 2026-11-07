"""Program identity and constants (A-6 unit 1). Moved verbatim from the entry, section 0, plus
POST_RUN_MASS_TOLERANCE (formerly in the self-test section)."""
from __future__ import annotations

import numpy as np

PROGRAM_NAME = "GRIBS"
PROGRAM_VERSION = "0.6.2-alpha"
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
#: condensed), EOS acts on the gas phase through Psi = R_eff*T0 (standard).
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
BACKEND_REFERENCE_CONSTANT = "reference_constant"
#: The complete set of selectable backends.  ``thermochemistry.backend`` is a
#: required configuration key whose value must be one of these; the program never
#: picks a backend on the user's behalf.
KNOWN_BACKENDS = (
    BACKEND_CEA_PYTHON,
    BACKEND_CEA_LEGACY_EXECUTABLE,
    BACKEND_REFERENCE_CONSTANT,
)
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
#: against (see the README, "Validation status").  Older majors do not provide the
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
#: Semantic cache identities are checked more tightly than ODE-level tests.
_MASS_BASIS_INVARIANT_REL_TOL = 1.0e-13
_MASS_BASIS_PHASE_SUM_ATOL = 1.0e-13
#: Scope thresholds for the A-5.1-6 gas-mass-fraction model classification.
MODEL_STATUS_VALIDATED_YG_THRESHOLD = 0.99
MODEL_STATUS_OUT_OF_SCOPE_YG_THRESHOLD = 0.65
#: Upper Yg_min bound of the high-condensate band (10/07 A-5.1-1 target band
#: Yg 0.65-0.80) for which reproduction of reported or experimental data
#: (gate G4b) has not been done.  It only selects the warning
#: W_HIGH_CONDENSATE_G4B_UNVALIDATED; model_status is unchanged.
HIGH_CONDENSATE_G4B_UNVALIDATED_YG_THRESHOLD = 0.80
G0 = 9.80665                      # standard gravity [m/s^2]
_TRAPZ = np.trapezoid if hasattr(np, "trapezoid") else np.trapz

_PHASE_TAGS = ("(cr)", "(L)", "(l)", "(gr)", "(s)", "(a)", "(b)", "(I)", "(II)",
               "(III)", "(IV)", "(V)")


#: Tolerance of the post-simulation mass-conservation tests.  The integrated
#: mass quadratures inherit the ODE error control (default rtol = 1e-10), so a relative
#: consistency of 1e-6 is far above the numerical noise floor while still
#: catching structural mass-basis errors (which are O(1) or O(Yc)).
POST_RUN_MASS_TOLERANCE = 1.0e-6
