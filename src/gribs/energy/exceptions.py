"""Errors of the energy-mode state recovery (A-6 stage 2, C5; contract 164 section 2b-FIX-1 (d))."""
from gribs.errors import GribsError


class StateRecoveryError(GribsError):
    """The state (p, T, phase fraction) could not be recovered from (E_cv, V_free, B_e): a failed equilibrium solve,
    start points that disagree, incomplete phase-transition data, or a state outside the model's range."""


class NoRootOutsidePlateauError(StateRecoveryError):
    """No root of the energy residual outside a declared plateau: the energy is out of the searched temperature range,
    or the residual jumps where no phase transition is declared.  A plateau of a declared transition is never reported
    this way (contract 2b-FIX-1 (d)-3)."""


class LatentHeatMismatchError(StateRecoveryError):
    """The latent heat L_input that the constant-cp data of the two phases of a declared transition give differs from the L_calc
    that the equilibrium provider's own plateau gives (the check at t = 0, C16), by more than
    energy_model.condensed.latent_heat_check_tolerance_rel.  An inconsistency of the input data of t = 0 (like InitialFillError):
    it depends on the provider's data as well as on the configuration, so it is checked at run time and not by --validate-config."""

    code = "E_ENERGY_LATENT_HEAT_MISMATCH"


class ExitOutOfRangeError(StateRecoveryError):
    """The numerical isentropic flow has no continuous solution at the area ratio of the nozzle: the mass flux G stops decreasing (the
    fold of the supersonic branch, at the exit of a phase-transition plateau) before G* / G reaches the area ratio (A-6 stage 2, U-H25-5
    group H-I; decision 109, K-H25-36, K-H25-61).  Raised from the right-hand side or an event function, so a run that meets the fold
    at any time of the integration stops (it leaves no result); the command line ends with code 4 like any StateRecoveryError."""

    code = "E_ENERGY_EXIT_OUT_OF_RANGE"


class InitialFillError(GribsError):
    """The initial fill at t = 0 is rejected: the element basis does not cover an element, the weights are negative or not
    unique, the composition differs from the main-propellant products (an air fill, for example), or the equilibrium of
    the fill could not be solved (contract 164 section 2b-FIX-2)."""
