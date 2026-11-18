"""Errors of the energy-mode state recovery (A-6 stage 2, C5; contract 164 section 2b-FIX-1 (d))."""
from gribs.errors import GribsError


class StateRecoveryError(GribsError):
    """The state (p, T, phase fraction) could not be recovered from (E_cv, V_free, B_e): a failed equilibrium solve,
    start points that disagree, incomplete phase-transition data, or a state outside the model's range."""


class NoRootOutsidePlateauError(StateRecoveryError):
    """No root of the energy residual outside a declared plateau: the energy is out of the searched temperature range,
    or the residual jumps where no phase transition is declared.  A plateau of a declared transition is never reported
    this way (contract 2b-FIX-1 (d)-3)."""


class InitialFillError(GribsError):
    """The initial fill at t = 0 is rejected: the element basis does not cover an element, the weights are negative or not
    unique, the composition differs from the main-propellant products (an air fill, for example), or the equilibrium of
    the fill could not be solved (contract 164 section 2b-FIX-2)."""
