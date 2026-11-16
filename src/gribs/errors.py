"""Exception hierarchy (A-6 unit 1). Moved verbatim from the entry, section 0b."""
from __future__ import annotations

from typing import Any, Dict, Optional


class GribsError(Exception):
    """Base class for all GRIBS exceptions."""


class ConfigurationError(GribsError):
    """JSON configuration is missing, malformed, or inconsistent."""


class ThermochemistryError(GribsError):
    """A thermochemistry backend could not produce physical gas properties."""


class ThermochemistryCacheError(ThermochemistryError):
    """A cached chamber-property table exists but is unusable.

    ``code`` and ``context`` are optional so existing cache failures retain
    their current exception behavior while selected failures can be consumed
    by scripts without parsing prose.
    """

    def __init__(self, message: str, *, code: Optional[str] = None,
                 context: Optional[Dict[str, Any]] = None):
        super().__init__(message)
        self.code = code
        self.context = dict(context or {})


class ThermochemistryInvariantError(ThermochemistryError):
    """A thermochemistry table violates a required mass-basis invariant."""

    code = "E_THERMO_MASS_BASIS_INVARIANT_FAILED"

    def __init__(self, message: str, *, context: Optional[Dict[str, Any]] = None):
        super().__init__(message)
        self.context = dict(context or {})


class GammaBelowUnityError(ThermochemistryError):
    """A tabulated thermochemistry isentropic exponent is not greater than one."""

    code = "E_GAMMA_S_BELOW_UNITY"

    def __init__(self, message: str, *, context: Optional[Dict[str, Any]] = None):
        super().__init__(message)
        self.context = dict(context or {})


class BackendUnavailableError(ThermochemistryError):
    """A selected thermochemistry backend cannot be initialized."""


class CEADatabaseError(ThermochemistryError):
    """The bundled NASA CEA database could not be pinned in this process.

    Raised with ``code = "E_CEA_DATABASE_ALREADY_INITIALIZED"`` when the
    ``cea`` module of this process was already initialized with a different
    ``thermo.lib`` / ``trans.lib`` before GRIBS took it (P5-1).
    """

    def __init__(self, message: str, *, code: str,
                 context: Optional[Dict[str, Any]] = None):
        super().__init__(message)
        self.code = code
        self.context = dict(context or {})
