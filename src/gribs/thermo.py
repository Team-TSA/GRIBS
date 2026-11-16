"""Thermochemistry (A-6 unit 3b): backends, the chamber-property cache, the backend state and the solver-facing
accessors. Moved verbatim from the entry: section 3 and, from post-processing, CEA_DENSITY_NOTE,
N_GAS_PER_KG_TOTAL_DEFINITION and _sonic_column_status (report 113 placement). One line differs: a relative
fcea2 executable path is looked up under _PROGRAM_DIR, the directory of the entry, as before the move.

The backend state is read and written only by _thermo_state() and _set_thermo() (ADR_A6_PACKAGE_SPLIT decision 2,
guard T10); it is one state for every loaded copy of the entry. Tests patch the names of this module here, e.g.
"gribs.thermo.thermochemistry"; the entry reads them as gribs.thermo.<name>.
"""
from __future__ import annotations

import hashlib
import importlib
import importlib.metadata as importlib_metadata
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.interpolate import PchipInterpolator

from gribs.constants import (
    PROGRAM_VERSION,
    SCHEMA_VERSION,
    TWO_PHASE_HOMOGENEOUS,
    TWO_PHASE_SINGLE_LEGACY,
    BACKEND_CEA_PYTHON,
    BACKEND_CEA_LEGACY_EXECUTABLE,
    BACKEND_REFERENCE_CONSTANT,
    KNOWN_BACKENDS,
    REMOVED_BACKENDS,
    CEA_PROJECT_URL,
    CEA_INSTALL_COMMAND,
    CEA_MIN_API_VERSION,
    CEA_TESTED_API_VERSION,
    R_UNIVERSAL,
    _PHASE_FRACTION_ROUNDING_TOL,
    _MASS_BASIS_INVARIANT_REL_TOL,
    _MASS_BASIS_PHASE_SUM_ATOL,
    G0,
)
from gribs.errors import (
    ConfigurationError,
    ThermochemistryError,
    ThermochemistryCacheError,
    ThermochemistryInvariantError,
    GammaBelowUnityError,
    BackendUnavailableError,
    CEADatabaseError,
)
from gribs.util import _report_value_text
from gribs.config import Config, require_finite_config

#: The program directory: this file is <dir>/src/gribs/thermo.py and the entry is in <dir> (the entry puts <dir>/src
#: on sys.path). A relative fcea2 executable path was resolved against the entry's directory before the move.
_PROGRAM_DIR = Path(__file__).resolve().parents[2]


# ---- section 3 of the entry: thermochemistry backends ----
#: Actionable message required whenever the official package is missing.
MISSING_CEA_PACKAGE_MESSAGE = f"""The official NASA CEA Python package is required for the
'{BACKEND_CEA_PYTHON}' backend.

Install it with:
    {CEA_INSTALL_COMMAND}

Official project:
    {CEA_PROJECT_URL}"""

_TABLE_REQUIRED_KEYS = ("pressure_Pa", "temperature_K", "gamma_s",
                        "equilibrium_molecular_weight_kg_kmol", "effective_gas_constant_J_kgK")
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

    def thermo_dump_table(self) -> Dict[str, List[float]]:
        """Return deterministic chamber-property columns for --dump-thermo."""
        table = getattr(self, "table", None)
        if not isinstance(table, dict):
            raise ThermochemistryError(
                f"Backend {self.name!r} does not provide a thermochemistry dump."
            )

        keys = (
            "pressure_Pa",
            "temperature_K",
            "gamma_s",
            "equilibrium_molecular_weight_kg_kmol",
            "effective_gas_constant_J_kgK",
            "gas_mass_fraction",
            "condensed_mass_fraction",
            "psi_J_kg",
        )
        return {
            key: list(table.get(key, []))
            for key in keys
        }

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
    interpolators in z = ln(p).  The derivatives of Psi = R_eff*T0 (two-phase
    mode) and Theta = Rg*T0 (single-phase legacy mode) are taken from the
    analytic derivative of the interpolators and converted with

        dPsi/dp = dPsi/dln(p) / p

    so that the ODE pressure equation never sees finite-difference noise and CEA is
    never called from inside solve_ivp.

    Cache schema 4 (v0.5.0-alpha, A-5.1-4) tabulates the total
    equilibrium-mixture mass basis: the molecular weight is
    M_CEA = 1 / EqSolution.n and the gas constant is R_eff = Ru / M_CEA.
    Schema 3 caches (obsolete gas-phase mass basis) and schema 2 caches
    (v0.4.3-alpha, no phase split) are rejected and never silently reused.

    Corrected mass basis (cache schema 4):

        M_CEA(p) = 1 / EqSolution.n      [kg/kmol, total mixture]
        R_eff(p) = Ru / M_CEA(p) = Yg(p) * R_gas(p)
        p/(rho*T) = R_eff

    The chamber EOS uses Psi = R_eff * T_used with no additional Yg factor.
    """

    cache_schema = 4
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
        self.cache_write_durability = None      # set only when fsync of the cache directory is unsupported

        if (c.cea_cache_enabled and not c.cea_rebuild_cache
                and self.cache_path.is_file()):
            table = self._read_cache()
            self.cache_used = True
        else:
            if c.cea_cache_enabled and not c.cea_rebuild_cache:
                legacy_cache = self._matching_legacy_schema3_cache()
                if legacy_cache is not None:
                    old_path, old_table = legacy_cache
                    self._raise_old_thermo_cache(
                        old_path,
                        old_table.get("schema"),
                        detection="matching_schema3_cache_key",
                        legacy_cache_key=self._legacy_schema3_cache_key(),
                    )
            table = self._build_table()
            self._validate_mass_basis_table(table, table_source="generated")
            if c.cea_cache_enabled:
                try:                                           # atomic cache write, see _cache_atomic_write_text
                    self.cache_root.mkdir(parents=True, exist_ok=True)
                except OSError as exc:
                    raise _cache_write_failed(_CacheWriteError(
                        "create_directory", "create_directory", self.cache_root, exc, False,
                        {"status": "not_created", "created": False})) from exc
                try:
                    self.cache_write_durability = _cache_atomic_write_text(
                        self.cache_path, json.dumps(table, indent=2, sort_keys=False) + "\n")
                except _CacheWriteError as exc:
                    raise _cache_write_failed(exc) from exc
                if self.cache_write_durability is not None:
                    print("WARNING: chamber-property cache written, but fsync of its directory is not supported "
                          f"here ({self.cache_write_durability['error']}); its durability after a crash or power "
                          "loss is not guaranteed", file=sys.stderr)
                self.cache_written = True
        self.table = table
        self._load_table(table)
        self.property_validation = self._validate_table()

    # -- cache key -----------------------------------------------------------
    #: Integer version of the mass-basis semantics.  Incremented when the
    #: physical meaning of M_CEA, R_eff, Psi or the nozzle mass-flow model
    #: changes.  Disjoint from ``cache_schema`` which tracks the file format.
    mass_basis_semantics_version = 1

    def cache_key_payload(self) -> Dict[str, Any]:
        """Everything that can change the tabulated properties.

        Solid constituent densities and the packing fraction affect rho_p only and
        are deliberately excluded (they are still reported in the provenance), so
        changing them does not invalidate the gas-property cache.
        """
        return {
            "cache_schema": self.cache_schema,
            "mass_basis_semantics_version": self.mass_basis_semantics_version,
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
            "two_phase_mode": str(self.c.two_phase_mode),
            "molecular_weight_definition": (
                "M_CEA = 1 / EqSolution.n [kg-total/kmol-gas]; "
                "R_eff = Ru / M_CEA = Yg * R_gas; "
                "Psi = R_eff * T_used (no additional Yg factor)"),
            "options": self.option_signature(),
        }

    def option_signature(self) -> Dict[str, Any]:
        return {}

    def _legacy_option_signature(self) -> Dict[str, Any]:
        """Option signature of the frozen schema-3 key (the current one unless a backend overrides it)."""
        return self.option_signature()

    def cache_key(self) -> str:
        payload = self.cache_key_payload()
        blob = json.dumps(payload, sort_keys=True, default=str).encode("utf-8")
        return hashlib.sha256(blob).hexdigest()[:24]

    def _legacy_schema3_cache_key(self) -> str:
        """Recompute the exact schema-3 key for this backend/configuration.

        Schema 4 intentionally changes the key, so a same-backend schema-3
        cache can otherwise be missed and silently rebuilt. This frozen payload
        mirrors the last schema-3 identity and lets us check only the one old
        cache that could have matched the current thermochemistry inputs.
        """
        payload = {
            "cache_schema": 3,
            "backend": self.name,
            "reactants": [[r.name, float(r.wt_percent), float(r.temperature_K)]
                          for r in self.reactants],
            "equilibrium": self.equilibrium_formulation,
            "pressure_min_Pa": float(self.c.p_fit_min),
            "pressure_max_Pa": float(self.c.p_fit_max),
            "pressure_points": int(self.c.cea_pressure_points),
            "temperature_efficiency": float(self.c.eta_T0),
            "phase_fractions_tabulated": bool(self.supports_two_phase),
            "molecular_weight_definition": (
                "gas_phase_M "
                + "(gas-phase M; total MW is diagnostic only since "
                + "v0.5.0-alpha)"
            ),
            "options": self._legacy_option_signature(),
        }
        blob = json.dumps(payload, sort_keys=True, default=str).encode("utf-8")
        return hashlib.sha256(blob).hexdigest()[:24]

    def _legacy_schema3_cache_path(self) -> Path:
        return self.cache_root / (
            f"{self.table_stem}_{self._legacy_schema3_cache_key()}.json")

    def _matching_legacy_schema3_cache(
        self,
    ) -> Optional[Tuple[Path, Dict[str, Any]]]:
        """Return the matching old cache, without flagging unrelated entries."""
        path = self._legacy_schema3_cache_path()
        if not path.is_file():
            return None
        try:
            table = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError):
            return None
        if not isinstance(table, dict) or table.get("schema") != 3:
            return None
        legacy_key = self._legacy_schema3_cache_key()
        if table.get("backend") not in (None, self.name):
            return None
        if table.get("cache_key") not in (None, legacy_key):
            return None
        return path, table

    def _raise_old_thermo_cache(
        self,
        cache_path: Path,
        schema_actual: Any,
        *,
        detection: str,
        legacy_cache_key: Optional[str] = None,
    ) -> None:
        detail = (
            "uses obsolete CEA mass-basis semantics"
            if schema_actual == 3 else "is obsolete"
        )
        context = {
            "backend": self.name,
            "cache_path": str(cache_path),
            "schema_actual": schema_actual,
            "schema_expected": self.cache_schema,
            "detection": detection,
            "cache_key_current": self.key,
        }
        if legacy_cache_key is not None:
            context["cache_key_schema3"] = legacy_cache_key
        raise ThermochemistryCacheError(
            f"E_OLD_THERMO_CACHE: Chamber-property cache schema "
            f"{schema_actual} {detail} and is never reused: {cache_path}\n"
            f"Expected schema: {self.cache_schema}; actual schema: "
            f"{schema_actual}.\n"
            "Rebuild it with --no-cache (or set 'rebuild_cache' to true), "
            "or remove the obsolete cache file.",
            code="E_OLD_THERMO_CACHE",
            context=context,
        )

    def _raise_old_cache_semantics_version(self, actual: Any) -> None:
        context = {
            "backend": self.name,
            "cache_path": str(self.cache_path),
            "field": "mass_basis_semantics_version",
            "value_actual": actual,
            "value_expected": self.mass_basis_semantics_version,
            "cache_key": self.key,
        }
        raise ThermochemistryCacheError(
            "E_OLD_THERMO_CACHE: chamber-property cache mass-basis semantics "
            f"version is {actual!r}; expected "
            f"{self.mass_basis_semantics_version}: {self.cache_path}\n"
            "Rebuild it with --no-cache (or set 'rebuild_cache' to true).",
            code="E_OLD_THERMO_CACHE",
            context=context,
        )

    def _raise_mass_basis_invariant(
        self,
        table_source: str,
        invariant: str,
        detail: str,
        **context_values: Any,
    ) -> None:
        context = {
            "backend": self.name,
            "cache_path": str(self.cache_path),
            "cache_key": self.key,
            "table_source": table_source,
            "invariant": invariant,
            **context_values,
        }
        raise ThermochemistryInvariantError(
            f"E_THERMO_MASS_BASIS_INVARIANT_FAILED: {detail}\n"
            f"Backend: {self.name}; cache: {self.cache_path}",
            context=context,
        )

    def _validate_gamma_s_above_unity(
        self,
        pressure: Sequence[float],
        gamma_s: Sequence[float],
        *,
        table_source: str,
    ) -> None:
        """Report non-admissible tabulated or interpolated gamma_s values."""
        p = np.asarray(pressure, dtype=float)
        gamma = np.asarray(gamma_s, dtype=float)
        if p.ndim != 1 or p.size == 0 or gamma.shape != p.shape:
            raise ThermochemistryError(
                f"Cannot validate gamma_s for {table_source} table: pressure "
                "and gamma_s must be non-empty, same-shaped one-dimensional arrays."
            )
        if not np.all(np.isfinite(gamma)):
            index = int(np.flatnonzero(~np.isfinite(gamma))[0])
            raise ThermochemistryError(
                f"Non-finite gamma_s at pressure-grid index {index} "
                f"({p[index]:.8g} Pa) in {table_source} table."
            )

        invalid = np.flatnonzero(gamma <= 1.0)
        if invalid.size == 0:
            return

        split_at = np.flatnonzero(np.diff(invalid) > 1) + 1
        intervals = [
            [float(p[int(run[0])]), float(p[int(run[-1])])]
            for run in np.split(invalid, split_at)
        ]
        invalid_gamma = gamma[invalid]
        context = {
            "backend": self.name,
            "cache_path": str(self.cache_path),
            "cache_key": self.key,
            "table_source": table_source,
            "cause": "The current isentropic property/nozzle relations require gamma_s > 1.",
            "gamma_s_requirement": "> 1",
            "gamma_s_threshold_exclusive": 1.0,
            "invalid_point_count": int(invalid.size),
            "invalid_sample_index_range": [int(invalid[0]), int(invalid[-1])],
            "invalid_pressure_intervals_Pa": intervals,
            "invalid_gamma_s_range": [
                float(np.min(invalid_gamma)), float(np.max(invalid_gamma))
            ],
            "pressure_scan_range_Pa": [float(p[0]), float(p[-1])],
            "configured_pressure_range_Pa": [
                float(getattr(self.c, "p_fit_min", p[0])),
                float(getattr(self.c, "p_fit_max", p[-1])),
            ],
        }
        interval_text = "; ".join(
            f"{low:.8g}–{high:.8g} Pa" for low, high in intervals
        )
        raise GammaBelowUnityError(
            "E_GAMMA_S_BELOW_UNITY: the current property/nozzle model requires "
            f"gamma_s > 1; the {table_source} check found {invalid.size} "
            f"sample(s) with gamma_s <= 1 over "
            f"{interval_text}. Do not clip these values; the selected "
            "equilibrium state lies outside the model's admissible gamma range.",
            context=context,
        )

    def _validate_mass_basis_table(
        self,
        table: Dict[str, Any],
        *,
        table_source: str,
    ) -> None:
        """Check persisted/build-table identities without altering valid values."""
        semantics_version_present = "mass_basis_semantics_version" in table
        semantics_version = table.get("mass_basis_semantics_version")
        if semantics_version_present:
            if (isinstance(semantics_version, bool)
                    or not isinstance(semantics_version, int)
                    or semantics_version != self.mass_basis_semantics_version):
                self._raise_old_cache_semantics_version(semantics_version)

        try:
            pressure = np.asarray(table["pressure_Pa"], dtype=float)
        except (KeyError, TypeError, ValueError, OverflowError) as exc:
            self._raise_mass_basis_invariant(
                table_source,
                "pressure_grid_numeric",
                f"pressure_Pa cannot be interpreted as a numeric array: {exc}",
            )
        if pressure.ndim != 1 or pressure.size == 0:
            self._raise_mass_basis_invariant(
                table_source,
                "pressure_grid_shape",
                "pressure_Pa must be a non-empty one-dimensional array.",
                actual_shape=list(pressure.shape),
            )
        if not np.all(np.isfinite(pressure)):
            index = int(np.flatnonzero(~np.isfinite(pressure))[0])
            self._raise_mass_basis_invariant(
                table_source,
                "pressure_grid_finite",
                f"pressure_Pa contains a non-finite value at index {index}.",
                index=index,
                pressure_Pa=repr(float(pressure[index])),
            )
        bad_pressure = pressure <= 0.0
        if np.any(bad_pressure):
            index = int(np.flatnonzero(bad_pressure)[0])
            self._raise_mass_basis_invariant(
                table_source,
                "pressure_grid_positive",
                "pressure_Pa must be positive at every grid point.",
                index=index,
                pressure_Pa=float(pressure[index]),
                expected="> 0",
            )
        nonincreasing = np.flatnonzero(np.diff(pressure) <= 0.0)
        if nonincreasing.size:
            index = int(nonincreasing[0]) + 1
            self._raise_mass_basis_invariant(
                table_source,
                "pressure_grid_strictly_increasing",
                f"pressure_Pa is not strictly increasing at index {index}.",
                index=index,
                pressure_Pa=float(pressure[index]),
                previous_pressure_Pa=float(pressure[index - 1]),
            )

        def numeric_column(name: str, *, required: bool = True):
            if name not in table:
                if required:
                    self._raise_mass_basis_invariant(
                        table_source,
                        f"column_present:{name}",
                        f"required mass-basis column {name!r} is missing.",
                        column=name,
                    )
                return None
            try:
                values = np.asarray(table[name], dtype=float)
            except (TypeError, ValueError, OverflowError) as exc:
                self._raise_mass_basis_invariant(
                    table_source,
                    f"column_numeric:{name}",
                    f"column {name!r} is not numeric: {exc}",
                    column=name,
                )
            if values.shape != pressure.shape:
                self._raise_mass_basis_invariant(
                    table_source,
                    f"column_shape:{name}",
                    f"column {name!r} does not match the pressure grid shape.",
                    column=name,
                    actual_shape=list(values.shape),
                    expected_shape=list(pressure.shape),
                )
            if not np.all(np.isfinite(values)):
                index = int(np.flatnonzero(~np.isfinite(values))[0])
                self._raise_mass_basis_invariant(
                    table_source,
                    f"column_finite:{name}",
                    f"column {name!r} contains a non-finite value at index {index}.",
                    column=name,
                    index=index,
                    pressure_Pa=float(pressure[index]),
                    value=repr(float(values[index])),
                )
            return values

        temperature = numeric_column("temperature_K") * float(self.c.eta_T0)
        gamma = numeric_column("gamma_s")
        molecular_weight = numeric_column("equilibrium_molecular_weight_kg_kmol")
        gas_constant = numeric_column("effective_gas_constant_J_kgK")
        for name, values in (
            ("temperature_K", temperature),
            ("equilibrium_molecular_weight_kg_kmol", molecular_weight),
            ("effective_gas_constant_J_kgK", gas_constant),
        ):
            bad = values <= 0.0
            if np.any(bad):
                index = int(np.flatnonzero(bad)[0])
                self._raise_mass_basis_invariant(
                    table_source,
                    f"positive:{name}",
                    f"{name} must be positive at every pressure-grid point.",
                    column=name,
                    index=index,
                    pressure_Pa=float(pressure[index]),
                    value=float(values[index]),
                    expected="> 0",
                )
        self._validate_gamma_s_above_unity(
            pressure,
            gamma,
            table_source=table_source,
        )

        with np.errstate(over="ignore", invalid="ignore"):
            rm_product = gas_constant * molecular_weight
            rm_relative_error = np.abs(rm_product / R_UNIVERSAL - 1.0)
        if not np.all(np.isfinite(rm_relative_error)):
            self._raise_mass_basis_invariant(
                table_source,
                "R_eff_times_M_CEA_equals_Ru",
                "R_eff * M_CEA produced a non-finite value.",
                tolerance=_MASS_BASIS_INVARIANT_REL_TOL,
            )
        worst_rm = int(np.argmax(rm_relative_error))
        if float(rm_relative_error[worst_rm]) > _MASS_BASIS_INVARIANT_REL_TOL:
            self._raise_mass_basis_invariant(
                table_source,
                "R_eff_times_M_CEA_equals_Ru",
                "R_eff * M_CEA differs from the universal gas constant.",
                index=worst_rm,
                pressure_Pa=float(pressure[worst_rm]),
                actual=float(rm_product[worst_rm]),
                expected=R_UNIVERSAL,
                relative_error=float(rm_relative_error[worst_rm]),
                tolerance=_MASS_BASIS_INVARIANT_REL_TOL,
            )

        if "n_gas_per_kg_total" in table or semantics_version_present:
            n_gas = numeric_column("n_gas_per_kg_total")
            with np.errstate(over="ignore", invalid="ignore"):
                mn_product = molecular_weight * n_gas
                mn_relative_error = np.abs(mn_product - 1.0)
            nonfinite_mn = np.flatnonzero(~np.isfinite(mn_relative_error))
            if nonfinite_mn.size:
                index = int(nonfinite_mn[0])
                self._raise_mass_basis_invariant(
                    table_source,
                    "M_CEA_times_n_gas_per_kg_total_equals_one",
                    "M_CEA * n_gas_per_kg_total produced a non-finite value.",
                    index=index,
                    pressure_Pa=float(pressure[index]),
                    actual=repr(float(mn_product[index])),
                    tolerance=_MASS_BASIS_INVARIANT_REL_TOL,
                )
            worst_mn = int(np.argmax(mn_relative_error))
            if float(mn_relative_error[worst_mn]) > _MASS_BASIS_INVARIANT_REL_TOL:
                self._raise_mass_basis_invariant(
                    table_source,
                    "M_CEA_times_n_gas_per_kg_total_equals_one",
                    "M_CEA * n_gas_per_kg_total differs from one.",
                    index=worst_mn,
                    pressure_Pa=float(pressure[worst_mn]),
                    actual=float(mn_product[worst_mn]),
                    expected=1.0,
                    relative_error=float(mn_relative_error[worst_mn]),
                    tolerance=_MASS_BASIS_INVARIANT_REL_TOL,
                )

        psi = numeric_column("psi_J_kg", required=self.supports_two_phase)
        if psi is not None:
            with np.errstate(over="ignore", invalid="ignore"):
                psi_expected = gas_constant * temperature
            nonfinite_psi_expected = np.flatnonzero(~np.isfinite(psi_expected))
            if nonfinite_psi_expected.size:
                index = int(nonfinite_psi_expected[0])
                self._raise_mass_basis_invariant(
                    table_source,
                    "psi_equals_R_eff_times_T_used",
                    "R_eff * T_used produced a non-finite value.",
                    index=index,
                    pressure_Pa=float(pressure[index]),
                    expected=repr(float(psi_expected[index])),
                    tolerance=_MASS_BASIS_INVARIANT_REL_TOL,
                )
            psi_relative_error = np.abs(psi - psi_expected) / np.maximum(
                np.abs(psi_expected), 1.0e-300)
            worst_psi = int(np.argmax(psi_relative_error))
            if float(psi_relative_error[worst_psi]) > _MASS_BASIS_INVARIANT_REL_TOL:
                self._raise_mass_basis_invariant(
                    table_source,
                    "psi_equals_R_eff_times_T_used",
                    "psi_J_kg differs from R_eff * T_used.",
                    index=worst_psi,
                    pressure_Pa=float(pressure[worst_psi]),
                    actual=float(psi[worst_psi]),
                    expected=float(psi_expected[worst_psi]),
                    relative_error=float(psi_relative_error[worst_psi]),
                    tolerance=_MASS_BASIS_INVARIANT_REL_TOL,
                )

        has_yg = "gas_mass_fraction" in table
        has_yc = "condensed_mass_fraction" in table
        if self.supports_two_phase and not (has_yg and has_yc):
            missing = [name for name, present in (
                ("gas_mass_fraction", has_yg),
                ("condensed_mass_fraction", has_yc),
            ) if not present]
            self._raise_mass_basis_invariant(
                table_source,
                "phase_fraction_columns_present",
                f"two-phase table is missing phase-fraction columns {missing}.",
                missing_columns=missing,
            )
        if has_yg != has_yc:
            self._raise_mass_basis_invariant(
                table_source,
                "phase_fraction_columns_pair",
                "gas_mass_fraction and condensed_mass_fraction must appear together.",
                gas_mass_fraction_present=has_yg,
                condensed_mass_fraction_present=has_yc,
            )
        if has_yg and has_yc:
            yg = numeric_column("gas_mass_fraction")
            yc = numeric_column("condensed_mass_fraction")
            bad_phase = (yg <= 0.0) | (yg > 1.0) | (yc < 0.0) | (yc >= 1.0)
            if np.any(bad_phase):
                index = int(np.flatnonzero(bad_phase)[0])
                self._raise_mass_basis_invariant(
                    table_source,
                    "phase_fraction_physical_ranges",
                    "require 0 < Yg <= 1 and 0 <= Yc < 1.",
                    index=index,
                    pressure_Pa=float(pressure[index]),
                    gas_mass_fraction=float(yg[index]),
                    condensed_mass_fraction=float(yc[index]),
                )
            phase_sum_error = np.abs(yg + yc - 1.0)
            worst_phase = int(np.argmax(phase_sum_error))
            if float(phase_sum_error[worst_phase]) > _MASS_BASIS_PHASE_SUM_ATOL:
                self._raise_mass_basis_invariant(
                    table_source,
                    "gas_plus_condensed_mass_fraction_equals_one",
                    "gas_mass_fraction + condensed_mass_fraction differs from one.",
                    index=worst_phase,
                    pressure_Pa=float(pressure[worst_phase]),
                    actual=float(yg[worst_phase] + yc[worst_phase]),
                    expected=1.0,
                    absolute_error=float(phase_sum_error[worst_phase]),
                    tolerance=_MASS_BASIS_PHASE_SUM_ATOL,
                )

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
        schema = table.get("schema")
        if schema is None:
            raise ThermochemistryCacheError(
                f"Chamber-property cache has no schema identity: "
                f"{self.cache_path}\n"
                "The A-5.1-4 CEA mass-basis correction requires an explicit "
                "integer cache schema.\n"
                "Rebuild it: set 'rebuild_cache' to true once, or delete the "
                "file, or run with --no-cache.")
        if isinstance(schema, bool) or not isinstance(schema, int):
            raise ThermochemistryCacheError(
                f"Chamber-property cache schema must be an integer, got "
                f"{schema!r}: {self.cache_path}\n"
                "Rebuild it: set 'rebuild_cache' to true once, or delete the "
                "file, or run with --no-cache.")
        if schema < self.cache_schema:
            self._raise_old_thermo_cache(
                self.cache_path,
                schema,
                detection="current_cache_path_schema_mismatch",
            )
        if schema != self.cache_schema:
            raise ThermochemistryCacheError(
                f"Chamber-property cache schema {schema} is incompatible with "
                f"required schema {self.cache_schema}: {self.cache_path}\n"
                "Rebuild it: set 'rebuild_cache' to true once, or delete the "
                "file, or run with --no-cache.")
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
        self._validate_mass_basis_table(table, table_source="cache")
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
        R = np.asarray(table["effective_gas_constant_J_kgK"], dtype=float)
        M = np.asarray(table["equilibrium_molecular_weight_kg_kmol"], dtype=float)
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
        if not np.all(np.isfinite(g)):
            raise ThermochemistryError(
                "The chamber-property table contains non-finite gamma_s values.")
        self._validate_gamma_s_above_unity(
            p,
            g,
            table_source="loaded",
        )

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
        # Psi = R_eff * T_used with T_used = eta_T0 * T_CEA.  Interpolated as
        # its own PCHIP column so that the analytic derivative dPsi/dp used by
        # the pressure ODE is consistent with the interpolated Psi itself
        # (never a finite difference, never a product of separately
        # interpolated factors).  R_eff already carries Yg, so no second
        # Yg factor may appear here.
        #
        # A-5.1-4 unified Psi and Theta: both are R_eff * T_used.
        # The separate theta interpolator is removed; theta_and_derivative()
        # delegates to psi_and_derivative().
        psi = R * T
        self.psi_tab = psi
        self.psii = PchipInterpolator(z, psi, extrapolate=True)
        self.dpsii = self.psii.derivative()
        self.Ygi = PchipInterpolator(z, Yg, extrapolate=True)

    def _validate_table(self) -> Dict[str, Any]:
        """Numerical sanity report for the table (recorded in summary.json)."""
        p = np.asarray(self.table["pressure_Pa"], dtype=float)
        T_used = (np.asarray(self.table["temperature_K"], dtype=float) * self.c.eta_T0)
        R = np.asarray(self.table["effective_gas_constant_J_kgK"], dtype=float)
        theta = R * T_used
        two_phase = (self.c.two_phase_mode == TWO_PHASE_HOMOGENEOUS
                     and self.supports_two_phase)
        if two_phase:
            Yg = np.asarray(self.table["gas_mass_fraction"], dtype=float)
            Yc = np.asarray(self.table["condensed_mass_fraction"], dtype=float)
            psi = theta
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
        # After A-5.1-4, Theta == Psi, so the Theta check uses the same
        # interpolator.  The separate factor_theta is retained for the
        # metadata/report only.
        factor_theta = factor
        # Psi consistency with R_eff*T at the table nodes (self-test 3 preview)
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
            "equilibrium_molecular_weight_min_kg_kmol":
                float(np.min(self.table["equilibrium_molecular_weight_kg_kmol"])),
            "equilibrium_molecular_weight_max_kg_kmol":
                float(np.max(self.table["equilibrium_molecular_weight_kg_kmol"])),
            # Deprecated aliases (misnamed: the values are M_CEA on the
            # total-mixture mass basis, not the gas-phase molecular weight).
            # Kept for existing summary readers; see deprecated_keys.
            "gas_phase_molecular_weight_min_kg_kmol":
                float(np.min(self.table["equilibrium_molecular_weight_kg_kmol"])),
            "gas_phase_molecular_weight_max_kg_kmol":
                float(np.max(self.table["equilibrium_molecular_weight_kg_kmol"])),
            "deprecated_keys": {
                "gas_phase_molecular_weight_min_kg_kmol":
                    "equilibrium_molecular_weight_min_kg_kmol",
                "gas_phase_molecular_weight_max_kg_kmol":
                    "equilibrium_molecular_weight_max_kg_kmol",
            },
            "gas_mass_fraction_min": float(np.min(Yg)),
            "gas_mass_fraction_max": float(np.max(Yg)),
            "condensed_mass_fraction_min": float(np.min(Yc)),
            "condensed_mass_fraction_max": float(np.max(Yc)),
            "theta_min_J_kg": float(np.min(theta)),
            "theta_max_J_kg": float(np.max(theta)),
            "psi_min_J_kg": float(np.min(psi)),
            "psi_max_J_kg": float(np.max(psi)),
            "psi_vs_ReffT_max_rel_error": psi_err,
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
                "Cached Psi column is inconsistent with R_eff*T0 (max rel. err "
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
        if not (_finite_positive(R) and _finite_positive(T) and math.isfinite(g)):
            raise ThermochemistryError(
                f"Non-physical interpolated gas properties at p = {p:.6g} Pa: "
                f"R = {R:.6g}, T0 = {T:.6g}, gamma_s = {g:.6g}.")
        if g <= 1.0:
            self._validate_gamma_s_above_unity(
                [pe],
                [g],
                table_source="interpolated",
            )
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
        """Psi = R_eff*T0 [J/kg] and its ANALYTIC derivative dPsi/dp [J/(kg Pa)].

        Theta and Psi are the same tabulated quantity R_eff*T_used since
        A-5.1-4, so every mode (including single_phase_legacy, where Yg := 1
        elsewhere) evaluates the same PCHIP column here.  This restores the
        pre-A-6 legacy path, whose separate Theta interpolator was built from
        the identical R*T nodes; delegating to theta_and_derivative() here
        recursed without bound (A-6 commit 3c24fbc regression).
        """
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
        """Theta = R_eff*T0 [J/kg] and dTheta/dp.

        After A-5.1-4, Theta == Psi (both are R_eff * T_used).
        This method delegates to psi_and_derivative() for backward
        compatibility with the single-phase legacy pathway.
        """
        return self.psi_and_derivative(p)

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
                "quantities": ["Rg", "T0", "gamma_s", "Yg", "Psi = R_eff*T0",
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
            "specific_gas_constant_definition": (
                "R_eff(p) = Ru / M_CEA(p) = Yg(p) * R_gas(p); "
                "M_CEA = 1 / EqSolution.n [kg-total/kmol-gas] on the "
                "total equilibrium mixture mass basis"),
            "phase_fraction_definition": ("Yg(p) = sum of the gas-species mass "
                                          "fractions, Yc(p) = sum of the "
                                          "condensed-species mass fractions, from "
                                          "the official CEA phase ordering "
                                          "(gas species first, then condensed)"),
            "mass_basis_semantics_version": self.mass_basis_semantics_version,
            "M_definition": "M_CEA = 1 / EqSolution.n [kg-total/kmol-gas]",
            "MW_definition": "sol.MW = 1 / sum_all(nj) [kg-total/kmol-all] (diagnostic only, never used in EOS)",
            "R_definition": "R_eff = Ru / M_CEA [J/(kg-total K)]",
            "psi_definition": "Psi = R_eff * T_used [J/kg-total], with T_used = eta_T0 * T_CEA; no additional Yg factor",
            "density_definition": "rho = p * M_CEA / (Ru * T) [kg-total/m^3] (total mixture density from CEA)",
            "gamma_s_definition": "equilibrium mixture isentropic exponent (dimensionless, from CEA)",
            "two_phase_capabilities": {
                "equilibrium_phase_fractions": bool(self.supports_two_phase),
                "particle_slip": False,
                "finite_rate_condensation": False,
                "condensed_volume": False,
                "phase_density_lag": False,
            },
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
            _PROGRAM_DIR / p,
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
            "equilibrium_molecular_weight_kg_kmol": M.tolist(),
            "effective_gas_constant_J_kgK": R.tolist(),
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
W_EXTERNAL_CEA_DATABASE_IGNORED = "W_EXTERNAL_CEA_DATABASE_IGNORED"
E_CEA_DATABASE_ALREADY_INITIALIZED = "E_CEA_DATABASE_ALREADY_INITIALIZED"


def import_cea_without_auto_init():
    """Import ``cea`` without letting the import load a database (P5-1).

    ``import cea`` (3.3.x) initializes the library at import time from the first
    ``thermo.lib`` / ``trans.lib`` it finds: the current working directory, then
    ``CEA_DATA_DIR``, then the package data.  That choice is invisible to GRIBS,
    and a damaged ``trans.lib`` found there ends the process inside the import,
    even with transport off.  ``CEA_SKIP_INIT`` is read only by the first import
    in a process, so it is set just for that import and restored afterwards
    (child processes do not inherit it).  The database is then loaded
    explicitly by :func:`pin_bundled_cea_database`.  When ``cea`` is already
    imported, the existing module is returned unchanged.
    """
    if "cea" in sys.modules:
        return importlib.import_module("cea")
    previous = os.environ.get("CEA_SKIP_INIT")
    os.environ["CEA_SKIP_INIT"] = "1"
    try:
        return importlib.import_module("cea")
    finally:
        if previous is None:
            os.environ.pop("CEA_SKIP_INIT", None)
        else:
            os.environ["CEA_SKIP_INIT"] = previous


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _external_cea_database_files(bundled_dir: Path) -> List[Dict[str, str]]:
    """Database files that ``import cea`` would have preferred to the package data.

    Only their existence is checked; they are never opened.
    """
    places: List[Tuple[str, Path]] = [("cwd", Path.cwd())]
    env_dir = os.environ.get("CEA_DATA_DIR")
    if env_dir:
        places.append(("CEA_DATA_DIR", Path(env_dir)))
    try:
        bundled_real = os.path.realpath(bundled_dir)
    except OSError:                                # pragma: no cover - defensive
        bundled_real = str(bundled_dir)
    found: List[Dict[str, str]] = []
    seen = set()
    for location, directory in places:
        for name in ("thermo.lib", "trans.lib"):
            path = directory / name
            try:
                if not path.is_file():
                    continue
                real = os.path.realpath(path)
            except OSError:
                continue
            if os.path.dirname(real) == bundled_real or real in seen:
                continue
            seen.add(real)
            found.append({"location": location, "name": name, "path": str(path.absolute())})
    return found


def pin_bundled_cea_database(cea_module: Any, transport: bool) -> Dict[str, Any]:
    """Load the database shipped inside the official package and return its provenance (P5-1).

    * ``thermo.lib`` of the package is loaded with ``cea.init_thermo``;
      ``trans.lib`` of the package only when ``transport`` is true (with
      transport off no ``trans.lib`` is read at all).
    * ``thermo.lib`` / ``trans.lib`` in the working directory or in
      ``CEA_DATA_DIR`` are ignored and listed in ``ignored_external_files``
      (reported as ``W_EXTERNAL_CEA_DATABASE_IGNORED``).
    * When ``cea`` was already initialized before GRIBS took it, the same calls
      are made: cea accepts them silently only if the loaded file is the
      package file, and raises ``ValueError`` otherwise.  GRIBS then stops with
      ``E_CEA_DATABASE_ALREADY_INITIALIZED`` instead of using the other
      database.

    The provenance holds absolute paths for display; the cache key uses only
    the content digests (:meth:`CEAPythonBackend.option_signature`).
    """
    bundled_dir = Path(os.path.dirname(os.path.abspath(cea_module.__file__))) / "data"
    thermo_path = bundled_dir / "thermo.lib"
    trans_path = bundled_dir / "trans.lib"
    needed = [thermo_path] + ([trans_path] if transport else [])
    missing = [str(p) for p in needed if not p.is_file()]
    if missing:
        raise BackendUnavailableError(
            "The database files of the official NASA CEA Python package are missing:\n"
            + "\n".join(f"  {p}" for p in missing)
            + f"\nReinstall the package:\n    {CEA_INSTALL_COMMAND}")
    ignored = _external_cea_database_files(bundled_dir)
    already = bool(cea_module.is_initialized())
    loads = [("thermo.lib", cea_module.init_thermo, thermo_path)]
    if transport:
        loads.append(("trans.lib", cea_module.init_trans, trans_path))
    for name, loader, path in loads:
        try:
            loader(str(path))
        except ValueError as exc:
            context = {"database_file": name, "required_path": str(path),
                       "process_already_initialized": already,
                       "ignored_external_files": ignored,
                       "cea_error": f"{type(exc).__name__}: {exc}"}
            if already:
                raise CEADatabaseError(
                    f"{E_CEA_DATABASE_ALREADY_INITIALIZED}: the 'cea' module of this process was "
                    f"already initialized with a {name} other than the one of the official package "
                    f"({path}); cea cannot load another database in the same process, and GRIBS does "
                    f"not use a database it has not selected.\n"
                    "This happens when 'cea' was imported before GRIBS while a thermo.lib / trans.lib "
                    "was present in the working directory or in CEA_DATA_DIR (cea prints the loaded "
                    "path on standard output). Run GRIBS in a new process, or import cea only "
                    "through GRIBS.",
                    code=E_CEA_DATABASE_ALREADY_INITIALIZED, context=context) from exc
            raise BackendUnavailableError(
                f"The official NASA CEA Python package could not load its own {name} "
                f"({path}): {type(exc).__name__}: {exc}") from exc
    thermo_sha = _sha256_file(thermo_path)
    trans_sha = _sha256_file(trans_path) if transport else None
    if ignored:
        listing = ", ".join(f"{f['path']} ({f['location']})" for f in ignored)
        print(f"WARNING [{W_EXTERNAL_CEA_DATABASE_IGNORED}]: ignored {listing}; using the database "
              f"of the official CEA package ({thermo_path})", file=sys.stderr)
    return {
        "selection": "bundled",
        "thermo_lib_path": str(thermo_path),
        "thermo_lib_size_bytes": thermo_path.stat().st_size,
        "thermo_lib_sha256": thermo_sha,
        "trans_lib_loaded_by_gribs": bool(transport),
        "trans_lib_path": str(trans_path) if transport else None,
        "trans_lib_size_bytes": trans_path.stat().st_size if transport else None,
        "trans_lib_sha256": trans_sha,
        "process_already_initialized": already,
        "ignored_external_files": ignored,
    }


def import_official_cea_module():
    """Import and validate the official NASA CEA Python package.

    The module is imported without its automatic database initialization
    (:func:`import_cea_without_auto_init`); the backend loads the package
    database explicitly with :func:`pin_bundled_cea_database`.

    Raises
    ------
    BackendUnavailableError
        with an actionable message when the package is missing, when a different
        package named ``cea`` is installed, or when the installed official API
        version is not supported.
    """
    try:
        module = import_cea_without_auto_init()
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
    """Pre-P5-1 identity of the package database (16-digit digest prefixes and paths).

    Used only to recompute the frozen schema-3 cache key
    (:meth:`CEAPythonBackend._legacy_option_signature`).  The current cache key
    and provenance come from :func:`pin_bundled_cea_database`.
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
    * Molecular weight: ``EqSolution.M`` is M_CEA = 1 / EqSolution.n, the
      molecular weight of the TOTAL equilibrium mixture in kg/kmol
      (= g/mol) on the total-mixture mass basis - i.e. the same quantity the
      legacy ``plot m`` column provided.  The v0.4.3 option ``total_MW`` is
      REJECTED in v0.4.4: ``EqSolution.MW`` is tabulated as a diagnostic
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
        self.cea_database = pin_bundled_cea_database(self.cea, bool(c.cea_py_transport))
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
                "  The chamber equation of state and the nozzle model use "
                "M_CEA = EqSolution.M = 1/EqSolution.n [kg-total/kmol-gas] "
                "(= M_gas/Yg, total-mixture mass per kmol of gas), so that "
                "R_eff = Ru/M_CEA = Yg*R_gas. EqSolution.MW also counts "
                "condensed moles and is physically inconsistent there "
                "whenever condensed products exist.\n"
                "  Set molecular_weight = 'gas_phase_M' (compatibility token "
                "name; it selects M_CEA).  EqSolution.MW is still tabulated "
                "and reported as the diagnostic "
                "'total_molecular_weight_kg_kmol'.")
        if c.cea_py_molecular_weight != "gas_phase_M":
            raise ConfigurationError(
                "thermochemistry.cea_python.molecular_weight must be "
                "'gas_phase_M' (the only accepted value since v0.5.0-alpha).")
        if c.cea_py_products_from_reactants and c.cea_py_product_species:
            raise ConfigurationError(
                "thermochemistry.cea_python: 'product_species' must be empty when "
                "'products_from_reactants' is true.")

    def option_signature(self) -> Dict[str, Any]:
        """Backend part of the cache key: content digests only, no absolute paths (P5-1)."""
        c = self.c
        return {
            "cea_package_version": str(getattr(self.cea, "__version__", "unknown")),
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
            "equilibrium_type": "HP",
            "pressure_unit_for_api": "bar",
            "enthalpy_convention": "h/R [K]",
            "database": "official package thermo.lib (trans.lib only with transport)",
            "database_identity": {
                "selection": self.cea_database["selection"],
                "thermo_lib_sha256": self.cea_database["thermo_lib_sha256"],
                "trans_lib_sha256": self.cea_database["trans_lib_sha256"],
            },
        }

    def _legacy_option_signature(self) -> Dict[str, Any]:
        """The pre-P5-1 option signature, frozen for the schema-3 cache key."""
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
            "database_identity": _cea_database_identity(self.cea),
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
        n_gas_out: List[float] = []
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

        for grid_index, p_pa in enumerate(pressures):
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
            M = float(solution.M)          # M_CEA = 1/n [kg-total/kmol-gas] (= M_gas/Yg)
            MW = float(solution.MW)        # total molecular weight (diagnostic)
            gamma = float(solution.gamma_s)
            if not _finite_positive(T):
                raise ThermochemistryError(
                    f"CEA returned a non-physical chamber temperature "
                    f"T = {T!r} K at p = {p_pa:.8g} Pa.")
            if not _finite_positive(M):
                raise ThermochemistryError(
                    f"CEA returned a non-physical equilibrium molecular weight "
                    f"M_CEA = EqSolution.M = {M!r} kg-total/kmol-gas at "
                    f"p = {p_pa:.8g} Pa.")
            if not _finite_positive(MW):
                raise ThermochemistryError(
                    f"CEA returned a non-physical total molecular weight "
                    f"MW = {MW!r} kg/kmol at p = {p_pa:.8g} Pa.")
            if not math.isfinite(gamma):
                raise ThermochemistryError(
                    f"CEA returned a non-finite isentropic exponent "
                    f"gamma_s = {gamma!r} at p = {p_pa:.8g} Pa.")
            if gamma <= 1.0:
                # Keep table generation fail-fast but report the scanned range.
                self._validate_gamma_s_above_unity(
                    pressures[:grid_index + 1],
                    g_out + [gamma],
                    table_source="generated",
                )

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
            n_gas_out.append(float(solution.n))
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
        # Psi = R_eff * (eta_T0 * T_CEA); stored for a cache round-trip
        # consistency check - it is recomputed from R, T at load time.
        psi_out = [r * c.eta_T0 * t for r, t in zip(R_out, T_out)]
        return {
            "schema": self.cache_schema,
            "backend": self.name,
            "cache_key": self.key,
            "mass_basis_semantics_version": self.mass_basis_semantics_version,
            "M_definition": "M_CEA = EqSolution.M = 1 / EqSolution.n [kg-total/kmol-gas]",
            "R_definition": "R_eff = Ru / M_CEA [J/(kg-total K)]",
            "psi_definition": (
                "Psi = R_eff * T_used [J/kg-total], T_used = eta_T0 * T_CEA; "
                "no additional Yg factor"),
            "source": (f"official NASA CEA Python package {getattr(cea, '__version__', '?')} "
                       f"(HP equilibrium)"),
            "cea_version": str(getattr(cea, "__version__", "unknown")),
            "cea_library_version": str(cea.lib_version()),
            "cea_module_path": str(getattr(cea, "__file__", "unknown")),
            "pressure_Pa": pressures.tolist(),
            "temperature_K": T_out,
            "gamma_s": g_out,
            "equilibrium_molecular_weight_kg_kmol": M_out,
            "n_gas_per_kg_total": n_gas_out,
            "effective_gas_constant_J_kgK": R_out,
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
            "n_gas_per_kg_total_definition": (
                "EqSolution.n: total gas-phase moles per kg of TOTAL equilibrium "
                "products [kmol-gas/kg-total]; M_CEA = 1/n"),
            "molecular_weight_definition": (
                "gas_phase_M compatibility token (name kept; not the gas-phase "
                "mean molecular weight); EqSolution.M is M_CEA on the "
                "total equilibrium mixture mass basis; EqSolution.MW is "
                "tabulated as a diagnostic only"),
            "product_species_count": int(product_mixture.num_species),
            "num_gas": int(solver.num_gas),
            "num_condensed": int(solver.num_condensed),
            "num_elements": int(solver.num_elements),
            "diagnostics": {
                "density_kg_m3": rho_out,
                "density_note": CEA_DENSITY_NOTE,
                "cp_eq_kJ_kgK": cp_out,
                "sonic_velocity_m_s": a_out,
                "sonic_velocity_status": (
                    "UNAVAILABLE" if all(not math.isfinite(v) for v in a_out)
                    else ("AVAILABLE" if all(math.isfinite(v) for v in a_out)
                          else "PARTIAL")),
                "sonic_velocity_note": ("read from EqSolution.sonic_velocity when "
                                        "that attribute exists; CEA 3.3.4 EqSolution "
                                        "has no such attribute, so the column is NaN "
                                        "and status is UNAVAILABLE. RocketSolution "
                                        "station values are not used here"),
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
            "molecular_weight_definition": (
                "EqSolution.M is M_CEA = 1 / EqSolution.n on the total "
                "equilibrium mixture mass basis (kg/kmol); EqSolution.MW is "
                "diagnostic only"),
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
            "database": ("thermo.lib of the official package, loaded explicitly; trans.lib of the "
                         "package only with transport; thermo.lib / trans.lib in the working "
                         "directory or CEA_DATA_DIR are ignored"),
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
        md["semantic_definitions"] = self.semantic_definitions()
        return md

    def semantic_definitions(self) -> Dict[str, Any]:
        """Current definitions of the CEA diagnostic columns (cleanup P3a).

        Reported by the running code, so summary.json / --dump-thermo carry
        the current definitions even when the table was read from a cache
        file whose own diagnostic text differs (same cache key; the cache
        file is never rewritten).  Only observable facts about the cache are
        reported (which fields are missing, whether the stored density_note
        equals the current one); the cache's creation date is NOT inferred.
        """
        diag = self.table.get("diagnostics", {})
        if not isinstance(diag, dict):
            diag = {}
        n_grid = len(self.table.get("pressure_Pa", []) or [])
        sonic_status, sonic_reason = _sonic_column_status(
            diag.get("sonic_velocity_m_s", None), n_grid)
        checked = ("sonic_velocity_status", "sonic_velocity_note", "density_note")
        missing = [k for k in checked if k not in diag]
        if "n_gas_per_kg_total_definition" not in self.table:
            missing.append("n_gas_per_kg_total_definition")
        return {
            "density_definition": CEA_DENSITY_NOTE,
            "n_gas_per_kg_total_definition": N_GAS_PER_KG_TOTAL_DEFINITION,
            "sonic_velocity_status": sonic_status,
            "sonic_velocity_status_reason": sonic_reason,
            "table_source": "cache" if self.cache_used else "generated_this_run",
            "table_fields_missing": missing,
            "table_density_note_matches_current": (
                diag.get("density_note") == CEA_DENSITY_NOTE),
            "note": ("definitions above are from the running code; a cache file "
                     "keeps the diagnostic text written by its writer"),
        }

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
# 3.5 Constant-property analytic reference backend
# ------------------------------------------------------------------------------
class ReferenceConstantBackend(ThermochemistryBackend):
    """Deterministic constant-property backend for Phase A-5 references.

    This backend performs no CEA calculation and creates no property cache.
    The gas is single-phase, with Yg = 1 and Yc = 0.  R, T0, gamma, Theta,
    and Psi are independent of pressure, so their pressure derivatives vanish
    exactly.  The backend is intended for analytic Level 0/1 verification, not
    as a combustion-chemistry model.
    """

    name = BACKEND_REFERENCE_CONSTANT
    uses_cea = False
    supports_two_phase = True
    equilibrium_formulation = (
        "single-phase constant-property analytic reference"
    )

    def __init__(
        self,
        c: Config,
        reactants: Sequence[ReactantSpec],
        packing_fraction: float,
    ):
        self.c = c
        self.reactants = list(reactants)
        self.packing_fraction = float(packing_fraction)
        self.ideal_mixture_density = ideal_mixture_density(self.reactants)
        self.bulk_mixture_density = (
            self.ideal_mixture_density * self.packing_fraction
        )

        self.pmin = float(c.p_fit_min)
        self.pmax = float(c.p_fit_max)
        self._temperature = float(c.reference_temperature_K) * float(c.eta_T0)
        self._gamma = float(c.reference_gamma)
        self._molecular_weight = float(
            c.reference_molecular_weight_kg_kmol
        )
        self._gas_constant = R_UNIVERSAL / self._molecular_weight
        self._theta = self._gas_constant * self._temperature

        self.cache_path = None
        self.cache_used = False
        self.cache_written = False
        # The analytic reference is a single gas phase at every pressure
        # (phase_fractions() returns (1, 0)); the phase-fraction range is
        # stated explicitly so the two-phase applicability check never has
        # to assume a default.
        self.property_validation = {
            "passed": True,
            "constant_properties": True,
            "pressure_independent": True,
            "gas_mass_fraction_min": 1.0,
            "gas_mass_fraction_max": 1.0,
            "condensed_mass_fraction_min": 0.0,
            "condensed_mass_fraction_max": 0.0,
        }

    def props(self, pressure_pa: float) -> Tuple[float, float, float]:
        return self._gas_constant, self._temperature, self._gamma

    def theta_and_derivative(
        self,
        pressure_pa: float,
    ) -> Tuple[float, float]:
        return self._theta, 0.0

    def psi_and_derivative(
        self,
        pressure_pa: float,
    ) -> Tuple[float, float]:
        return self._theta, 0.0

    def phase_fractions(self, pressure_pa: float) -> Tuple[float, float]:
        return 1.0, 0.0

    def is_extrapolated(self, pressure_pa: float) -> bool:
        return not (self.pmin <= float(pressure_pa) <= self.pmax)

    @property
    def outside_range_policy(self) -> str:
        return self.c.property_policy

    def metadata(self) -> Dict[str, Any]:
        return {
            "backend": self.name,
            "uses_cea": False,
            "constant_properties": True,
            "analytic_reference_backend": True,
            "temperature_K": self._temperature,
            "gamma": self._gamma,
            "molecular_weight_kg_kmol": self._molecular_weight,
            "effective_gas_constant_J_kgK": self._gas_constant,
            "gas_mass_fraction": 1.0,
            "condensed_mass_fraction": 0.0,
            "equilibrium_formulation": self.equilibrium_formulation,
            "pressure_range_Pa": [self.pmin, self.pmax],
            "outside_range_policy": self.outside_range_policy,
            "cache_enabled": False,
            "ideal_mixture_density_kg_m3": self.ideal_mixture_density,
            "packing_fraction": self.packing_fraction,
            "bulk_propellant_density_kg_m3": self.bulk_mixture_density,
        }

    def thermo_dump_table(self) -> Dict[str, List[float]]:
        """Return a two-point representation of the analytic constants."""
        pressures = [self.pmin, self.pmax]
        return {
            "pressure_Pa": pressures,
            "temperature_K": [self._temperature, self._temperature],
            "gamma_s": [self._gamma, self._gamma],
            "equilibrium_molecular_weight_kg_kmol": [
                self._molecular_weight,
                self._molecular_weight,
            ],
            "effective_gas_constant_J_kgK": [
                self._gas_constant,
                self._gas_constant,
            ],
            "gas_mass_fraction": [1.0, 1.0],
            "condensed_mass_fraction": [0.0, 0.0],
            "psi_J_kg": [self._theta, self._theta],
        }


# ------------------------------------------------------------------------------
# 3.6 Backend factory and solver-facing accessors
# ------------------------------------------------------------------------------
# NOTE: _THERMO is a module-level global accessed by gas_props(), psi_and_deriv(),
# phase_fractions(), theta_and_deriv(), is_extrapolated() — all called from deep
# inside the solver (make_rhs → run_model).  Removing it requires threading a
# ``backend`` parameter through 64+ call sites.  The _ThermoSwap context manager
# now provides safe temporary replacement for self_tests().  Full removal is
# tracked as A-6 deferred (future PR).
_THERMO: Optional[ThermochemistryBackend] = None


def _thermo_state() -> Optional[ThermochemistryBackend]:
    """The installed backend, or None before initialization.

    _thermo_state() and _set_thermo() are the only code that reads or writes _THERMO (ADR_A6_PACKAGE_SPLIT
    decision 2, A-6 unit 3a); every other reader, the tests and the conftest state guard go through them.
    """
    return _THERMO


def _set_thermo(backend: Optional[ThermochemistryBackend]) -> Optional[ThermochemistryBackend]:
    """Install `backend` (None clears it) and return the backend that was installed before."""
    global _THERMO
    previous = _THERMO
    _THERMO = backend
    return previous


class _ThermoSwap:
    """Context manager for temporarily replacing the global thermo backend.

    Usage::

        with _ThermoSwap(synthetic_backend):
            # _THERMO is synthetic_backend here
        # _THERMO is restored to its previous value here

    This replaces the manual ``saved = _THERMO; _THERMO = x; ...; _THERMO = saved``
    pattern and guarantees restoration even if an exception occurs.
    """

    def __init__(self, replacement: ThermochemistryBackend):
        self._replacement = replacement
        self._saved: Optional[ThermochemistryBackend] = None

    def __enter__(self):
        self._saved = _set_thermo(self._replacement)
        return self._replacement

    def __exit__(self, exc_type, exc_val, exc_tb):
        _set_thermo(self._saved)
        return False  # do not suppress exceptions


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
    if c.thermo_backend == BACKEND_REFERENCE_CONSTANT:
        return ReferenceConstantBackend(c, reactants, packing_fraction)
    if c.thermo_backend in REMOVED_BACKENDS:
        raise ConfigurationError(
            f"thermochemistry.backend: {c.thermo_backend!r} - "
            f"{REMOVED_BACKENDS[c.thermo_backend]}")
    raise ConfigurationError(
        f"thermochemistry.backend: unknown backend {c.thermo_backend!r}. "
        f"Set it explicitly to one of: {', '.join(KNOWN_BACKENDS)}.")


class _CacheWriteError(Exception):
    """Internal: failing phase / step of a cache write, plus the outcome of the temporary-file cleanup."""

    def __init__(self, phase, step, path, original, replaced, cleanup):
        super().__init__(f"{phase}/{step}: {type(original).__name__}: {original}")
        self.phase, self.step, self.path = phase, step, str(path)
        self.original, self.replaced, self.cleanup = original, bool(replaced), dict(cleanup or {})


_CACHE_CLEANUP_TEXT = {
    "not_created": "no temporary file had been created",
    "removed": "this run's temporary file was removed",
    "already_absent": "this run's temporary file was no longer present",
    "left_untouched_not_this_runs_file": ("the temporary name no longer referred to this run's file; whatever is "
                                          "there now was left untouched"),
    "not_applicable": "the temporary file itself had been renamed onto the cache name",
    "left_unidentified": ("this run created a temporary file but could not identify it, so it was NOT removed; "
                          "it may remain"),
    "present": "this run's temporary name was still present (no removal was attempted)",
}


def _cache_write_failed(err):
    c = err.cleanup or {}
    if c.get("status") == "failed":
        cleanup = (f"this run's temporary file could NOT be removed ({c.get('error')}); it was left at "
                   f"{c.get('path')}")
    elif c.get("status") == "left_unidentified":
        cleanup = f"{_CACHE_CLEANUP_TEXT['left_unidentified']} (at {c.get('path')})"
    else:
        cleanup = _CACHE_CLEANUP_TEXT.get(c.get("status"), "temporary-file state unknown")
    if err.phase == "adopt_existing":
        state = ("this run did not change the cache name: an identical cache published by another run was adopted, "
                 "but this run could not finish its own steps (temp removal / directory sync / close)")
    elif err.replaced:
        state = ("the cache name was changed by this run (published or replaced), so any previous cache is gone; "
                 "the new file's identity or durability could not be confirmed")
    elif err.phase == "create_directory":
        state = "the cache directory could not be created; no cache file was written"
    else:
        state = "this run did not replace the cache file (an existing cache, if any, was not modified by this run)"
    exc = err.original
    return ThermochemistryCacheError(
        f"Chamber-property cache could not be written ({err.phase}, {err.step}): {err.path}\n"
        f"Original error: {type(exc).__name__}: {exc}\n{state}; {cleanup}.",
        code="E_THERMO_CACHE_WRITE_FAILED",
        context={"phase": err.phase, "step": err.step, "path": err.path,
                 "error": f"{type(exc).__name__}: {exc}", "errno": getattr(exc, "errno", None),
                 "replaced": err.replaced, "temp_cleanup": c})


def _cache_rename_noreplace(src, dst, dir_fd) -> None:
    """renameat2(dir_fd, src, dir_fd, dst, RENAME_NOREPLACE) through the C library (Linux, glibc 2.28 or later).
    Raises FileExistsError if dst exists, OSError(ENOSYS) if the call is not available, OSError(errno) otherwise."""
    import ctypes
    import errno
    try:
        fn = ctypes.CDLL(None, use_errno=True).renameat2
    except (OSError, AttributeError) as exc:
        raise OSError(errno.ENOSYS, "renameat2 is not available in this C library") from exc
    fn.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_uint]
    fn.restype = ctypes.c_int
    if fn(dir_fd, os.fsencode(src), dir_fd, os.fsencode(dst), 1) != 0:   # 1 = RENAME_NOREPLACE
        e = ctypes.get_errno()
        raise OSError(e, os.strerror(e), dst)


def _cache_write_checkpoint(stage, dir_fd, temp_name):
    """No-op test seam called at the stages 'temp_created', 'attributes_set', 'before_replace' and
    'after_replace'; tests replace it to synchronise writers or inject failures.  It does nothing in a run."""


def _cache_atomic_write_text(path, text: str):
    """Write a chamber-property cache file atomically.  Returns None, or a dict when fsync of the directory is
    unsupported (ENOTSUP /
    EOPNOTSUPP): durability not guaranteed.  Raises _CacheWriteError otherwise.
    Assumes a trusted cache directory (writable only by trusted principals); the name checks narrow, but do not
    close, rename races.  Concurrent writers of one key are assumed to write identical bytes."""
    import errno
    import stat as _stat
    import uuid
    UNSUP = {errno.ENOTSUP, getattr(errno, "EOPNOTSUPP", errno.ENOTSUP)}
    NOLINK = {errno.EPERM, errno.ENOTSUP, getattr(errno, "EOPNOTSUPP", errno.ENOTSUP), errno.EMLINK}
    RD = os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK | os.O_NOCTTY | os.O_CLOEXEC
    path = Path(path)
    target = Path(os.path.realpath(path)) if path.is_symlink() else path
    parent, name = target.parent, target.name
    tmp_name = f".{name}.{os.getpid()}.{uuid.uuid4().hex[:12]}.tmp"
    st = {"step": "open_directory", "dfd": None, "fd": None, "tmp_id": None, "created": False,
          "mode": None, "use_existing": False}
    data = text.encode("utf-8")

    def ident(x):
        return (x.st_dev, x.st_ino)

    def close_fd():
        if st["fd"] is not None:
            fd, st["fd"] = st["fd"], None
            try:
                os.close(fd)
            except OSError:
                pass

    def cleanup():
        where = str(parent / tmp_name)
        if not st["created"]:
            close_fd()
            return {"status": "not_created", "created": False}
        if st["tmp_id"] is None and st["fd"] is not None:
            try:                                         # identify it through our own descriptor
                st["tmp_id"] = ident(os.fstat(st["fd"]))
            except OSError:
                pass
        close_fd()
        if st["tmp_id"] is None:
            return {"status": "left_unidentified", "created": True, "path": where,
                    "note": "created by this run but could not be identified, so it was not removed"}
        try:
            cur = os.stat(tmp_name, dir_fd=st["dfd"], follow_symlinks=False)
        except FileNotFoundError:
            return {"status": "already_absent", "created": True, "path": where}
        except OSError as e:
            return {"status": "failed", "created": True, "path": where, "error": f"{type(e).__name__}: {e}"}
        if ident(cur) != st["tmp_id"] or not _stat.S_ISREG(cur.st_mode):
            return {"status": "left_untouched_not_this_runs_file", "created": True, "path": where}
        try:
            os.unlink(tmp_name, dir_fd=st["dfd"])
        except OSError as e:
            return {"status": "failed", "created": True, "path": where, "error": f"{type(e).__name__}: {e}"}
        return {"status": "removed", "created": True, "path": where}

    def close_dir():
        if st["dfd"] is not None:
            d, st["dfd"] = st["dfd"], None
            os.close(d)

    def open_regular(nm):
        """Open nm in the pinned directory only if it is a regular file: lstat first (a FIFO / device / symlink is
        rejected without being opened), then a non-blocking, non-following open, then fstat again."""
        pre = os.stat(nm, dir_fd=st["dfd"], follow_symlinks=False)
        if not _stat.S_ISREG(pre.st_mode):
            raise OSError(errno.EINVAL, f"not a regular file (mode {oct(pre.st_mode)})", str(parent / nm))
        fd = os.open(nm, RD, dir_fd=st["dfd"])
        try:
            fs = os.fstat(fd)
            if not _stat.S_ISREG(fs.st_mode):
                raise OSError(errno.EINVAL, "not a regular file", str(parent / nm))
        except BaseException:
            os.close(fd)
            raise
        return fd, fs

    def same_bytes(fd):
        chunks, size = [], 0
        while size <= len(data):
            b = os.read(fd, 1 << 16)
            if not b:
                break
            chunks.append(b)
            size += len(b)
        return b"".join(chunks) == data

    # ---------------- before publishing: nothing visible changes under the cache name -------------------------
    try:
        st["dfd"] = os.open(str(parent), os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC)
        dfd = st["dfd"]
        st["step"] = "inspect_existing"
        old, old_x = None, {}
        try:
            ofd, old = open_regular(name)               # attributes come from the file actually opened
        except FileNotFoundError:
            ofd = None
        if ofd is not None:
            try:
                old_x = {k: os.getxattr(ofd, k) for k in os.listxattr(ofd)}
            finally:
                os.close(ofd)
            if old.st_uid != os.getuid():
                raise PermissionError(errno.EPERM, "cache file owned by another user; its owner cannot be kept "
                                      "by an atomic replacement", str(target))
            if not os.access(name, os.W_OK, dir_fd=dfd):
                raise PermissionError(errno.EACCES, "Permission denied", str(target))
        st["step"] = "create_temp"
        st["fd"] = os.open(tmp_name, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW | os.O_CLOEXEC,
                           0o600 if old is not None else 0o666, dir_fd=dfd)
        st["created"] = True
        fd = st["fd"]
        st["tmp_id"] = ident(os.fstat(fd))
        _cache_write_checkpoint("temp_created", dfd, tmp_name)
        if old is not None:                              # attributes on the fd, BEFORE any content
            st["step"] = "set_attributes"
            if os.fstat(fd).st_gid != old.st_gid:
                os.fchown(fd, -1, old.st_gid)
            for k in set(os.listxattr(fd)) - set(old_x):  # e.g. inherited default ACL
                os.removexattr(fd, k)
            for k, v in old_x.items():
                try:
                    same = os.getxattr(fd, k) == v
                except OSError:
                    same = False
                if not same:
                    os.setxattr(fd, k, v)
            if "system.posix_acl_access" not in old_x:
                os.fchmod(fd, _stat.S_IMODE(old.st_mode))
            st["step"] = "verify_attributes"
            now = os.fstat(fd)
            got = {k: os.getxattr(fd, k) for k in os.listxattr(fd)}
            bad = sorted(set(got) ^ set(old_x)) + sorted(k for k in set(got) & set(old_x) if got[k] != old_x[k])
            if bad or now.st_gid != old.st_gid or _stat.S_IMODE(now.st_mode) != _stat.S_IMODE(old.st_mode):
                raise OSError(errno.EPERM, "attributes of the existing cache could not be reproduced "
                              f"(xattrs {bad}, gid {now.st_gid}/{old.st_gid}, mode "
                              f"{oct(_stat.S_IMODE(now.st_mode))}/{oct(_stat.S_IMODE(old.st_mode))})", str(target))
        _cache_write_checkpoint("attributes_set", dfd, tmp_name)
        st["step"] = "write"
        view = memoryview(data); off = 0
        while off < len(data):
            off += os.write(fd, view[off:])
        st["step"] = "fsync_file"
        os.fsync(fd)
        st["step"] = "close_temp"
        st["fd"] = None
        os.close(fd)
        _cache_write_checkpoint("before_replace", dfd, tmp_name)
        st["step"] = "verify_temp_name"
        cur = os.stat(tmp_name, dir_fd=dfd, follow_symlinks=False)
        if ident(cur) != st["tmp_id"] or not _stat.S_ISREG(cur.st_mode):
            raise OSError(errno.ESTALE, "the temporary name no longer refers to this run's file "
                          "(the cache directory was modified concurrently)", str(parent / tmp_name))
        st["step"] = "check_readable"                    # never publish a cache this user cannot read back
        try:
            rfd = os.open(tmp_name, RD, dir_fd=dfd)
        except PermissionError as e:
            raise PermissionError(errno.EACCES, "the new cache file would not be readable by this user "
                                  "(umask / default ACL / copied mode); it was not published",
                                  str(parent / tmp_name)) from e
        try:
            if ident(os.fstat(rfd)) != st["tmp_id"]:
                raise OSError(errno.ESTALE, "the temporary name no longer refers to this run's file",
                              str(parent / tmp_name))
        finally:
            os.close(rfd)
        if old is None:                                  # new cache: publish without clobbering
            def adopt_or_refuse():                       # another writer published first
                st["step"] = "verify_concurrent_existing"
                efd, _ = open_regular(name)
                try:
                    identical = same_bytes(efd)
                finally:
                    os.close(efd)
                if not identical:
                    raise OSError(errno.EEXIST, "a cache with different content appeared concurrently; it was "
                                  "left in place", str(target))
                st["use_existing"] = True
            st["step"] = "publish_new"
            try:
                os.link(tmp_name, name, src_dir_fd=dfd, dst_dir_fd=dfd, follow_symlinks=False)
                st["mode"] = "link"
            except FileExistsError:
                adopt_or_refuse()
            except OSError as e:
                if e.errno not in NOLINK:
                    raise
                st["step"] = "publish_new_noreplace"     # no hard links: rename WITHOUT replacing, or fail
                try:
                    _cache_rename_noreplace(tmp_name, name, dfd)
                    st["mode"] = "rename_noreplace"
                except FileExistsError:
                    adopt_or_refuse()
                except OSError as e2:
                    raise OSError(e2.errno or errno.ENOTSUP,
                                  "no non-clobbering way to publish a new cache on this file system "
                                  f"(link: {e}; renameat2 RENAME_NOREPLACE: {e2}); nothing was published",
                                  str(target)) from e2
        else:
            st["step"] = "replace"
            os.rename(tmp_name, name, src_dir_fd=dfd, dst_dir_fd=dfd)
            st["mode"] = "rename"
    except BaseException as exc:
        step = st["step"]
        c = cleanup()
        try:
            close_dir()
        except OSError:
            pass
        if isinstance(exc, OSError):
            raise _CacheWriteError("before_replace", step, target, exc, False, c) from exc
        if hasattr(exc, "add_note"):
            exc.add_note(f"cache temp cleanup: {c}")
        raise
    if st["use_existing"]:                               # identical cache already there: adopt it
        c = {"status": "present", "created": True, "path": str(parent / tmp_name)}
        try:
            st["step"] = "remove_temp"
            c = cleanup()
            if c.get("status") != "removed":
                raise OSError(errno.EIO, "an identical cache was already published, but this run's temporary "
                              f"file could not be removed ({c.get('status')})", str(parent / tmp_name))
            st["step"] = "fsync_directory"               # make the adopted name durable from our side too
            warn = None
            try:
                os.fsync(st["dfd"])
            except OSError as e:
                if e.errno not in UNSUP:
                    raise
                warn = {"reason": "directory_fsync_unsupported", "errno": e.errno,
                        "error": f"{type(e).__name__}: {e}", "directory": str(parent),
                        "replaced": False, "adopted_existing": True, "durability": "not_guaranteed"}
            st["step"] = "close_directory"
            close_dir()
            return warn
        except BaseException as exc:
            step = st["step"]
            try:
                close_dir()
            except OSError:
                pass
            if isinstance(exc, OSError):
                raise _CacheWriteError("adopt_existing", step, target, exc, False, c) from exc
            if hasattr(exc, "add_note"):
                exc.add_note(f"cache write interrupted while adopting an identical cache (step {step}); the "
                             f"cache name was not changed by this run; temp cleanup: {c}")
            raise
    # ---------------- after publishing: the cache name holds the new file ------------------------------------
    where = str(parent / tmp_name)

    def drop_temp_link():                                # remove only the temp NAME; never the cache name
        try:
            cur = os.stat(tmp_name, dir_fd=st["dfd"], follow_symlinks=False)
        except FileNotFoundError:
            return {"status": "already_absent", "created": True, "path": where}, None
        except OSError as e:
            return {"status": "failed", "created": True, "path": where, "error": f"{type(e).__name__}: {e}"}, e
        if ident(cur) != st["tmp_id"]:
            return {"status": "left_untouched_not_this_runs_file", "created": True, "path": where}, None
        try:
            os.unlink(tmp_name, dir_fd=st["dfd"])
        except OSError as e:
            return {"status": "failed", "created": True, "path": where, "error": f"{type(e).__name__}: {e}"}, e
        return {"status": "removed", "created": True, "path": where}, None

    if st["mode"] == "link":                             # the cache and the temp name share one inode
        done = {"status": "present", "created": True, "path": where}
    else:
        done = {"status": "not_applicable", "created": True}
    try:
        st["step"] = "after_rename"
        _cache_write_checkpoint("after_replace", st["dfd"], tmp_name)
        if st["mode"] == "link":                         # drop the temporary second name
            st["step"] = "remove_temp_name"
            done, err = drop_temp_link()
            if err is not None:
                raise err
        st["step"] = "verify_target"
        vfd, cur = open_regular(name)                    # never blocks on a FIFO / device
        try:
            if ident(cur) != st["tmp_id"] and not same_bytes(vfd):
                raise OSError(errno.ESTALE, "the cache name refers to a different file with different content "
                              "(the cache directory was modified concurrently)", str(target))
        finally:
            os.close(vfd)
        st["step"] = "fsync_directory"
        warn = None
        try:
            os.fsync(st["dfd"])
        except OSError as e:
            if e.errno not in UNSUP:
                raise
            warn = {"reason": "directory_fsync_unsupported", "errno": e.errno,
                    "error": f"{type(e).__name__}: {e}", "directory": str(parent),
                    "replaced": True, "durability": "not_guaranteed"}
        st["step"] = "close_directory"
        close_dir()
        return warn
    except BaseException as exc:
        step = st["step"]
        if done.get("status") == "present" and st["dfd"] is not None:
            done, _ = drop_temp_link()                   # cache name untouched
        try:
            close_dir()
        except OSError:
            pass
        if isinstance(exc, OSError):
            raise _CacheWriteError("after_replace", step, target, exc, True, done) from exc
        if hasattr(exc, "add_note"):
            exc.add_note(f"cache write interrupted after publishing (step {step}); the cache name already holds "
                         f"the new file; directory descriptor closed; temp cleanup: {done}")
        raise


def initialize_thermochemistry(c: Config, default_cache_dir: Path) -> ThermochemistryBackend:
    """Build (or load from cache) the chamber-property backend and its density."""
    require_finite_config(c)   # finiteness only; the range rules stay in validate()
    cache_root = Path(c.cea_cache_dir).expanduser() if c.cea_cache_dir else Path(default_cache_dir)
    backend = build_backend(c, c.reactants_spec_list, c.packing_fraction, cache_root)
    _set_thermo(backend)
    # Additive-volume mixture density x packing fraction.  A CEA combustion-product
    # (gas) density is never used as the solid-propellant density.
    c.rho_p = backend.bulk_mixture_density
    return backend


def thermochemistry() -> ThermochemistryBackend:
    backend = _thermo_state()
    if backend is None:
        raise RuntimeError("Thermochemistry has not been initialized.")
    return backend


def gas_props(p0: float, c: Config) -> Tuple[float, float, float]:
    """Return GAS-PHASE chamber Rg [J/(kg K)], equilibrium T0 [K], gamma_s [-]."""
    return thermochemistry().props(p0)


def psi_and_deriv(p0: float, c: Config) -> Tuple[float, float]:
    """Return Psi = R_eff*T0 [J/kg] and its analytic dPsi/dp [J/(kg Pa)].

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


# ---- formerly in post-processing (report 113: diagnostic definitions used by the CEA backend) ----
# Cleanup P3a: current diagnostic definitions (single source of truth)
CEA_DENSITY_NOTE = ("EqSolution.density = p*M_CEA/(Ru*T) "
                    "[kg-total/m^3]: total mixture mass per gas-phase "
                    "volume (condensed-phase volume neglected); equals "
                    "the gas-phase ideal density p/(R_gas*T) divided "
                    "by Yg. It equals the two-phase mixture density "
                    "only under the negligible condensed-volume "
                    "assumption. Diagnostic only")
N_GAS_PER_KG_TOTAL_DEFINITION = (
    "EqSolution.n: total gas-phase moles per kg of TOTAL equilibrium "
    "products [kmol-gas/kg-total]; M_CEA = 1/n")


def _sonic_column_status(column: Any, n_grid: int) -> Tuple[str, str]:
    """Classify a stored sonic-velocity column without trusting it.

    UNAVAILABLE : missing, empty, or every entry NaN/null (no value obtained)
    AVAILABLE   : list of the grid length, every entry a finite positive number
    INVALID     : anything else (wrong type/length, non-numeric, not
                  convertible to float, non-positive, infinite, or a mix of
                  values and missing entries).  Never raises.
    """
    if column is None:
        return "UNAVAILABLE", "column missing"
    if not isinstance(column, list):
        return "INVALID", f"column is {type(column).__name__}, not a list"
    if len(column) == 0:
        return "UNAVAILABLE", "column empty"
    vals: List[Optional[float]] = []
    for v in column:
        if v is None:
            vals.append(None)
        elif isinstance(v, bool) or not isinstance(v, (int, float)):
            return "INVALID", f"non-numeric entry {_report_value_text(v, 80)}"
        else:
            try:
                vals.append(float(v))
            except Exception as exc:  # noqa: BLE001 - untrusted cache data; classification must not raise
                return "INVALID", (f"entry not convertible to float "
                                   f"({type(exc).__name__}): {_report_value_text(v, 80)}")
    missing = [v is None or math.isnan(v) for v in vals]
    if all(missing):
        return "UNAVAILABLE", "all entries NaN/null"
    if len(vals) != n_grid:
        return "INVALID", f"length {len(vals)} != pressure grid {n_grid}"
    if any(missing):
        return "INVALID", "mix of values and NaN/null entries"
    if any((not math.isfinite(v)) or v <= 0.0 for v in vals):
        return "INVALID", "non-finite or non-positive entry"
    return "AVAILABLE", "all entries finite and positive"