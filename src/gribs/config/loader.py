"""JSON configuration (A-6 unit 9). Moved verbatim from section 12 of the entry: the v0.4+ schema readers and
checks, the schema migrations (v0.3.0, v0.4.3, v0.4.4), the Configuration record, configuration_from_document
and load_json_configuration.

Tests patch load_json_configuration here ("gribs.config.loader.load_json_configuration"); gribs.cli reads it
as gribs.config.loader.load_json_configuration and the entry re-exports every name, so a patch on the entry
would not reach it. One line differs from the entry text: a relative path of the configuration is resolved
against _PROGRAM_DIR (gribs.thermo: the directory that holds src/, i.e. the entry's directory) instead of the
directory of the entry file itself.
"""
from __future__ import annotations

import copy
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from gribs.config import Config, _is_nonfinite_number
from gribs.constants import (
    BACKEND_CEA_LEGACY_EXECUTABLE,
    BACKEND_CEA_PYTHON,
    BACKEND_REFERENCE_CONSTANT,
    KNOWN_BACKENDS,
    KNOWN_TWO_PHASE_MODES,
    LEGACY_BACKEND_ALIASES,
    LEGACY_SCHEMA_VERSION,
    PRE_TWO_PHASE_SCHEMA_VERSION,
    PRE_V050_SCHEMA_VERSION,
    REMOVED_BACKENDS,
    SCHEMA_VERSION,
    TWO_PHASE_CONDENSED_VOLUME_CHOICES,
    TWO_PHASE_HOMOGENEOUS,
    TWO_PHASE_NOZZLE_ENTRAINMENT_CHOICES,
    TWO_PHASE_SINGLE_LEGACY,
)
from gribs.errors import ConfigurationError
from gribs.thermo import _PROGRAM_DIR, ReactantSpec, ideal_mixture_density


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
              "solver", "output", "energy_model")
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
    return p if p.is_absolute() else _PROGRAM_DIR / p


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
                             exclusive_minimum=0.0, maximum=1.0)),
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
    # A-6 stage 2: the optional energy_model section; its settings are an attribute, not a Config field, so that
    # dataclasses.asdict(c) (resolved_configuration) is the same for every legacy configuration.
    if "energy_model" in root:
        from gribs.energy.config_parse import parse_energy_model
        c.energy = parse_energy_model(root["energy_model"])  # type: ignore[attr-defined]

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
