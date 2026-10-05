#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Run the Phase A-3 numerical-convergence and event-policy sweep.

This harness adapts the historical investigation preserved under
Specification/ to the current GRIBS v0.5.0-alpha implementation.

It intentionally loads the current single-file GRIBS module dynamically.
Migration into gribs/numerics/sweep.py belongs to the later package
reorganization work.
"""

from __future__ import annotations

import argparse
import copy
import csv
import importlib.metadata
import importlib.util
import json
import math
import platform
import subprocess
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_GRIBS_PATH = REPOSITORY_ROOT / "GRIBS_v0.5.0-alpha.py"
DEFAULT_CONFIG_PATH = REPOSITORY_ROOT / "gribs_config.json"
DEFAULT_SWEEP_CONFIG_PATH = (
    REPOSITORY_ROOT / "tools" / "convergence_sweep_config.json"
)
DEFAULT_OUTPUT_DIR = (
    REPOSITORY_ROOT / "verification" / "a3-convergence" / "latest"
)
DEFAULT_CACHE_DIR = REPOSITORY_ROOT / "results" / "thermo_cache"

NUMERIC_METRICS = (
    "p_max_Pa",
    "burn_time_s",
    "F_max_N",
    "impulse_total_Ns",
    "total_mass_consistency_error",
    "propellant_mass_balance_error",
)

DIAGNOSTIC_METRICS = (
    "burn_stop_reason",
    "blowdown_stop_reason",
    "property_extrapolation",
    "regimes_visited",
)

CSV_PATH_NAME = "convergence_sweep.csv"
MARKDOWN_PATH_NAME = "convergence_sweep.md"
METADATA_PATH_NAME = "convergence_sweep_metadata.json"


@dataclass(frozen=True)
class SweepCase:
    """One convergence or event-policy sensitivity case."""

    setting: str
    display_value: str
    attribute: str
    value: Any
    comparison_kind: str
    is_reference: bool = False


def format_case_value(value: Any, format_spec: str) -> str:
    """Format a configured case value for display and output."""
    if format_spec == "s":
        if not isinstance(value, str):
            raise ValueError(
                "Display format 's' requires a string value."
            )
        return value

    if format_spec == "d":
        if isinstance(value, bool) or not isinstance(value, int):
            raise ValueError(
                "Display format 'd' requires an integer value."
            )
        return format(value, "d")

    try:
        return format(value, format_spec)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"Unable to format value {value!r} with "
            f"display_format {format_spec!r}."
        ) from exc


def load_sweep_configuration(path: Path) -> dict[str, Any]:
    """Load and validate the Phase A-3 sweep configuration."""
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        raise FileNotFoundError(
            f"Sweep configuration not found: {path}"
        ) from None
    except json.JSONDecodeError as exc:
        raise ValueError(
            f"Invalid sweep configuration JSON at {path}: {exc}"
        ) from exc

    if not isinstance(document, dict):
        raise ValueError(
            "Sweep configuration root must be a JSON object."
        )

    if document.get("schema_version") != "1.0":
        raise ValueError(
            "Sweep configuration schema_version must be '1.0'."
        )

    sweeps = document.get("sweeps")
    if not isinstance(sweeps, list) or not sweeps:
        raise ValueError(
            "Sweep configuration 'sweeps' must be a non-empty array."
        )

    seen_settings: set[str] = set()

    for index, sweep in enumerate(sweeps, start=1):
        location = f"sweeps[{index - 1}]"

        if not isinstance(sweep, dict):
            raise ValueError(f"{location} must be an object.")

        required = {
            "setting",
            "attribute",
            "comparison_kind",
            "values",
            "reference",
            "display_format",
        }
        missing = sorted(required - set(sweep))
        if missing:
            raise ValueError(
                f"{location} is missing required keys: "
                + ", ".join(missing)
            )

        setting = sweep["setting"]
        attribute = sweep["attribute"]
        comparison_kind = sweep["comparison_kind"]
        values = sweep["values"]
        reference = sweep["reference"]
        display_format = sweep["display_format"]

        if not isinstance(setting, str) or not setting:
            raise ValueError(
                f"{location}.setting must be a non-empty string."
            )

        if setting in seen_settings:
            raise ValueError(
                f"Duplicate sweep setting: {setting!r}."
            )
        seen_settings.add(setting)

        if not isinstance(attribute, str) or not attribute:
            raise ValueError(
                f"{location}.attribute must be a non-empty string."
            )

        if comparison_kind not in {
            "convergence",
            "event_sensitivity",
        }:
            raise ValueError(
                f"{location}.comparison_kind must be "
                "'convergence' or 'event_sensitivity'."
            )

        if not isinstance(values, list) or not values:
            raise ValueError(
                f"{location}.values must be a non-empty array."
            )

        if reference not in values:
            raise ValueError(
                f"{location}.reference must occur in values."
            )

        if not isinstance(display_format, str):
            raise ValueError(
                f"{location}.display_format must be a string."
            )

        display_values = [
            format_case_value(value, display_format)
            for value in values
        ]
        if len(display_values) != len(set(display_values)):
            raise ValueError(
                f"{location} has duplicate formatted values."
            )

    metrics = document.get("metrics")
    if not isinstance(metrics, dict):
        raise ValueError(
            "Sweep configuration 'metrics' must be an object."
        )

    relative = metrics.get("relative")
    absolute_limit = metrics.get("absolute_limit")
    diagnostic = metrics.get("diagnostic")

    if not isinstance(relative, dict) or not relative:
        raise ValueError(
            "metrics.relative must be a non-empty object."
        )

    if not isinstance(absolute_limit, dict) or not absolute_limit:
        raise ValueError(
            "metrics.absolute_limit must be a non-empty object."
        )

    if not isinstance(diagnostic, list) or not diagnostic:
        raise ValueError(
            "metrics.diagnostic must be a non-empty array."
        )

    numeric_metrics = set(NUMERIC_METRICS)
    configured_numeric_metrics = (
        set(relative) | set(absolute_limit)
    )
    if configured_numeric_metrics != numeric_metrics:
        raise ValueError(
            "Configured relative and absolute-limit metrics must "
            "exactly match the harness numeric metrics."
        )

    if set(diagnostic) != set(DIAGNOSTIC_METRICS):
        raise ValueError(
            "Configured diagnostic metrics must exactly match "
            "the harness diagnostic metrics."
        )

    for metric, specification in relative.items():
        if not isinstance(specification, dict):
            raise ValueError(
                f"metrics.relative.{metric} must be an object."
            )
        tolerance = specification.get("tolerance")
        unit = specification.get("unit")
        if (
            isinstance(tolerance, bool)
            or not isinstance(tolerance, (int, float))
            or tolerance <= 0.0
        ):
            raise ValueError(
                f"metrics.relative.{metric}.tolerance must be "
                "a positive number."
            )
        if not isinstance(unit, str) or not unit:
            raise ValueError(
                f"metrics.relative.{metric}.unit must be "
                "a non-empty string."
            )

    for metric, specification in absolute_limit.items():
        if not isinstance(specification, dict):
            raise ValueError(
                f"metrics.absolute_limit.{metric} must be an object."
            )
        maximum = specification.get("maximum")
        unit = specification.get("unit")
        if (
            isinstance(maximum, bool)
            or not isinstance(maximum, (int, float))
            or maximum <= 0.0
        ):
            raise ValueError(
                f"metrics.absolute_limit.{metric}.maximum must be "
                "a positive number."
            )
        if not isinstance(unit, str) or not unit:
            raise ValueError(
                f"metrics.absolute_limit.{metric}.unit must be "
                "a non-empty string."
            )

    acceptance = document.get("acceptance")
    if not isinstance(acceptance, dict):
        raise ValueError(
            "Sweep configuration 'acceptance' must be an object."
        )

    if not isinstance(acceptance.get("enforce"), bool):
        raise ValueError(
            "acceptance.enforce must be a boolean."
        )

    applied_kinds = acceptance.get(
        "apply_relative_tolerances_to"
    )
    if not isinstance(applied_kinds, list):
        raise ValueError(
            "acceptance.apply_relative_tolerances_to must be an array."
        )

    valid_kinds = {"convergence", "event_sensitivity"}
    if any(kind not in valid_kinds for kind in applied_kinds):
        raise ValueError(
            "acceptance.apply_relative_tolerances_to contains an "
            "unknown comparison kind."
        )

    return document


def build_cases(
    sweep_configuration: dict[str, Any] | None = None,
) -> list[SweepCase]:
    """Build cases from the configured Phase A-3 sweep matrix."""
    if sweep_configuration is None:
        sweep_configuration = load_sweep_configuration(
            DEFAULT_SWEEP_CONFIG_PATH
        )

    cases: list[SweepCase] = []

    for sweep in sweep_configuration["sweeps"]:
        reference = sweep["reference"]
        display_format = sweep["display_format"]

        for value in sweep["values"]:
            cases.append(
                SweepCase(
                    setting=sweep["setting"],
                    display_value=format_case_value(
                        value,
                        display_format,
                    ),
                    attribute=sweep["attribute"],
                    value=value,
                    comparison_kind=sweep["comparison_kind"],
                    is_reference=(value == reference),
                )
            )

    if not cases:
        raise ValueError(
            "Sweep configuration produced no cases."
        )

    return cases


def load_gribs(path: Path) -> Any:
    """Load the current single-file GRIBS implementation."""
    if not path.is_file():
        raise FileNotFoundError(f"GRIBS module not found: {path}")

    module_name = "gribs_a3_convergence_target"
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Unable to create import specification for {path}")

    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def run_git(
    arguments: Sequence[str],
    repository_root: Path,
) -> str | None:
    """Return stripped Git output, or None when unavailable."""
    try:
        completed = subprocess.run(
            ["git", *arguments],
            cwd=repository_root,
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None

    return completed.stdout.strip()


def package_version(distribution_name: str) -> str | None:
    """Return an installed distribution version when available."""
    try:
        return importlib.metadata.version(distribution_name)
    except importlib.metadata.PackageNotFoundError:
        return None


def cea_version() -> str | None:
    """Return the detected CEA Python package version."""
    for distribution_name in ("cea", "nasa-cea"):
        version = package_version(distribution_name)
        if version is not None:
            return version
    return None


def relative_difference(value: Any, reference: Any) -> float | None:
    """Return absolute relative difference, with safe zero handling."""
    try:
        numeric_value = float(value)
        numeric_reference = float(reference)
    except (TypeError, ValueError):
        return None

    if not (math.isfinite(numeric_value) and math.isfinite(numeric_reference)):
        return None

    if numeric_reference == 0.0:
        return 0.0 if numeric_value == 0.0 else math.inf

    return abs(numeric_value - numeric_reference) / abs(numeric_reference)


def normalize_diagnostic(value: Any) -> Any:
    """Convert diagnostic values to deterministic CSV-friendly values."""
    if isinstance(value, (list, tuple, set)):
        return json.dumps(
            sorted(str(item) for item in value),
            ensure_ascii=False,
        )
    if isinstance(value, dict):
        return json.dumps(
            value,
            ensure_ascii=False,
            sort_keys=True,
        )
    return value


def run_case(
    gribs: Any,
    base_config: Any,
    case: SweepCase,
    cache_dir: Path,
) -> dict[str, Any]:
    """Execute one case using a fresh configuration and backend."""
    config = copy.deepcopy(base_config)
    setattr(config, case.attribute, case.value)

    started = time.perf_counter()
    gribs.initialize_thermochemistry(config, cache_dir)
    result = gribs.run_model(config)
    history = gribs.sample(result, config)
    summary = gribs.summarize(result, history, config)
    wall_time_s = time.perf_counter() - started

    row: dict[str, Any] = {
        "setting": case.setting,
        "value": case.display_value,
        "comparison_kind": case.comparison_kind,
        "is_reference": case.is_reference,
        "status": "PASS",
        "wall_time_s": wall_time_s,
        "error_type": "",
        "error_message": "",
    }

    for metric in NUMERIC_METRICS:
        row[metric] = summary.get(metric)

    for metric in DIAGNOSTIC_METRICS:
        row[metric] = normalize_diagnostic(summary.get(metric))

    return row


def failure_row(case: SweepCase, exc: Exception) -> dict[str, Any]:
    """Create a result row for a failed case."""
    row: dict[str, Any] = {
        "setting": case.setting,
        "value": case.display_value,
        "comparison_kind": case.comparison_kind,
        "is_reference": case.is_reference,
        "status": "FAIL",
        "wall_time_s": None,
        "error_type": type(exc).__name__,
        "error_message": str(exc),
    }

    for metric in NUMERIC_METRICS:
        row[metric] = None

    for metric in DIAGNOSTIC_METRICS:
        row[metric] = None

    return row


def add_comparisons(rows: list[dict[str, Any]]) -> None:
    """Add relative differences and diagnostic equality to result rows."""
    references: dict[str, dict[str, Any]] = {}

    for row in rows:
        if row["is_reference"] and row["status"] == "PASS":
            references[row["setting"]] = row

    for row in rows:
        reference = references.get(row["setting"])
        row["reference_value"] = (
            reference["value"] if reference is not None else None
        )

        for metric in NUMERIC_METRICS:
            column = f"{metric}_relative_difference"
            if row["status"] != "PASS" or reference is None:
                row[column] = None
            else:
                row[column] = relative_difference(
                    row.get(metric),
                    reference.get(metric),
                )

        for metric in DIAGNOSTIC_METRICS:
            column = f"{metric}_matches_reference"
            if row["status"] != "PASS" or reference is None:
                row[column] = None
            else:
                row[column] = row.get(metric) == reference.get(metric)


def evaluate_acceptance(
    rows: list[dict[str, Any]],
    sweep_configuration: dict[str, Any],
) -> dict[str, Any]:
    """Evaluate configured tolerances without changing computed results."""
    metrics = sweep_configuration["metrics"]
    relative_specifications = metrics["relative"]
    absolute_specifications = metrics["absolute_limit"]
    acceptance = sweep_configuration["acceptance"]

    applied_kinds = set(
        acceptance["apply_relative_tolerances_to"]
    )

    evaluated_rows = 0
    passed_rows = 0
    failed_rows = 0
    unevaluated_rows = 0
    violation_count = 0

    for row in rows:
        violations: list[str] = []
        checks_performed = 0

        if row["status"] != "PASS":
            row["acceptance_evaluated"] = True
            row["acceptance_pass"] = False
            row["acceptance_violations"] = "case_execution_failed"
            evaluated_rows += 1
            failed_rows += 1
            violation_count += 1

            for metric, specification in (
                relative_specifications.items()
            ):
                row[f"{metric}_tolerance"] = (
                    specification["tolerance"]
                )
                row[f"{metric}_within_tolerance"] = None

            for metric, specification in (
                absolute_specifications.items()
            ):
                row[f"{metric}_maximum"] = (
                    specification["maximum"]
                )
                row[f"{metric}_within_limit"] = None

            continue

        relative_is_applicable = (
            row["comparison_kind"] in applied_kinds
        )

        for metric, specification in (
            relative_specifications.items()
        ):
            tolerance = float(specification["tolerance"])
            difference = row.get(
                f"{metric}_relative_difference"
            )

            row[f"{metric}_tolerance"] = tolerance

            if not relative_is_applicable or difference is None:
                row[f"{metric}_within_tolerance"] = None
                continue

            checks_performed += 1
            within_tolerance = (
                math.isfinite(float(difference))
                and float(difference) <= tolerance
            )
            row[f"{metric}_within_tolerance"] = within_tolerance

            if not within_tolerance:
                violations.append(
                    f"{metric}: relative difference "
                    f"{float(difference):.16g} exceeds "
                    f"{tolerance:.16g}"
                )

        for metric, specification in (
            absolute_specifications.items()
        ):
            maximum = float(specification["maximum"])
            value = row.get(metric)

            row[f"{metric}_maximum"] = maximum

            try:
                numeric_value = float(value)
            except (TypeError, ValueError):
                within_limit = False
            else:
                within_limit = (
                    math.isfinite(numeric_value)
                    and abs(numeric_value) <= maximum
                )

            checks_performed += 1
            row[f"{metric}_within_limit"] = within_limit

            if not within_limit:
                violations.append(
                    f"{metric}: absolute value "
                    f"{value!r} exceeds {maximum:.16g}"
                )

        acceptance_evaluated = checks_performed > 0
        acceptance_pass = (
            not violations if acceptance_evaluated else None
        )

        row["acceptance_evaluated"] = acceptance_evaluated
        row["acceptance_pass"] = acceptance_pass
        row["acceptance_violations"] = "; ".join(violations)

        if not acceptance_evaluated:
            unevaluated_rows += 1
        elif acceptance_pass:
            evaluated_rows += 1
            passed_rows += 1
        else:
            evaluated_rows += 1
            failed_rows += 1
            violation_count += len(violations)

    return {
        "enforced": bool(acceptance["enforce"]),
        "evaluated_rows": evaluated_rows,
        "passed_rows": passed_rows,
        "failed_rows": failed_rows,
        "unevaluated_rows": unevaluated_rows,
        "violation_count": violation_count,
        "all_evaluated_rows_pass": failed_rows == 0,
    }


def csv_columns() -> list[str]:
    """Return the stable CSV column order."""
    columns = [
        "setting",
        "value",
        "comparison_kind",
        "is_reference",
        "reference_value",
        "status",
        "wall_time_s",
    ]

    for metric in NUMERIC_METRICS:
        columns.extend(
            [
                metric,
                f"{metric}_relative_difference",
            ]
        )

    columns.extend(
        [
            "acceptance_evaluated",
            "acceptance_pass",
            "acceptance_violations",
        ]
    )

    columns.extend(
        [
            "p_max_Pa_tolerance",
            "p_max_Pa_within_tolerance",
            "burn_time_s_tolerance",
            "burn_time_s_within_tolerance",
            "F_max_N_tolerance",
            "F_max_N_within_tolerance",
            "impulse_total_Ns_tolerance",
            "impulse_total_Ns_within_tolerance",
            "total_mass_consistency_error_maximum",
            "total_mass_consistency_error_within_limit",
            "propellant_mass_balance_error_maximum",
            "propellant_mass_balance_error_within_limit",
        ]
    )

    for metric in DIAGNOSTIC_METRICS:
        columns.extend(
            [
                metric,
                f"{metric}_matches_reference",
            ]
        )

    columns.extend(["error_type", "error_message"])
    return columns


def write_csv(rows: list[dict[str, Any]], path: Path) -> None:
    """Write complete raw and comparison results."""
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=csv_columns())
        writer.writeheader()
        writer.writerows(rows)


def format_number(value: Any, precision: int = 6) -> str:
    """Format a numeric value for compact Markdown output."""
    if value is None:
        return ""

    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return str(value)

    if math.isnan(numeric):
        return "nan"
    if math.isinf(numeric):
        return "inf"
    return f"{numeric:.{precision}g}"


def write_markdown(
    rows: list[dict[str, Any]],
    metadata: dict[str, Any],
    acceptance_summary: dict[str, Any],
    path: Path,
) -> None:
    """Write a human-readable convergence report."""
    lines = [
        "# Phase A-3 Numerical Convergence Sweep",
        "",
        "## Reproducibility",
        "",
        f"- UTC execution time: `{metadata['executed_at_utc']}`",
        f"- GRIBS version: `{metadata['gribs_version']}`",
        f"- Schema version: `{metadata['schema_version']}`",
        f"- Python: `{metadata['python_version']}`",
        f"- CEA package: `{metadata['cea_package_version']}`",
        f"- Git commit: `{metadata['git_commit']}`",
        f"- Git working tree clean before outputs: "
        f"`{metadata['git_working_tree_clean_before_outputs']}`",
        f"- Configuration: `{metadata['configuration_path']}`",
        f"- Case count: `{metadata['case_count']}`",
        "",
        "## Acceptance summary",
        "",
        f"- Enforcement enabled: "
        f"`{acceptance_summary['enforced']}`",
        f"- Evaluated rows: "
        f"`{acceptance_summary['evaluated_rows']}`",
        f"- Passed rows: "
        f"`{acceptance_summary['passed_rows']}`",
        f"- Failed rows: "
        f"`{acceptance_summary['failed_rows']}`",
        f"- Unevaluated rows: "
        f"`{acceptance_summary['unevaluated_rows']}`",
        f"- Total violations: "
        f"`{acceptance_summary['violation_count']}`",
        f"- All evaluated rows pass: "
        f"`{acceptance_summary['all_evaluated_rows_pass']}`",
        "",
        "Acceptance failures are reported even when enforcement is "
        "disabled. With enforcement disabled, tolerance violations do "
        "not change the process exit status.",
        "",
        "## Results",
        "",
        "| setting | value | kind | reference | status | "
        "p_max [Pa] | rel. diff | burn time [s] | rel. diff | "
        "F_max [N] | rel. diff | impulse [N s] | rel. diff | "
        "mass consistency | rel. diff | propellant balance | rel. diff | wall [s] |",
        "|---|---:|---|---:|---|---:|---:|---:|---:|---:|---:|"
        "---:|---:|---:|---:|---:|---:|---:|",
    ]

    for row in rows:
        lines.append(
            "| "
            + " | ".join(
                [
                    str(row["setting"]),
                    str(row["value"]),
                    str(row["comparison_kind"]),
                    str(row.get("reference_value", "")),
                    str(row["status"]),
                    format_number(row.get("p_max_Pa")),
                    format_number(
                        row.get("p_max_Pa_relative_difference"),
                        precision=3,
                    ),
                    format_number(row.get("burn_time_s")),
                    format_number(
                        row.get("burn_time_s_relative_difference"),
                        precision=3,
                    ),
                    format_number(row.get("F_max_N")),
                    format_number(
                        row.get("F_max_N_relative_difference"),
                        precision=3,
                    ),
                    format_number(row.get("impulse_total_Ns")),
                    format_number(
                        row.get("impulse_total_Ns_relative_difference"),
                        precision=3,
                    ),
                    format_number(
                        row.get("total_mass_consistency_error"),
                        precision=3,
                    ),
                    format_number(
                        row.get(
                            "total_mass_consistency_error_relative_difference"
                        ),
                        precision=3,
                    ),
                    format_number(
                        row.get("propellant_mass_balance_error"),
                        precision=3,
                    ),
                    format_number(
                        row.get(
                            "propellant_mass_balance_error_relative_difference"
                        ),
                        precision=3,
                    ),
                    format_number(row.get("wall_time_s"), precision=3),
                ]
            )
            + " |"
        )

    lines.extend(
        [
            "",
            "## Diagnostic comparison",
            "",
            "| setting | value | burn stop | matches reference | "
            "blowdown stop | matches reference | extrapolation | "
            "matches reference | regimes | matches reference |",
            "|---|---:|---|---|---|---|---|---|---|---|",
        ]
    )

    for row in rows:
        lines.append(
            "| "
            + " | ".join(
                [
                    str(row["setting"]),
                    str(row["value"]),
                    str(row.get("burn_stop_reason", "")),
                    str(
                        row.get(
                            "burn_stop_reason_matches_reference",
                            "",
                        )
                    ),
                    str(row.get("blowdown_stop_reason", "")),
                    str(
                        row.get(
                            "blowdown_stop_reason_matches_reference",
                            "",
                        )
                    ),
                    str(row.get("property_extrapolation", "")),
                    str(
                        row.get(
                            "property_extrapolation_matches_reference",
                            "",
                        )
                    ),
                    str(row.get("regimes_visited", "")),
                    str(
                        row.get(
                            "regimes_visited_matches_reference",
                            "",
                        )
                    ),
                ]
            )
            + " |"
        )

    acceptance_failures = [
        row
        for row in rows
        if row.get("acceptance_pass") is False
    ]

    lines.extend(["", "## Acceptance violations", ""])

    if not acceptance_failures:
        lines.append("No acceptance violations.")
    else:
        for row in acceptance_failures:
            lines.append(
                f"- `{row['setting']}={row['value']}`: "
                f"{row['acceptance_violations']}"
            )

    failures = [row for row in rows if row["status"] != "PASS"]
    lines.extend(["", "## Case failures", ""])

    if not failures:
        lines.append("No case failures.")
    else:
        for row in failures:
            lines.append(
                f"- `{row['setting']}={row['value']}`: "
                f"`{row['error_type']}: {row['error_message']}`"
            )

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_metadata(
    gribs: Any,
    args: argparse.Namespace,
    case_count: int,
    working_tree_status_before_outputs: str | None,
) -> dict[str, Any]:
    """Build reproducibility metadata."""
    return {
        "executed_at_utc": datetime.now(timezone.utc).isoformat(),
        "program": getattr(gribs, "PROGRAM_NAME", "GRIBS"),
        "gribs_version": getattr(gribs, "PROGRAM_VERSION", "unknown"),
        "schema_version": getattr(gribs, "SCHEMA_VERSION", "unknown"),
        "python_version": platform.python_version(),
        "python_implementation": platform.python_implementation(),
        "platform": platform.platform(),
        "cea_package_version": cea_version(),
        "git_commit": run_git(
            ["rev-parse", "HEAD"],
            REPOSITORY_ROOT,
        ),
        "git_branch": run_git(
            ["branch", "--show-current"],
            REPOSITORY_ROOT,
        ),
        "git_working_tree_status_before_outputs": (
            working_tree_status_before_outputs
        ),
        "git_working_tree_clean_before_outputs": (
            working_tree_status_before_outputs == ""
        ),
        "command_line": [sys.executable, *sys.argv],
        "gribs_module_path": str(args.gribs.resolve()),
        "configuration_path": str(args.config.resolve()),
        "sweep_configuration_path": str(
            args.sweep_config.resolve()
        ),
        "sweep_configuration": args.sweep_configuration,
        "output_directory": str(args.output_dir.resolve()),
        "cache_directory": str(args.cache_dir.resolve()),
        "case_count": case_count,
        "selected_case_numbers": args.case_numbers,
        "numeric_metrics": list(NUMERIC_METRICS),
        "diagnostic_metrics": list(DIAGNOSTIC_METRICS),
    }


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Run the 20-case Phase A-3 GRIBS numerical-convergence "
            "and event-policy sweep."
        )
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG_PATH,
        help="GRIBS JSON configuration file.",
    )
    parser.add_argument(
        "--gribs",
        type=Path,
        default=DEFAULT_GRIBS_PATH,
        help="Current single-file GRIBS implementation.",
    )
    parser.add_argument(
        "--sweep-config",
        type=Path,
        default=DEFAULT_SWEEP_CONFIG_PATH,
        help=(
            "JSON file defining sweep values, metrics, and "
            "acceptance settings."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for CSV, Markdown, and metadata outputs.",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=DEFAULT_CACHE_DIR,
        help="Thermochemistry cache directory.",
    )
    parser.add_argument(
        "--list-cases",
        action="store_true",
        help="List the 20 cases without running GRIBS.",
    )
    parser.add_argument(
        "--case",
        type=int,
        action="append",
        dest="case_numbers",
        metavar="NUMBER",
        help=(
            "Run only the numbered case. May be specified multiple times. "
            "Case numbers are shown by --list-cases."
        ),
    )
    return parser


def list_cases(cases: list[SweepCase]) -> None:
    """Print the sweep matrix."""
    print(f"{'number':>6}  {'setting':<22} {'value':<10} {'kind'}")
    for number, case in enumerate(cases, start=1):
        reference = " [reference]" if case.is_reference else ""
        print(
            f"{number:6d}  {case.setting:<22} "
            f"{case.display_value:<10} "
            f"{case.comparison_kind}{reference}"
        )


def main(argv: Sequence[str] | None = None) -> int:
    """Run the complete sweep."""
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        sweep_configuration = load_sweep_configuration(
            args.sweep_config
        )
    except (OSError, ValueError) as exc:
        parser.error(str(exc))

    args.sweep_configuration = sweep_configuration
    cases = build_cases(sweep_configuration)

    if args.list_cases:
        list_cases(cases)
        return 0

    if args.case_numbers:
        invalid = sorted(
            {
                number
                for number in args.case_numbers
                if number < 1 or number > len(cases)
            }
        )
        if invalid:
            parser.error(
                "case numbers must be between 1 and "
                f"{len(cases)}; invalid: "
                + ", ".join(str(number) for number in invalid)
            )

        selected_numbers = set(args.case_numbers)
        cases = [
            case
            for number, case in enumerate(cases, start=1)
            if number in selected_numbers
        ]

    gribs = load_gribs(args.gribs)
    loaded = gribs.load_json_configuration(args.config)
    base_config = loaded.config

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.cache_dir.mkdir(parents=True, exist_ok=True)

    working_tree_status_before_outputs = run_git(
        ["status", "--short"],
        REPOSITORY_ROOT,
    )

    print(f"Running {len(cases)} Phase A-3 cases")
    print(f"Configuration: {args.config}")
    print(f"Output directory: {args.output_dir}")
    print()

    rows: list[dict[str, Any]] = []

    for number, case in enumerate(cases, start=1):
        label = f"{case.setting}={case.display_value}"
        print(f"[{number:02d}/{len(cases):02d}] {label}", flush=True)

        try:
            row = run_case(
                gribs,
                base_config,
                case,
                args.cache_dir,
            )
            print(
                f"  PASS in {row['wall_time_s']:.3f} s",
                flush=True,
            )
        except Exception as exc:
            row = failure_row(case, exc)
            print(
                f"  FAIL: {type(exc).__name__}: {exc}",
                file=sys.stderr,
                flush=True,
            )

        rows.append(row)

    add_comparisons(rows)
    acceptance_summary = evaluate_acceptance(
        rows,
        sweep_configuration,
    )

    metadata = build_metadata(
        gribs,
        args,
        len(cases),
        working_tree_status_before_outputs,
    )

    csv_path = args.output_dir / CSV_PATH_NAME
    markdown_path = args.output_dir / MARKDOWN_PATH_NAME
    metadata_path = args.output_dir / METADATA_PATH_NAME

    metadata["acceptance_summary"] = acceptance_summary

    write_csv(rows, csv_path)
    write_markdown(
        rows,
        metadata,
        acceptance_summary,
        markdown_path,
    )
    metadata_path.write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    print()
    print("Output files")
    print(f"  {csv_path}")
    print(f"  {markdown_path}")
    print(f"  {metadata_path}")

    failed = [row for row in rows if row["status"] != "PASS"]
    if failed:
        print(
            f"{len(failed)} of {len(rows)} cases failed.",
            file=sys.stderr,
        )
        return 1

    print(f"All {len(rows)} cases passed.")

    if acceptance_summary["failed_rows"]:
        print(
            f"Acceptance evaluation: "
            f"{acceptance_summary['failed_rows']} row(s) failed "
            f"with {acceptance_summary['violation_count']} "
            f"violation(s)."
        )
    else:
        print("Acceptance evaluation: all evaluated rows passed.")

    if (
        acceptance_summary["enforced"]
        and acceptance_summary["failed_rows"]
    ):
        print(
            "Acceptance enforcement is enabled; returning failure.",
            file=sys.stderr,
        )
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
