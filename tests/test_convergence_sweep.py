"""Tests for the Phase A-3 convergence harness."""

from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import pytest


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
HARNESS_PATH = REPOSITORY_ROOT / "tools" / "convergence_sweep.py"


def load_harness():
    """Load the standalone A-3 harness as a test module."""
    module_name = "gribs_a3_convergence_harness_test"
    spec = importlib.util.spec_from_file_location(
        module_name,
        HARNESS_PATH,
    )

    assert spec is not None
    assert spec.loader is not None

    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def harness():
    return load_harness()


def test_build_cases_returns_expected_20_case_matrix(harness):
    cases = harness.build_cases()

    assert len(cases) == 20

    counts = {}
    for case in cases:
        counts[case.setting] = counts.get(case.setting, 0) + 1

    assert counts == {
        "rtol": 6,
        "atol_p": 4,
        "max_step_burn": 4,
        "cea_pressure_points": 4,
        "unchoked_policy": 2,
    }


def test_each_setting_has_exactly_one_reference_case(harness):
    cases = harness.build_cases()
    settings = {case.setting for case in cases}

    for setting in settings:
        references = [
            case
            for case in cases
            if case.setting == setting and case.is_reference
        ]
        assert len(references) == 1, setting


def test_reference_values_are_the_documented_values(harness):
    cases = harness.build_cases()

    references = {
        case.setting: case.display_value
        for case in cases
        if case.is_reference
    }

    assert references == {
        "rtol": "1e-11",
        "atol_p": "1e-04",
        "max_step_burn": "0.005",
        "cea_pressure_points": "161",
        "unchoked_policy": "switch",
    }


def test_comparison_kinds_separate_convergence_and_event_sensitivity(
    harness,
):
    cases = harness.build_cases()

    convergence_settings = {
        case.setting
        for case in cases
        if case.comparison_kind == "convergence"
    }
    event_settings = {
        case.setting
        for case in cases
        if case.comparison_kind == "event_sensitivity"
    }

    assert convergence_settings == {
        "rtol",
        "atol_p",
        "max_step_burn",
        "cea_pressure_points",
    }
    assert event_settings == {"unchoked_policy"}


@pytest.mark.parametrize(
    ("value", "reference", "expected"),
    [
        (10.0, 10.0, 0.0),
        (11.0, 10.0, 0.1),
        (9.0, 10.0, 0.1),
        (0.0, 0.0, 0.0),
        (1.0, 0.0, math.inf),
    ],
)
def test_relative_difference(harness, value, reference, expected):
    actual = harness.relative_difference(value, reference)

    if math.isinf(expected):
        assert math.isinf(actual)
    else:
        assert actual == pytest.approx(expected)


@pytest.mark.parametrize(
    ("value", "reference"),
    [
        (None, 1.0),
        (1.0, None),
        ("not-a-number", 1.0),
        (math.nan, 1.0),
        (1.0, math.inf),
    ],
)
def test_relative_difference_rejects_unusable_values(
    harness,
    value,
    reference,
):
    assert harness.relative_difference(value, reference) is None


def test_add_comparisons_uses_group_reference(harness):
    rows = [
        {
            "setting": "rtol",
            "value": "1e-06",
            "comparison_kind": "convergence",
            "is_reference": False,
            "status": "PASS",
            "p_max_Pa": 90.0,
            "burn_time_s": 10.0,
            "F_max_N": 18.0,
            "impulse_total_Ns": 40.0,
            "total_mass_consistency_error": 2e-6,
            "propellant_mass_balance_error": 4e-8,
            "burn_stop_reason": "burnout",
            "blowdown_stop_reason": "ambient",
            "property_extrapolation": False,
            "regimes_visited": '["choked", "subsonic"]',
        },
        {
            "setting": "rtol",
            "value": "1e-11",
            "comparison_kind": "convergence",
            "is_reference": True,
            "status": "PASS",
            "p_max_Pa": 100.0,
            "burn_time_s": 10.0,
            "F_max_N": 20.0,
            "impulse_total_Ns": 50.0,
            "total_mass_consistency_error": 1e-6,
            "propellant_mass_balance_error": 2e-8,
            "burn_stop_reason": "burnout",
            "blowdown_stop_reason": "ambient",
            "property_extrapolation": False,
            "regimes_visited": '["choked", "subsonic"]',
        },
    ]

    harness.add_comparisons(rows)

    candidate, reference = rows

    assert candidate["reference_value"] == "1e-11"
    assert reference["reference_value"] == "1e-11"

    assert candidate["p_max_Pa_relative_difference"] == pytest.approx(
        0.1
    )
    assert candidate[
        "impulse_total_Ns_relative_difference"
    ] == pytest.approx(0.2)

    assert reference["p_max_Pa_relative_difference"] == 0.0

    assert candidate[
        "burn_stop_reason_matches_reference"
    ] is True
    assert candidate[
        "blowdown_stop_reason_matches_reference"
    ] is True


def test_add_comparisons_leaves_difference_empty_without_reference(
    harness,
):
    rows = [
        {
            "setting": "rtol",
            "value": "1e-06",
            "comparison_kind": "convergence",
            "is_reference": False,
            "status": "PASS",
            "p_max_Pa": 90.0,
            "burn_time_s": 10.0,
            "F_max_N": 18.0,
            "impulse_total_Ns": 40.0,
            "total_mass_consistency_error": 2e-6,
            "propellant_mass_balance_error": 4e-8,
            "burn_stop_reason": "burnout",
            "blowdown_stop_reason": "ambient",
            "property_extrapolation": False,
            "regimes_visited": '["choked", "subsonic"]',
        }
    ]

    harness.add_comparisons(rows)

    assert rows[0]["reference_value"] is None
    assert rows[0]["p_max_Pa_relative_difference"] is None
    assert rows[0][
        "burn_stop_reason_matches_reference"
    ] is None


def test_failure_row_records_exception(harness):
    case = harness.build_cases()[0]

    row = harness.failure_row(
        case,
        RuntimeError("synthetic failure"),
    )

    assert row["status"] == "FAIL"
    assert row["error_type"] == "RuntimeError"
    assert row["error_message"] == "synthetic failure"
    assert row["p_max_Pa"] is None


def test_load_default_sweep_configuration(harness):
    document = harness.load_sweep_configuration(
        harness.DEFAULT_SWEEP_CONFIG_PATH
    )

    assert document["schema_version"] == "1.0"
    assert len(document["sweeps"]) == 5
    assert sum(
        len(sweep["values"])
        for sweep in document["sweeps"]
    ) == 20
    assert document["acceptance"]["enforce"] is False


def test_build_cases_uses_supplied_configuration(harness):
    document = {
        "schema_version": "1.0",
        "sweeps": [
            {
                "setting": "rtol",
                "attribute": "rtol",
                "comparison_kind": "convergence",
                "values": [1e-6, 1e-9],
                "reference": 1e-9,
                "display_format": ".0e",
            }
        ],
        "metrics": {
            "relative": {
                "p_max_Pa": {
                    "unit": "Pa",
                    "tolerance": 1e-6,
                },
                "burn_time_s": {
                    "unit": "s",
                    "tolerance": 1e-6,
                },
                "F_max_N": {
                    "unit": "N",
                    "tolerance": 1e-6,
                },
                "impulse_total_Ns": {
                    "unit": "N s",
                    "tolerance": 1e-6,
                },
            },
            "absolute_limit": {
                "total_mass_consistency_error": {
                    "unit": "1",
                    "maximum": 1e-4,
                },
                "propellant_mass_balance_error": {
                    "unit": "1",
                    "maximum": 1e-7,
                },
            },
            "diagnostic": list(
                harness.DIAGNOSTIC_METRICS
            ),
        },
        "acceptance": {
            "enforce": False,
            "apply_relative_tolerances_to": [
                "convergence"
            ],
        },
    }

    cases = harness.build_cases(document)

    assert len(cases) == 2
    assert cases[0].display_value == "1e-06"
    assert cases[0].is_reference is False
    assert cases[1].display_value == "1e-09"
    assert cases[1].is_reference is True


def test_load_sweep_configuration_rejects_missing_reference(
    harness,
    tmp_path,
):
    path = tmp_path / "invalid.json"
    path.write_text(
        """{
          "schema_version": "1.0",
          "sweeps": [
            {
              "setting": "rtol",
              "attribute": "rtol",
              "comparison_kind": "convergence",
              "values": [1e-6],
              "reference": 1e-9,
              "display_format": ".0e"
            }
          ]
        }
        """,
        encoding="utf-8",
    )

    with pytest.raises(
        ValueError,
        match="reference must occur in values",
    ):
        harness.load_sweep_configuration(path)


def acceptance_configuration(harness, enforce=False):
    return {
        "metrics": {
            "relative": {
                "p_max_Pa": {
                    "unit": "Pa",
                    "tolerance": 1e-6,
                },
                "burn_time_s": {
                    "unit": "s",
                    "tolerance": 1e-6,
                },
                "F_max_N": {
                    "unit": "N",
                    "tolerance": 1e-6,
                },
                "impulse_total_Ns": {
                    "unit": "N s",
                    "tolerance": 1e-6,
                },
            },
            "absolute_limit": {
                "total_mass_consistency_error": {
                    "unit": "1",
                    "maximum": 1e-4,
                },
                "propellant_mass_balance_error": {
                    "unit": "1",
                    "maximum": 1e-7,
                },
            },
            "diagnostic": list(
                harness.DIAGNOSTIC_METRICS
            ),
        },
        "acceptance": {
            "enforce": enforce,
            "apply_relative_tolerances_to": [
                "convergence"
            ],
        },
    }


def acceptance_row(
    *,
    comparison_kind="convergence",
    relative_difference=5e-7,
    mass_error=1e-5,
    propellant_error=1e-8,
):
    row = {
        "setting": "rtol",
        "value": "1e-06",
        "comparison_kind": comparison_kind,
        "is_reference": False,
        "status": "PASS",
        "total_mass_consistency_error": mass_error,
        "propellant_mass_balance_error": propellant_error,
    }

    for metric in (
        "p_max_Pa",
        "burn_time_s",
        "F_max_N",
        "impulse_total_Ns",
    ):
        row[metric] = 1.0
        row[f"{metric}_relative_difference"] = (
            relative_difference
        )

    return row


def test_evaluate_acceptance_passes_within_limits(harness):
    rows = [acceptance_row()]
    configuration = acceptance_configuration(harness)

    summary = harness.evaluate_acceptance(
        rows,
        configuration,
    )

    assert rows[0]["acceptance_evaluated"] is True
    assert rows[0]["acceptance_pass"] is True
    assert rows[0]["acceptance_violations"] == ""

    assert summary == {
        "enforced": False,
        "evaluated_rows": 1,
        "passed_rows": 1,
        "failed_rows": 0,
        "unevaluated_rows": 0,
        "violation_count": 0,
        "all_evaluated_rows_pass": True,
    }


def test_evaluate_acceptance_reports_relative_violation(harness):
    rows = [
        acceptance_row(
            relative_difference=2e-6,
        )
    ]
    configuration = acceptance_configuration(harness)

    summary = harness.evaluate_acceptance(
        rows,
        configuration,
    )

    assert rows[0]["acceptance_pass"] is False
    assert rows[0]["p_max_Pa_within_tolerance"] is False
    assert "p_max_Pa" in rows[0]["acceptance_violations"]

    assert summary["failed_rows"] == 1
    assert summary["violation_count"] == 4
    assert summary["all_evaluated_rows_pass"] is False


def test_evaluate_acceptance_uses_absolute_residual_limits(
    harness,
):
    rows = [
        acceptance_row(
            mass_error=2e-4,
            propellant_error=2e-7,
        )
    ]
    configuration = acceptance_configuration(harness)

    summary = harness.evaluate_acceptance(
        rows,
        configuration,
    )

    assert rows[0][
        "total_mass_consistency_error_within_limit"
    ] is False
    assert rows[0][
        "propellant_mass_balance_error_within_limit"
    ] is False

    assert summary["failed_rows"] == 1
    assert summary["violation_count"] == 2


def test_event_sensitivity_skips_relative_tolerances(harness):
    rows = [
        acceptance_row(
            comparison_kind="event_sensitivity",
            relative_difference=1.0,
        )
    ]
    configuration = acceptance_configuration(harness)

    summary = harness.evaluate_acceptance(
        rows,
        configuration,
    )

    assert rows[0]["p_max_Pa_within_tolerance"] is None
    assert rows[0]["burn_time_s_within_tolerance"] is None
    assert rows[0]["acceptance_pass"] is True
    assert summary["passed_rows"] == 1


def test_acceptance_summary_records_enforcement_flag(harness):
    rows = [acceptance_row()]
    configuration = acceptance_configuration(
        harness,
        enforce=True,
    )

    summary = harness.evaluate_acceptance(
        rows,
        configuration,
    )

    assert summary["enforced"] is True
