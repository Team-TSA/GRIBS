"""P5-5 (F-16): six run-state keys at the top level of the summary, and the run-state warnings.

The keys are solver_completed, mass_consistency_passed, nominal_termination_reached, model_domain_passed,
warnings_present and results_usable (= the first four together).  They are added only to a run that produced a
summary; a run that raised has no summary and gets no keys.  The warnings W_NON_NOMINAL_TERMINATION,
W_SUBSONIC_BURNING and W_PROPERTY_TABLE_EXTRAPOLATED are added by summarize().  run_passed, the exit code and all
computed values are unchanged.  Leaving the property-table range is reported by a warning and is not a model-domain
violation; At > Ae is.  The end-to-end fixtures are golden g1 (reference_constant, no CEA) and variants of it.
"""
import json
import math
from pathlib import Path

import numpy as np
import pytest

import gribs.api
from gribs.simulation import postprocess
from gribs.simulation.postprocess import (
    RUN_STATE_KEYS,
    W_NON_NOMINAL_TERMINATION,
    W_PROPERTY_TABLE_EXTRAPOLATED,
    W_SUBSONIC_BURNING,
    W_THROAT_EXCEEDS_EXIT_AREA,
    add_run_state_keys,
    append_run_state_warnings,
    nominal_termination,
    run_verification_exit_code,
)

ROOT = Path(__file__).resolve().parents[1]
G1 = ROOT / "cases" / "golden" / "g1_reference_constant.json"
NEW_CODES = {W_NON_NOMINAL_TERMINATION, W_SUBSONIC_BURNING, W_PROPERTY_TABLE_EXTRAPOLATED}


def _codes(summary):
    return [w.get("code") for w in summary["warnings"]]


def _summary(**over):
    s = {
        "pressure_min_Pa": 1.0e5, "pressure_max_Pa": 2.0e6,
        "burn_stop_reason": "burnout", "blowdown_stop_reason": "not-run",
        "post_run_consistency": {"passed": True}, "ledger": {"passed": True},
        "event_times": {"burnout": [0.5], "throat_exceeds_exit_s": []},
        "warnings": [],
    }
    s.update(over)
    return s


@pytest.mark.fast
def test_key_names_and_order():
    assert RUN_STATE_KEYS == ("solver_completed", "mass_consistency_passed", "nominal_termination_reached",
                              "model_domain_passed", "warnings_present", "results_usable")


@pytest.mark.fast
@pytest.mark.parametrize("burn, blow, expected", [
    ("burnout", "not-run", True),        # blowdown disabled
    ("burnout", "ambient", True),
    ("burnout", "time-limit", False),
    ("burnout", "choking-loss", False),
    ("ambient", "not-run", False),
    ("choking-loss", "not-run", False),
    ("time-limit", "not-run", False),
])
def test_nominal_termination(burn, blow, expected):
    assert nominal_termination(burn, blow) is expected


@pytest.mark.fast
def test_keys_of_a_nominal_run_are_true_bools():
    s = _summary()
    add_run_state_keys(s)
    assert {k: s[k] for k in RUN_STATE_KEYS} == {
        "solver_completed": True, "mass_consistency_passed": True, "nominal_termination_reached": True,
        "model_domain_passed": True, "warnings_present": False, "results_usable": True}
    assert all(type(s[k]) is bool for k in RUN_STATE_KEYS)


@pytest.mark.fast
@pytest.mark.parametrize("over, false_keys", [
    ({"burn_stop_reason": "time-limit"}, {"nominal_termination_reached"}),
    ({"blowdown_stop_reason": "time-limit"}, {"nominal_termination_reached"}),
    ({"event_times": {"throat_exceeds_exit_s": [0.1]}}, {"model_domain_passed"}),
    ({"warnings": [{"code": "E_OUT_OF_VALIDATED_RANGE"}]}, {"model_domain_passed"}),
    ({"ledger": {"passed": False}}, {"mass_consistency_passed"}),
    ({"post_run_consistency": {"passed": False}}, {"mass_consistency_passed"}),
    ({"pressure_max_Pa": float("nan")}, {"solver_completed"}),
    ({"pressure_min_Pa": float("inf")}, {"solver_completed"}),
])
def test_each_failed_state_makes_results_unusable(over, false_keys):
    s = _summary(**over)
    add_run_state_keys(s)
    assert {k for k in RUN_STATE_KEYS[:4] if s[k] is False} == false_keys
    assert s["results_usable"] is False


@pytest.mark.fast
def test_warnings_present_does_not_change_usability():
    s = _summary(warnings=[{"code": W_PROPERTY_TABLE_EXTRAPOLATED}, {"code": W_SUBSONIC_BURNING}])
    add_run_state_keys(s)
    assert s["warnings_present"] is True
    assert s["model_domain_passed"] is True and s["results_usable"] is True


def _history(phase, regime, extrap):
    n = len(phase)
    return {"t": np.linspace(0.0, 1.0, n), "phase": np.array(phase, dtype=object),
            "regime": np.array(regime, dtype=object), "extrap": np.array(extrap, dtype=bool)}


@pytest.mark.fast
def test_warnings_from_the_history():
    h = _history(["burn", "burn", "burn", "blowdown"], ["choked", "subsonic", "choked", "subsonic"],
                 [True, False, False, True])
    s = _summary(property_policy="clamp", property_range_min_Pa=2.0e5, property_range_max_Pa=1.0e7,
                 property_extrapolation=True, property_extrapolation_points=2, property_extrapolation_fraction=0.5)
    append_run_state_warnings(s, h)
    append_run_state_warnings(s, h)                    # idempotent
    assert _codes(s) == [W_SUBSONIC_BURNING, W_PROPERTY_TABLE_EXTRAPOLATED]
    sub = s["warnings"][0]["context"]
    assert sub["subsonic_burn_samples"] == 1 and sub["burn_samples"] == 3
    assert sub["first_subsonic_burn_time_s"] == pytest.approx(1.0 / 3.0)
    ext = s["warnings"][1]["context"]
    assert ext["policy"] == "clamp" and ext["points"] == 2 and ext["fraction"] == 0.5


@pytest.mark.fast
def test_no_warning_for_blowdown_only_subsonic_or_inside_table():
    h = _history(["burn", "burn", "blowdown"], ["choked", "choked", "subsonic"], [False, False, False])
    s = _summary(blowdown_stop_reason="ambient", property_extrapolation=False)
    append_run_state_warnings(s, h)
    assert s["warnings"] == []


@pytest.mark.fast
def test_non_nominal_warning():
    h = _history(["burn", "burn"], ["choked", "choked"], [False, False])
    s = _summary(burn_stop_reason="time-limit", property_extrapolation=False)
    append_run_state_warnings(s, h)
    assert _codes(s) == [W_NON_NOMINAL_TERMINATION]
    assert s["warnings"][0]["context"]["burn_stop_reason"] == "time-limit"


@pytest.mark.fast
def test_new_warnings_do_not_change_the_verdict_or_exit_code():
    s = _summary(warnings=[{"code": c} for c in sorted(NEW_CODES | {W_THROAT_EXCEEDS_EXIT_AREA})])
    assert postprocess.run_verdict_passed(s) is True
    assert run_verification_exit_code(s) == 0


def _run(tmp_path, name, edit=None):
    doc = json.loads(G1.read_text(encoding="utf-8"))
    if edit:
        edit(doc)
    doc["output"]["directory"] = str(tmp_path / f"out_{name}")
    path = tmp_path / f"{name}.json"
    path.write_text(json.dumps(doc), encoding="utf-8")
    return gribs.api.run_simulation(path)


def _keys(summary):
    return {k: summary[k] for k in RUN_STATE_KEYS}


def test_golden_g1_is_usable_with_a_table_range_warning(tmp_path):
    """g1 leaves the property-table pressure range; that is a warning, not a model-domain violation."""
    r = _run(tmp_path, "g1")
    s = r.summary
    assert s["property_extrapolation"] is True
    assert _keys(s) == {k: True for k in RUN_STATE_KEYS}
    assert _codes(s).count(W_PROPERTY_TABLE_EXTRAPOLATED) == 1
    assert W_NON_NOMINAL_TERMINATION not in _codes(s) and W_SUBSONIC_BURNING not in _codes(s)
    assert r.run_passed is True and run_verification_exit_code(s) == 0


def test_time_limit_is_not_nominal_and_keeps_exit_zero(tmp_path):
    def edit(doc):
        doc["solver"]["burning_time_limit_s"] = 0.1
    s = _run(tmp_path, "limit", edit).summary
    assert s["burn_stop_reason"] == "time-limit"
    assert s["nominal_termination_reached"] is False and s["results_usable"] is False
    assert s["solver_completed"] is True and s["mass_consistency_passed"] is True
    assert W_NON_NOMINAL_TERMINATION in _codes(s) and "W_INTEGRATION_TRUNCATED" in _codes(s)
    assert run_verification_exit_code(s) == 0


def test_throat_exceeding_exit_is_a_domain_violation(tmp_path):
    def edit(doc):
        doc["nozzle"]["expansion_ratio"] = 1.05
        doc["nozzle"]["throat_erosion"] = {"enabled": True, "rate_m_s_at_reference_pressure": 2.0e-5,
                                           "pressure_exponent": 0.0}
        doc["propellant"]["burn_law"]["coefficient_m_s"] = 0.001
    r = _run(tmp_path, "at_gt_ae", edit)
    s = r.summary
    assert s["event_times"]["throat_exceeds_exit_s"] and W_THROAT_EXCEEDS_EXIT_AREA in _codes(s)
    assert s["model_domain_passed"] is False and s["results_usable"] is False
    assert s["nominal_termination_reached"] is True and s["warnings_present"] is True
    assert r.run_passed is True and run_verification_exit_code(s) == 0


def test_subsonic_burning_is_a_warning_only(tmp_path):
    def edit(doc):
        doc["nozzle"]["initial_throat_radius_m"] = 0.012
        doc["environment"]["initial_chamber_pressure_Pa"] = 120000.0
        doc["propellant"]["burn_law"]["coefficient_m_s"] = 0.003
    s = _run(tmp_path, "subsonic", edit).summary
    assert "subsonic" in s["regimes_visited"] and s["burn_stop_reason"] == "burnout"
    assert W_SUBSONIC_BURNING in _codes(s)
    assert s["model_domain_passed"] is True and s["results_usable"] is True
    assert math.isfinite(s["p_max_Pa"])
