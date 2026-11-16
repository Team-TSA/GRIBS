"""P5-3 (F-06): a run whose eroding throat area exceeds the fixed exit area records the time and a warning.

Only the warning W_THROAT_EXCEEDS_EXIT_AREA and event_times["throat_exceeds_exit_s"] are added; the nozzle model,
the computed values and the run verdict are unchanged.  The fixtures are golden g1 (reference_constant, no CEA)
with throat erosion and the burn-rate coefficient lowered to 0.001 m/s so that the burn lasts past the crossing.
With pressure exponent 0 the throat radius grows linearly, Rt = Rt0 + c*t, so the crossing time is
t_x = Rt0*(sqrt(eps) - 1)/c.
"""
import json
import math
from pathlib import Path

import pytest

import gribs.api
import gribs.nozzle
from gribs.simulation.postprocess import W_THROAT_EXCEEDS_EXIT_AREA

ROOT = Path(__file__).resolve().parents[1]
G1 = ROOT / "cases" / "golden" / "g1_reference_constant.json"
EVENT = "throat_exceeds_exit_s"
BASE_EVENTS = {"burnout", "unchoked", "rechoked", "ambient"}


def _case(tmp_path, *, eps, rate):
    doc = json.loads(G1.read_text(encoding="utf-8"))
    doc["nozzle"]["expansion_ratio"] = eps
    doc["nozzle"]["throat_erosion"] = {"enabled": True, "rate_m_s_at_reference_pressure": rate, "pressure_exponent": 0.0}
    doc["propellant"]["burn_law"]["coefficient_m_s"] = 0.001
    doc["output"]["directory"] = str(tmp_path / "out")
    path = tmp_path / "case.json"
    path.write_text(json.dumps(doc), encoding="utf-8")
    return path, doc


def _codes(summary):
    return [w.get("code") for w in summary["warnings"]]


@pytest.mark.fast
def test_predicate_compares_radius_without_tolerance():
    R0 = 2.0e-3
    # The integrator's At (pi*Rt*Rt) is one ulp above eps*(pi*R0**2) for this radius with eps = 1.
    assert math.pi * R0 * R0 > 1.0 * (math.pi * R0 ** 2)
    assert not gribs.nozzle.throat_exceeds_exit_area(R0, R0, 1.0)
    assert gribs.nozzle.throat_exceeds_exit_area(math.nextafter(R0, 1.0), R0, 1.0)
    limit = gribs.nozzle.exit_area_throat_radius(R0, 1.05)
    assert not gribs.nozzle.throat_exceeds_exit_area(limit, R0, 1.05)
    assert gribs.nozzle.throat_exceeds_exit_area(math.nextafter(limit, 1.0), R0, 1.05)


def test_crossing_time_and_warning(tmp_path):
    path, doc = _case(tmp_path, eps=1.05, rate=2.0e-5)
    result = gribs.api.run_simulation(path)
    s = result.summary
    R0, rtol = doc["nozzle"]["initial_throat_radius_m"], doc["solver"]["relative_tolerance"]
    t_x = R0 * (math.sqrt(1.05) - 1.0) / 2.0e-5
    times = s["event_times"][EVENT]
    assert len(times) == 1
    assert abs(times[0] - t_x) <= 10.0 * rtol * t_x
    assert _codes(s).count(W_THROAT_EXCEEDS_EXIT_AREA) == 1
    w = next(w for w in s["warnings"] if w["code"] == W_THROAT_EXCEEDS_EXIT_AREA)
    assert w["context"]["time_s"] == times[0]
    assert w["context"]["phase"] == "burn"
    assert w["context"]["final_throat_to_exit_area_ratio"] > 1.0
    assert result.run_passed            # the warning does not change the run verdict


def test_eps_one_with_erosion_crosses_at_start(tmp_path):
    path, _ = _case(tmp_path, eps=1.0, rate=2.0e-5)
    s = gribs.api.run_simulation(path).summary
    assert s["event_times"][EVENT] == [0.0]
    assert _codes(s).count(W_THROAT_EXCEEDS_EXIT_AREA) == 1


def test_golden_g1_is_not_flagged(tmp_path):
    """eps = 1 without erosion: At equals Ae in the model; the one-ulp area difference must not be flagged."""
    doc = json.loads(G1.read_text(encoding="utf-8"))
    doc["output"]["directory"] = str(tmp_path / "out")
    path = tmp_path / "g1.json"
    path.write_text(json.dumps(doc), encoding="utf-8")
    s = gribs.api.run_simulation(path).summary
    assert s["event_times"][EVENT] == []
    assert set(s["event_times"]) == BASE_EVENTS | {EVENT}
    assert W_THROAT_EXCEEDS_EXIT_AREA not in _codes(s)
