"""CEA theoretical diagnostics must never break or falsify the run outputs.

The diagnostic is comparison-only.  When it fails, is incomplete or carries
invalid numbers, summary.json, summary.txt and the figure are still written,
the state is recorded explicitly, no CEA value is substituted, and the
simulation results and exit code are those of a run with a valid diagnostic.
"""
import copy
import importlib.util
import json
import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
MODULE_FILES = sorted(ROOT.glob("GRIBS_v*.py"))
assert len(MODULE_FILES) == 1, MODULE_FILES
spec = importlib.util.spec_from_file_location("gribs", MODULE_FILES[0])
gribs = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = gribs
spec.loader.exec_module(gribs)

G1_CONFIG = ROOT / "cases" / "golden" / "g1_reference_constant.json"

VALID = {
    "cea_theoretical_cstar_m_s": 1555.25,
    "cea_theoretical_cf": 1.23456,
    "cea_theoretical_isp_s": 190.1234,
    "cea_theoretical_isp_vacuum_s": 210.5678,
    "cea_theoretical_isp_exit_m_s": 1864.0,
    "cea_theoretical_chamber_pressure_Pa": 4.0e7,
    "cea_theoretical_chamber_temperature_K": 3000.0,
    "cea_theoretical_exit_definition": "throat (converging nozzle, Ae/At = 1)",
    "cea_theoretical_exit_mach": 1.0,
    "cea_theoretical_note": "comparison only",
}


# ---------------------------------------------------------------- unit level
@pytest.mark.fast
@pytest.mark.parametrize("diag,expected", [
    (None, ("not_requested", [], [])),
    (VALID, ("ok", [], [])),
    ({"cea_theoretical_error": "ThermochemistryError: boom"}, ("failed", [], [])),
    ({k: v for k, v in VALID.items() if k != "cea_theoretical_cf"},
     ("incomplete", ["cea_theoretical_cf"], [])),
    ({**VALID, "cea_theoretical_cstar_m_s": float("nan")},
     ("invalid", [], ["cea_theoretical_cstar_m_s"])),
    ({**VALID, "cea_theoretical_isp_s": float("inf")},
     ("invalid", [], ["cea_theoretical_isp_s"])),
    ({**VALID, "cea_theoretical_cf": True}, ("invalid", [], ["cea_theoretical_cf"])),
    ({**VALID, "cea_theoretical_cf": "1.2"}, ("invalid", [], ["cea_theoretical_cf"])),
    ({**VALID, "cea_theoretical_exit_definition": None},
     ("invalid", [], ["cea_theoretical_exit_definition"])),
    ({**VALID, "cea_theoretical_exit_mach": float("nan")},
     ("invalid", [], ["cea_theoretical_exit_mach"])),
    ({}, ("incomplete", list(gribs._CEA_DIAG_DISPLAY_NUMERIC_FIELDS
                             + gribs._CEA_DIAG_DISPLAY_TEXT_FIELDS), [])),
])
def test_cea_diagnostics_status_classification(diag, expected):
    status, missing, invalid = gribs.cea_diagnostics_status(diag)
    assert (status, missing, invalid) == expected


class _Backend:
    def __init__(self, behaviour):
        self.behaviour = behaviour

    def theoretical_rocket_performance(self, p, eps):
        if isinstance(self.behaviour, Exception):
            raise self.behaviour
        return copy.deepcopy(self.behaviour)


def _cfg():
    from types import SimpleNamespace
    return SimpleNamespace(cea_py_rocket_diagnostics=True, eps_nozzle=1.0)


@pytest.mark.fast
def test_producer_returns_success_unchanged():
    out = gribs._cea_theoretical_diagnostics(_cfg(), _Backend(VALID), {"p_mean_burn_Pa": 4.0e7})
    assert out == VALID
    assert "cea_theoretical_status" not in out


@pytest.mark.fast
def test_producer_records_failure_without_values(capsys):
    out = gribs._cea_theoretical_diagnostics(
        _cfg(), _Backend(gribs.ThermochemistryError("non-finite P")), {"p_mean_burn_Pa": 4.0e7})
    assert out["cea_theoretical_status"] == "failed"
    assert "non-finite P" in out["cea_theoretical_error"]
    assert not any(k.endswith(("_m_s", "_cf", "_s")) and k != "cea_theoretical_status" for k in out)
    assert "unavailable" in capsys.readouterr().out


@pytest.mark.fast
@pytest.mark.parametrize("bad,status,field,listkey", [
    ({**VALID, "cea_theoretical_cstar_m_s": float("nan")}, "invalid", "cea_theoretical_cstar_m_s",
     "cea_theoretical_invalid_fields"),
    ({k: v for k, v in VALID.items() if k != "cea_theoretical_isp_s"}, "incomplete",
     "cea_theoretical_isp_s", "cea_theoretical_missing_fields"),
])
def test_producer_annotates_invalid_and_incomplete(bad, status, field, listkey):
    out = gribs._cea_theoretical_diagnostics(_cfg(), _Backend(bad), {"p_mean_burn_Pa": 4.0e7})
    assert out["cea_theoretical_status"] == status
    assert field in out[listkey]
    # the received values are kept as they are (nothing substituted or dropped)
    for k, v in bad.items():
        assert (math.isnan(out[k]) and math.isnan(v)) if isinstance(v, float) and math.isnan(v) else out[k] == v


@pytest.mark.fast
def test_producer_not_requested():
    from types import SimpleNamespace
    assert gribs._cea_theoretical_diagnostics(
        SimpleNamespace(cea_py_rocket_diagnostics=False, eps_nozzle=1.0), _Backend(VALID), {}) is None
    assert gribs._cea_theoretical_diagnostics(_cfg(), object(), {}) is None


# ---------------------------------------------------------------- CLI level
SCENARIOS = {
    "ok": VALID,
    "failed": gribs.ThermochemistryError(
        "CEA returned a non-finite rocket-performance array 'P'"),
    "invalid": {**VALID, "cea_theoretical_cstar_m_s": float("nan")},
    "incomplete": {k: v for k, v in VALID.items() if k != "cea_theoretical_cf"},
}


def _run(tmp_path, monkeypatch, name):
    doc = json.loads(G1_CONFIG.read_text())
    out_dir = tmp_path / name
    doc["output"]["directory"] = str(out_dir)
    cfg_path = tmp_path / f"{name}.json"
    cfg_path.write_text(json.dumps(doc))
    behaviour = SCENARIOS[name]

    def fake(self, p, eps):
        if isinstance(behaviour, Exception):
            raise behaviour
        return copy.deepcopy(behaviour)

    monkeypatch.setattr(gribs.ReferenceConstantBackend, "theoretical_rocket_performance",
                        fake, raising=False)
    rc = gribs.main(["--config", str(cfg_path)])
    summary = json.loads((out_dir / "summary.json").read_text())
    text = (out_dir / "summary.txt").read_text()
    figures = list(out_dir.glob("*.png"))
    return rc, summary, text, figures


@pytest.mark.integration
def test_cli_outputs_for_all_diagnostic_states(tmp_path, monkeypatch, capsys):
    runs = {name: _run(tmp_path, monkeypatch, name) for name in SCENARIOS}
    rc_ok, sum_ok, text_ok, fig_ok = runs["ok"]
    assert fig_ok, "figure missing for the ok run"
    assert "cea_theoretical_status" not in sum_ok["cea_theoretical_diagnostics"]
    assert "cea_theoretical_cstar_m_s      : 1555.250" in text_ok
    assert "note: both sets are reported side by side" in text_ok

    for name in ("failed", "invalid", "incomplete"):
        rc, summary, text, figures = runs[name]
        # run verdict and simulation results are untouched by the diagnostic state
        assert rc == rc_ok, name
        assert summary["results"] == sum_ok["results"], name
        assert figures, f"figure missing for {name}"
        diag = summary["cea_theoretical_diagnostics"]
        assert diag["cea_theoretical_status"] == name
        assert f"status                         : {name}" in text
        assert "no CEA comparison value is substituted" in text
        assert "note: both sets are reported side by side" not in text
        assert "GRIBS effective c*" in text

    failed_diag = runs["failed"][1]["cea_theoretical_diagnostics"]
    assert "non-finite rocket-performance array" in failed_diag["cea_theoretical_error"]
    assert "cea_theoretical_cstar_m_s" not in failed_diag
    assert "cea_theoretical_cstar_m_s      :" not in runs["failed"][2]

    inv_text = runs["invalid"][2]
    assert "invalid fields                 : cea_theoretical_cstar_m_s" in inv_text
    assert "cea_theoretical_cstar_m_s      : invalid (nan)" in inv_text
    assert math.isnan(runs["invalid"][1]["cea_theoretical_diagnostics"]["cea_theoretical_cstar_m_s"])

    inc = runs["incomplete"]
    assert inc[1]["cea_theoretical_diagnostics"]["cea_theoretical_missing_fields"] == ["cea_theoretical_cf"]
    assert "cea_theoretical_cf             : unavailable (missing)" in inc[2]
