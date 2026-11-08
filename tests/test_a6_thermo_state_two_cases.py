"""A-6 unit 3a (ADR_A6_PACKAGE_SPLIT, report 112 addendum): two cases run one after another in one process, without
the conftest state guard and without restoring the backend, give the same results as each case run alone."""
import json, subprocess, sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
_BODIES = sorted(ROOT.glob("GRIBS_v*.py"))
if len(_BODIES) != 1:
    raise RuntimeError(f"expected exactly one GRIBS_v*.py, found {[b.name for b in _BODIES]}")
BODY = _BODIES[0]
CHILD = Path(__file__).with_name("_a6_two_cases_child.py")


def _run(tmp_path, *cases):
    work = tmp_path / "-".join(cases)
    work.mkdir()
    r = subprocess.run([sys.executable, str(CHILD), str(BODY), str(work), *cases],
                       capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stderr
    records = [json.loads(line) for line in r.stdout.splitlines() if line.startswith("{")]
    assert [x["case"] for x in records] == list(cases), r.stdout
    return records


@pytest.mark.integration
def test_two_cases_in_one_process_without_restoring_match_each_case_alone(tmp_path):
    alone = {case: _run(tmp_path, case)[0] for case in ("A", "B")}
    assert alone["A"]["digest"] != alone["B"]["digest"]  # the cases differ, so a leaked backend would show
    for case, rec in alone.items():
        assert rec["state_before"] is None and rec["installed"], rec
        assert rec["stop_burn"] == "burnout" and rec["n_burn"] > 100, rec
    for order in (("A", "B"), ("B", "A")):
        first, second = _run(tmp_path, *order)
        assert first["state_before"] is None
        assert second["state_before"] == first["case"]  # not restored: the first case's backend was still installed
        for rec in (first, second):
            assert rec["installed"], rec
            ref = alone[rec["case"]]
            assert (rec["digest"], rec["rho_p"], rec["n_burn"]) == (ref["digest"], ref["rho_p"], ref["n_burn"]), (order, rec)
