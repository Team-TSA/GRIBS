"""A-6 stage 2, C2: known answers of nozzle_state over the four flow regimes, and the re-export of _nozzle_flow.

tests/fixtures/nozzle_known_answers.json was generated from the UNREFACTORED nozzle_state (nozzle.py of c7bb092) with the
reference_constant backend of cases/golden/g1_reference_constant.json; this file compares every returned key with `==`
(floats are stored as JSON numbers, which round-trip exactly).  C2 only moves code, so nothing may differ.

K-01..K-04: one regime each (no-flow, subsonic, choked, separated).
K-05: static check that the entry imports _nozzle_flow from gribs.nozzle and that nozzle.py defines it.  That the entry
and the package hold the same object at run time is guaranteed by tests/test_a6_reexports.py; K-05 is a static aid.
"""
import ast
import json
import math
from pathlib import Path

import pytest

import gribs.config.loader
import gribs.nozzle
import gribs.thermo

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = json.loads((ROOT / "tests" / "fixtures" / "nozzle_known_answers.json").read_text(encoding="utf-8"))


def _backend_config(tmp_path):
    doc = json.loads((ROOT / "cases" / "golden" / "g1_reference_constant.json").read_text(encoding="utf-8"))
    c = gribs.config.loader.configuration_from_document(doc).config
    gribs.thermo.initialize_thermochemistry(c, tmp_path)   # the conftest guard restores the backend after the test
    return c


def _check_regime(regime, tmp_path):
    c = _backend_config(tmp_path)
    assert c.p_a == FIXTURE["p_a"]
    cases = [x for x in FIXTURE["cases"] if x["regime"] == regime]
    assert len(cases) == 3
    for case in cases:
        got = gribs.nozzle.nozzle_state(case["p0"], case["At"], case["Ae"], c)
        assert got["regime"] == regime
        assert set(got) == set(case["result"])
        for key, expected in case["result"].items():
            assert got[key] == expected, (regime, case["p0"], key, got[key], expected)
            assert not isinstance(expected, float) or math.isfinite(expected)


@pytest.mark.fast
def test_k01_no_flow(tmp_path):
    _check_regime("no-flow", tmp_path)


@pytest.mark.fast
def test_k02_subsonic(tmp_path):
    _check_regime("subsonic", tmp_path)


@pytest.mark.fast
def test_k03_choked(tmp_path):
    _check_regime("choked", tmp_path)


@pytest.mark.fast
def test_k04_separated(tmp_path):
    _check_regime("separated", tmp_path)


@pytest.mark.fast
def test_k05_entry_imports_nozzle_flow_from_gribs_nozzle_and_nozzle_defines_it():
    entries = sorted(ROOT.glob("GRIBS_v*.py"))
    assert len(entries) == 1
    sources = [node.module for node in ast.parse(entries[0].read_text(encoding="utf-8")).body
               if isinstance(node, ast.ImportFrom) and node.level == 0
               and "_nozzle_flow" in [a.name for a in node.names]]
    assert sources == ["gribs.nozzle"], sources      # imported, and from gribs.nozzle only
    nozzle = ast.parse((ROOT / "src" / "gribs" / "nozzle.py").read_text(encoding="utf-8"))
    assert [n.name for n in nozzle.body if isinstance(n, ast.FunctionDef) and n.name == "_nozzle_flow"] == ["_nozzle_flow"]
