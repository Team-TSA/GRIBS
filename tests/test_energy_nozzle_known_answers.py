"""A-6 stage 2, C2: known answers of nozzle_state over the four flow regimes, and the re-export of _nozzle_flow.

tests/fixtures/nozzle_known_answers.json was generated from the UNREFACTORED nozzle_state (nozzle.py of c7bb092) with the
reference_constant backend of cases/golden/g1_reference_constant.json; this file compares every returned key with `==`
(floats are stored as JSON numbers, which round-trip exactly).  C2 only moves code, so nothing may differ.

K-01..K-04: one regime each (no-flow, subsonic, choked, separated).
S-01..S-04 (C8): the flow from a chamber state (gribs.energy.nozzle_from_state).
K-05: static check that the entry imports _nozzle_flow from gribs.nozzle and that nozzle.py defines it.  That the entry
and the package hold the same object at run time is guaranteed by tests/test_a6_reexports.py; K-05 is a static aid.
"""
import ast
import dataclasses
import inspect
import json
import math
from pathlib import Path

import pytest

import gribs.config.loader
import gribs.nozzle
import gribs.thermo
from gribs.energy import balance
from gribs.energy.nozzle_from_state import nozzle_state_from_chamber

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


def _state_of_the_table(c, p0):
    """The chamber state that the property table gives at p0: (T0, R, gamma) from gas_props, Yg from phase_fractions."""
    R, T0, g = gribs.thermo.gas_props(p0, c)
    Yg, Yc = gribs.thermo.phase_fractions(p0, c)
    return R, T0, g, Yg, Yc


@pytest.mark.fast
def test_s01_the_mass_flow_is_proportional_to_p_over_sqrt_r_eff_t(tmp_path):
    c = _backend_config(tmp_path)
    At = math.pi * c.R_t0 ** 2
    R, T0, g, Yg, Yc = _state_of_the_table(c, 5.0e6)

    def mdot(p, T, R_eff):
        return nozzle_state_from_chamber(p, T, R_eff, g, Yg, At, At, c)["mdot_total"]       # eps = 1: choked at these pressures
    base = mdot(5.0e6, T0, R)
    expected = c.Cd * At * 5.0e6 * math.sqrt(g / (R * T0)) * (2.0 / (g + 1.0)) ** ((g + 1.0) / (2.0 * (g - 1.0)))
    assert base == pytest.approx(expected, rel=1e-12)
    assert mdot(10.0e6, T0, R) == pytest.approx(2.0 * base, rel=1e-12)                      # proportional to p
    assert mdot(5.0e6, 4.0 * T0, R) == pytest.approx(0.5 * base, rel=1e-12)                 # ~ 1 / sqrt(T)
    assert mdot(5.0e6, T0, 4.0 * R) == pytest.approx(0.5 * base, rel=1e-12)                 # ~ 1 / sqrt(R_eff)


@pytest.mark.fast
def test_s02_the_state_of_the_table_gives_what_the_legacy_nozzle_state_gives(tmp_path):
    c = _backend_config(tmp_path)
    assert len(FIXTURE["cases"]) == 12
    for case in FIXTURE["cases"]:                                                          # no-flow, subsonic, choked, separated
        R, T0, g, Yg, Yc = _state_of_the_table(c, case["p0"])
        new = nozzle_state_from_chamber(case["p0"], T0, R, g, Yg, case["At"], case["Ae"], c, Yc=Yc)
        old = gribs.nozzle.nozzle_state(case["p0"], case["At"], case["Ae"], c)
        assert new == old and new["regime"] == case["regime"], (case["regime"], case["p0"])


@pytest.mark.fast
def test_s03_the_pressure_argument_is_the_stagnation_pressure(tmp_path):
    assert "stagnation" in nozzle_state_from_chamber.__doc__
    c = _backend_config(tmp_path)
    At = math.pi * c.R_t0 ** 2
    R, T0, g, Yg, _ = _state_of_the_table(c, 5.0e6)
    p_static = 5.0e6
    p_stag = p_static * (1.0 + 0.5 * (g - 1.0) * 0.3 ** 2) ** (g / (g - 1.0))               # the stagnation pressure at M = 0.3
    flow_stag = nozzle_state_from_chamber(p_stag, T0, R, g, Yg, At, At, c)["mdot_total"]
    flow_static_taken_as_stag = nozzle_state_from_chamber(p_static, T0, R, g, Yg, At, At, c)["mdot_total"]
    assert flow_stag == pytest.approx(flow_static_taken_as_stag * p_stag / p_static, rel=1e-12)   # passing a static p would under-state the flow
    assert p_stag > p_static


@pytest.mark.fast
def test_s04_cd_scales_the_mass_flow_only(tmp_path):
    c = _backend_config(tmp_path)
    c8 = dataclasses.replace(c, Cd=0.8 * c.Cd)
    At = math.pi * c.R_t0 ** 2
    R, T0, g, Yg, _ = _state_of_the_table(c, 5.0e6)
    a = nozzle_state_from_chamber(5.0e6, T0, R, g, Yg, At, At, c)
    b = nozzle_state_from_chamber(5.0e6, T0, R, g, Yg, At, At, c8)
    assert b["mdot_total"] == pytest.approx(0.8 * a["mdot_total"], rel=1e-14)
    assert b["ve"] == a["ve"] and b["Me"] == a["Me"] and b["pe"] == a["pe"]                   # the exit state does not depend on Cd
    assert list(inspect.signature(balance.h0_out_from_chamber).parameters) == ["h_static", "v_inlet", "slip"]   # H0_out has no Cd
