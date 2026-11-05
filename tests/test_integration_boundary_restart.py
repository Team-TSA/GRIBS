"""Solver restart at the property-table limits and the integration safety limits.

Fast tests drive _segmented_solve() with a synthetic pressure ODE (no thermochemistry).
Integration tests run the g1 / g2 golden configurations (CEA for g2) in tmp_path.
"""
import importlib.util
import json
import math
import random
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
BODY = sorted(ROOT.glob("GRIBS_v*.py"))[0]
GOLD = ROOT / "cases" / "golden"
spec = importlib.util.spec_from_file_location("gribs_boundary_restart", BODY)
g = importlib.util.module_from_spec(spec)
sys.modules["gribs_boundary_restart"] = g
spec.loader.exec_module(g)

PMIN, PMAX = 1.0e5, 8.0e6


def cfg(policy="clamp", backend="cea_python"):
    return SimpleNamespace(property_policy=policy, thermo_backend=backend, p_fit_min=PMIN, p_fit_max=PMAX)


def oscillating(t, y):
    # p(t) = 4e6 - 5e6 cos(2 pi t): starts at -1e6 (below p_fit_min), crosses p_fit_min and p_fit_max upward,
    # then both downward, twice per unit time
    return np.array([5.0e6 * 2.0 * math.pi * math.sin(2.0 * math.pi * t)])


def p_exact(t):
    return 4.0e6 - 5.0e6 * np.cos(2.0 * np.pi * np.asarray(t))


def solve(c, t1=1.0, events=(), **kw):
    budget = g._IntegrationBudget()
    kw.setdefault("method", "LSODA"); kw.setdefault("rtol", 1e-10); kw.setdefault("atol", np.array([1e-3]))
    kw.setdefault("max_step", 0.01); kw.setdefault("dense_output", True)
    return g._segmented_solve(oscillating, (0.0, t1), np.array([p_exact(0.0)]), c, budget, events=events, **kw), budget


@pytest.mark.fast
def test_restarts_at_upper_and_lower_limits_in_both_directions():
    sol, _ = solve(cfg())
    got = [(r["limit"], r["direction"]) for r in sol.table_restarts]
    assert got == [("p_fit_min", "up"), ("p_fit_max", "up"), ("p_fit_max", "down"), ("p_fit_min", "down")]
    for r in sol.table_restarts:
        level = PMAX if r["limit"] == "p_fit_max" else PMIN
        assert abs(r["p0_Pa"] - level) <= 1e-3 * level
        assert abs(p_exact(r["t_s"]) - level) <= 1e-6 * level
    assert len(sol.segments) == 5 and sol.success and sol.t[-1] == pytest.approx(1.0)


@pytest.mark.fast
def test_joined_solution_is_continuous_and_accurate():
    sol, budget = solve(cfg())
    tt = np.linspace(0.0, 1.0, 2001)
    assert np.max(np.abs(sol.sol(tt)[0] - p_exact(tt))) <= 1e-3 * PMAX
    for r in sol.table_restarts:                       # value at and just around each restart time
        t = r["t_s"]
        assert abs(sol.sol(t)[0] - p_exact(t)) <= 1e-6 * PMAX
        assert abs(sol.sol(t - 1e-9)[0] - sol.sol(t + 1e-9)[0]) <= 1e-3 * PMAX
    assert np.all(np.diff(sol.t) > 0.0)                 # no duplicated restart point
    assert sol.nfev == sum(int(s.nfev) for s in sol.segments) and budget.rhs_calls >= sol.nfev


@pytest.mark.fast
@pytest.mark.parametrize("c", [cfg(backend="reference_constant"), cfg(policy="extrapolate")])
def test_inactive_without_clamped_table_limits(c):
    sol, _ = solve(c)
    assert sol.table_restarts == [] and len(sol.segments) == 1


@pytest.mark.fast
def test_user_events_are_merged_and_terminal_events_stop():
    def p_at_6e6(t, y):
        return y[0] - 6.0e6
    p_at_6e6.terminal = False
    sol, _ = solve(cfg(), events=[p_at_6e6])
    assert len(sol.t_events[0]) == 2                    # once up, once down, in different segments
    def p_down_6e6(t, y):
        return y[0] - 6.0e6
    p_down_6e6.terminal = True; p_down_6e6.direction = -1.0
    sol2, _ = solve(cfg(), events=[p_down_6e6])
    assert len(sol2.t_events[0]) == 1 and sol2.t[-1] == pytest.approx(sol2.t_events[0][0])
    assert [r["limit"] for r in sol2.table_restarts] == ["p_fit_min", "p_fit_max", "p_fit_max"]


@pytest.mark.fast
def test_rhs_budget_raises(monkeypatch):
    monkeypatch.setattr(g, "INTEGRATION_MAX_RHS_CALLS", 50)
    with pytest.raises(g.IntegrationBudgetError, match="E_INTEGRATION_BUDGET_EXCEEDED: more than 50 right-hand-side"):
        solve(cfg())


@pytest.mark.fast
def test_wall_budget_raises(monkeypatch):
    monkeypatch.setattr(g, "INTEGRATION_MAX_WALL_S", 0.0)
    with pytest.raises(g.IntegrationBudgetError, match="E_INTEGRATION_BUDGET_EXCEEDED: wall time above 0 s"):
        solve(cfg(), t1=50.0, max_step=1e-3)          # > 1024 RHS calls before the first wall check


@pytest.mark.fast
def test_restart_budget_raises(monkeypatch):
    monkeypatch.setattr(g, "INTEGRATION_MAX_TABLE_RESTARTS", 2)
    with pytest.raises(g.IntegrationBudgetError, match="more than 2 solver restarts"):
        solve(cfg())


@pytest.mark.fast
def test_budget_error_is_runtime_error_for_cli_exit_4():
    assert issubclass(g.IntegrationBudgetError, RuntimeError)
    assert g.IntegrationBudgetError.code == "E_INTEGRATION_BUDGET_EXCEEDED"
    assert (g.INTEGRATION_MAX_RHS_CALLS, g.INTEGRATION_MAX_WALL_S, g.INTEGRATION_MAX_TABLE_RESTARTS) == \
        (2_000_000, 300.0, 10_000)


# ----------------------------------------------------------------------------- integration (golden configurations)
def golden_doc(case, tmp_path, **changes):
    doc = json.loads(next(GOLD.glob(f"{case}_*.json")).read_text())
    doc["output"]["directory"] = str(tmp_path / "out")
    if "cea_python" in doc["thermochemistry"]:
        doc["thermochemistry"]["cea_python"]["cache_directory"] = str(tmp_path / "cache")
    for path, value in changes.items():
        node = doc
        keys = path.split("/")
        for k in keys[:-1]:
            node = node[k]
        node[keys[-1]] = value
    return doc


def run_doc(doc, tmp_path):
    c = g.configuration_from_document(doc).config
    g.initialize_thermochemistry(c, tmp_path / "cache")
    res = g.run_model(c)
    return c, res, g.sample(res, c)


@pytest.mark.integration
def test_g1_reference_constant_is_integrated_as_one_segment(tmp_path):
    _c, res, _h = run_doc(golden_doc("g1", tmp_path), tmp_path)
    assert res["integration_budget"]["table_restarts_burn"] == []
    assert len(res["sol_burn"].segments) == 1 and res["stop_burn"] == "burnout"


@pytest.mark.integration
def test_g2_restarts_once_at_the_upper_limit_and_completes(tmp_path):
    golden = json.loads((GOLD / "golden_values.json").read_text())["g2"]
    c, res, h = run_doc(golden_doc("g2", tmp_path), tmp_path)
    rs = res["integration_budget"]["table_restarts_burn"]
    assert [(r["limit"], r["direction"]) for r in rs] == [("p_fit_max", "up")]
    assert res["stop_burn"] == "burnout"
    assert res["t_events"]["burnout"][0] == pytest.approx(golden["burn_time_s"], rel=1e-9)
    assert float(np.max(h["p0"])) == pytest.approx(golden["p_max_Pa"], rel=1e-6)


@pytest.mark.integration
def test_g2_with_blowdown_restarts_downward_in_the_blowdown(tmp_path):
    _c, res, _h = run_doc(golden_doc("g2", tmp_path, **{"solver/blowdown/enabled": True}), tmp_path)
    down = [r for r in res["integration_budget"]["table_restarts_blowdown"] if r["limit"] == "p_fit_max"]
    assert down and down[0]["direction"] == "down"
    assert res["stop_blow"] in ("ambient", "time-limit")


@pytest.mark.integration
def test_g2_with_perturbed_gamma_table_completes(tmp_path):
    """The gamma_s table x (1 +- 1e-13) (seeded) stalled LSODA for > 600 s without the restart (r450-r452).

    Whether the old solver path stalls depends on the last digits of the CEA table of the environment; this test
    checks that the run completes in normal time with the restart.
    """
    doc = golden_doc("g2", tmp_path)
    doc["thermochemistry"]["cea_python"]["rocket_diagnostics"] = False
    cfg_build = tmp_path / "build.json"; cfg_build.write_text(json.dumps(doc))
    p = subprocess.run([sys.executable, str(BODY), "--config", str(cfg_build)], capture_output=True, text=True,
                       timeout=300)
    assert p.returncode == 0, p.stderr[-500:]
    cache = next((tmp_path / "cache").glob("*.json"))
    table = json.loads(cache.read_text())
    rng = random.Random("g2|G_eps1e-13_s2|2")
    table["gamma_s"] = [v * (1.0 + 1e-13 * rng.choice((-1.0, 1.0))) for v in table["gamma_s"]]
    cache.write_text(json.dumps(table))
    doc["thermochemistry"]["cea_python"]["rebuild_cache"] = False
    doc["output"]["directory"] = str(tmp_path / "out2")
    cfg_run = tmp_path / "run.json"; cfg_run.write_text(json.dumps(doc))
    p = subprocess.run([sys.executable, str(BODY), "--config", str(cfg_run)], capture_output=True, text=True,
                       timeout=120)
    assert p.returncode == 0, p.stderr[-500:]
    assert (tmp_path / "out2" / "summary.json").exists()
    line = [x for x in p.stderr.splitlines() if x.startswith("INTEGRATION_TABLE_RESTARTS ")]
    assert len(line) == 1
    rec = json.loads(line[0].split(" ", 1)[1])
    assert [(r["limit"], r["direction"]) for r in rec["table_restarts_burn"]] == [("p_fit_max", "up")]


@pytest.mark.integration
def test_cli_without_restart_writes_no_restart_line(tmp_path):
    doc = golden_doc("g1", tmp_path)
    cfg_run = tmp_path / "g1.json"; cfg_run.write_text(json.dumps(doc))
    p = subprocess.run([sys.executable, str(BODY), "--config", str(cfg_run)], capture_output=True, text=True, timeout=300)
    assert p.returncode == 0, p.stderr[-500:]
    assert "INTEGRATION_TABLE_RESTARTS" not in p.stderr
