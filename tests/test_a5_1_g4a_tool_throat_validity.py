"""G4a evidence tool (tools/a5_1_g4a_cea_comparison.py): CEA throat validity and aggregation.

The overall c* check must never pass on the remaining points when some CEA reference
points are invalid, and must never pass when no valid reference point exists.
Evidence tooling only; G4a stays not passed at 2e-3 (owner decision 2026-11-02).
"""
import importlib.util, math, sys
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("a5_1_g4a_cea_comparison", ROOT / "tools" / "a5_1_g4a_cea_comparison.py")
tool = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(tool)

P0 = 30e5


def _valid_throat(**over):
    kw = dict(converged=True, mach=1.0, area_ratio=1.0, cstar_cea=1500.0, p0_pa=P0, rho_t=2.0, u_t=P0 / (2.0 * 1500.0))
    kw.update(over)
    return tool.throat_validity(**kw)


def _row(name, pc, rel, valid=True):
    th = _valid_throat() if valid else _valid_throat(mach=0.90)
    return {"case": name, "pc_bar": pc, "cstar_rel_diff": rel, "cea_throat": th,
            "phase_sum_err": 0.0, "sonic_identity_max_rel_err_stations": 0.0, "density_identity_max_rel_err_stations": 0.0}


@pytest.mark.fast
def test_valid_throat_is_valid():
    th = _valid_throat()
    assert th["valid"] and th["invalid_reasons"] == [] and abs(th["own_throat_cstar_rel"]) < 1e-12


@pytest.mark.fast
@pytest.mark.parametrize("over, reason", [
    (dict(converged=False), "not_converged"),
    (dict(mach=0.90), "throat_mach_not_sonic"),          # observed: Mach_t ~0.90 with converged=True
    (dict(area_ratio=1.01), "throat_area_ratio_not_one"),
    (dict(cstar_cea=1500.0 * 1.01), "cstar_inconsistent_with_own_throat"),
    (dict(mach=float("nan")), "non_finite_station_value"),
])
def test_invalid_throat_conditions(over, reason):
    th = _valid_throat(**over)
    assert not th["valid"] and reason in th["invalid_reasons"]


@pytest.mark.fast
def test_converged_flag_alone_is_not_enough():
    th = _valid_throat(converged=True, mach=0.8978)
    assert th["converged"] and not th["valid"]


@pytest.mark.fast
def test_all_valid_and_within_passes():
    agg = tool.aggregate_cstar([_row("a", 10, 1e-3), _row("a", 30, -1.5e-3)])
    assert agg["passed"] and agg["status"] == "all_within"


@pytest.mark.fast
def test_valid_point_exceeding_fails():
    agg = tool.aggregate_cstar([_row("a", 10, 1e-3), _row("a", 30, 2.5e-3)])
    assert not agg["passed"] and agg["status"] == "exceeds"


@pytest.mark.fast
def test_invalid_reference_blocks_pass_even_if_remaining_points_are_within():
    rows = [_row("a", 10, 1e-3), _row("a", 30, 1e-3), _row("al8", 10, -8e-2, valid=False)]
    agg = tool.aggregate_cstar(rows)
    assert agg["n_within_among_valid"] == agg["n_valid_reference"] == 2
    assert not agg["passed"] and agg["status"] == "reference_invalid_points_present"
    assert agg["reference_invalid_points"] == ["al8@10"]
    res = tool.build_result(rows)
    assert res["checks_within_tolerance"]["cstar"] is False
    assert res["checks_within_tolerance"]["cstar_reference_valid"] is False


@pytest.mark.fast
def test_no_valid_reference_never_passes():
    rows = [_row("al8", 10, 1e-4, valid=False), _row("al8", 30, 1e-4, valid=False)]
    agg = tool.aggregate_cstar(rows)
    assert not agg["passed"] and agg["status"] == "no_valid_reference" and agg["max_abs_rel_diff_valid"] is None
    assert tool.build_result(rows)["checks_within_tolerance"]["cstar"] is False


@pytest.mark.fast
def test_empty_row_set_never_passes_any_check():
    agg = tool.aggregate_cstar([])
    assert not agg["passed"] and agg["status"] == "no_points"
    checks = tool.build_result([])["checks_within_tolerance"]
    assert not any(checks.values())


def _load_tool_fresh(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "tools" / "a5_1_g4a_cea_comparison.py")
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod


@pytest.mark.fast
def test_tool_records_the_body_it_executed():
    import hashlib
    body = sorted(ROOT.glob("GRIBS_v*.py"))[0].resolve()
    info = tool.GRIBS_MODULE_INFO
    assert Path(info["path"]) == body
    assert info["sha256"] == hashlib.sha256(body.read_bytes()).hexdigest()
    assert info["module_name"] != "gribs" and info["reused_existing_module"] is False
    assert Path(tool.gribs.__file__).resolve() == body
    assert tool.build_result([])["gribs_module"] == info


@pytest.mark.fast
def test_foreign_gribs_module_is_neither_reused_nor_replaced(tmp_path):
    import types
    fake = types.ModuleType("gribs"); fake.__file__ = str(tmp_path / "other_tree" / "GRIBS_v0.0.0.py")
    saved = sys.modules.get("gribs")
    sys.modules["gribs"] = fake
    try:
        fresh = _load_tool_fresh("a5_1_g4a_tool_foreign_check")
        assert fresh.gribs is not fake
        assert Path(fresh.GRIBS_MODULE_INFO["path"]) == sorted(ROOT.glob("GRIBS_v*.py"))[0].resolve()
        assert sys.modules["gribs"] is fake          # not replaced
    finally:
        if saved is None:
            sys.modules.pop("gribs", None)
        else:
            sys.modules["gribs"] = saved


def _cea_throat(reactants, pc):
    import cea, numpy as np
    names = [r[0] for r in reactants]
    w = np.array([r[1] for r in reactants], float); t = np.full(len(names), 298.15)
    rmix = cea.Mixture(names); pmix = cea.Mixture(names, products_from_reactants=True, omit=[], ions=False)
    hR = float(rmix.calc_property(cea.ENTHALPY, w, t)) / float(cea.R)
    rs = cea.RocketSolver(pmix, reactants=rmix, transport=False, trace=1e-10); sol = cea.RocketSolution(rs)
    rs.solve(sol, w, pc, supar=[4.0], hc=hR, iac=True); k = int(sol.num_pts)
    A = lambda n: np.asarray(getattr(sol, n), float)[:k]
    Ma, a, P, rho, ae = A("Mach"), A("sonic_velocity"), A("P"), A("density"), A("ae_at")
    it = int(np.argmin(np.abs(Ma - 1.0)))
    return tool.throat_validity(bool(sol.converged), float(Ma[it]), float(ae[it]), float(A("c_star")[0]),
                                float(P[0]) * 1e5, float(rho[it]), float(Ma[it] * a[it]))


@pytest.mark.integration
def test_real_cea_throat_failure_is_detected():
    # AP70/Al8/HTPB22 at 10 bar: RocketSolution reports converged with Mach_t ~0.90 (Al2O3 near melting at the throat)
    al8 = [("NH4CLO4(I)", 70.0), ("AL(cr)", 8.0), ("C4H6,butadiene", 22.0)]
    bad = _cea_throat(al8, 10.0)
    assert bad["converged"] and not bad["valid"] and "throat_mach_not_sonic" in bad["invalid_reasons"]
    good = _cea_throat(al8, 70.0)
    assert good["valid"], good
