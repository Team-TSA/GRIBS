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


def _the_body():
    bodies = sorted(ROOT.glob("GRIBS_v*.py"))
    if len(bodies) != 1:
        raise RuntimeError(f"expected exactly one GRIBS_v*.py, found {[b.name for b in bodies]}")
    return bodies[0].resolve()


@pytest.mark.fast
def test_tool_records_the_body_it_executed():
    import hashlib
    body = _the_body()
    info = tool.GRIBS_MODULE_INFO
    assert Path(info["path"]) == body
    # independent recomputation of the composite identity (ADR_A6_PACKAGE_SPLIT, High-1)
    pkg_dir = (ROOT / "src" / "gribs").resolve()
    entries = sorted((p.relative_to(pkg_dir.parent).as_posix().encode(), p.read_bytes())
                     for p in pkg_dir.rglob("*.py")
                     if not any(x == "__pycache__" or x.endswith(".egg-info") or x.startswith(".")
                                for x in p.relative_to(pkg_dir.parent).parts))
    digest = hashlib.sha256(b"".join(k + b"\0" + str(len(v)).encode() + b"\0" + v + b"\0" for k, v in entries)).hexdigest()
    entry_bytes = body.read_bytes()
    assert info["entry_sha256"] == hashlib.sha256(entry_bytes).hexdigest()
    assert Path(info["package_root"]) == pkg_dir
    assert info["package_tree_digest"] == digest
    assert info["sha256"] == hashlib.sha256(entry_bytes + b"\0" + digest.encode("ascii")).hexdigest()
    assert info["sha256_kind"] == "body_composite_v1"
    assert info["module_name"] != "gribs" and info["reused_existing_module"] is False
    assert Path(tool.gribs.__file__).resolve() == body
    assert tool.build_result([])["gribs_module"] == info


@pytest.mark.fast
def test_package_tree_digest_and_body_sha256_known_answer(tmp_path):
    """Known answer, computed outside Python (printf | sha256sum, report 115 / run r488): the stream is
    gribs/Z.py\\0 2\\0 z\\n\\0  gribs/__init__.py\\0 2\\0 A\\n\\0  gribs/b.py\\0 0\\0 \\0  gribs/sub/c.py\\0 4\\0 c=1\\n\\0
    (UTF-8 byte order puts "Z" before "_" before "b"); excluded: __pycache__, *.egg-info, names starting with ".",
    files that are not *.py, and anything outside the package directory."""
    pkg = tmp_path / "gribs"
    files = {"Z.py": b"z\n", "__init__.py": b"A\n", "b.py": b"", "sub/c.py": b"c=1\n",
             "__pycache__/x.py": b"x", "e.egg-info/f.py": b"f", ".hidden/y.py": b"y", ".z.py": b"z", "data.txt": b"t"}
    for rel, data in files.items():
        (pkg / rel).parent.mkdir(parents=True, exist_ok=True)
        (pkg / rel).write_bytes(data)
    (tmp_path / "outside.py").write_bytes(b"o")
    tree = "73df6a0a96d077efef921e99e29b32894fc2324c6c9eb881042d0593e56b0152"
    assert tool.package_tree_digest(pkg) == tree
    assert tool.body_sha256(b"E\n", tree) == "637415a6f6f900b34d310cb4377534cbc5784a35e8da1081d6387f33c86d629a"


@pytest.mark.fast
def test_check_package_root_rejects_none():
    with pytest.raises(RuntimeError, match="did not import"):
        tool._check_package_root(None, ROOT / "src" / "gribs")


@pytest.mark.fast
def test_check_package_root_rejects_another_path(tmp_path):
    import types
    other = tmp_path / "other" / "src" / "gribs"; other.mkdir(parents=True); (other / "__init__.py").write_text("")
    fake = types.ModuleType("gribs"); fake.__file__ = str(other / "__init__.py"); fake.__path__ = [str(other)]
    with pytest.raises(RuntimeError, match="not from"):
        tool._check_package_root(fake, ROOT / "src" / "gribs")


@pytest.mark.fast
def test_check_package_root_rejects_namespace_package():
    import types
    ns = types.ModuleType("gribs"); ns.__file__ = None; ns.__path__ = [str(ROOT / "src" / "gribs"), "/elsewhere/gribs"]
    with pytest.raises(RuntimeError, match="not a regular package"):
        tool._check_package_root(ns, ROOT / "src" / "gribs")


_FOREIGN_CHECK = r"""
import importlib.util, sys, types
from pathlib import Path
root = Path(sys.argv[1])
bodies = sorted(root.glob("GRIBS_v*.py"))
if len(bodies) != 1:
    raise RuntimeError(f"expected exactly one GRIBS_v*.py, found {[b.name for b in bodies]}")
body = bodies[0].resolve()
fake = types.ModuleType("gribs_entry"); fake.__file__ = str(root / "other_tree" / "GRIBS_v0.0.0.py")
sys.modules["gribs_entry"] = fake                      # a foreign entry under the tests' registration name
spec = importlib.util.spec_from_file_location("a5_1_g4a_tool_foreign_check", root / "tools" / "a5_1_g4a_cea_comparison.py")
tool = importlib.util.module_from_spec(spec); spec.loader.exec_module(tool)
assert tool.gribs is not fake                          # not reused
assert Path(tool.GRIBS_MODULE_INFO["path"]) == body    # this tree's single entry
assert sys.modules["gribs_entry"] is fake              # not replaced
pkg = sys.modules.get("gribs")                         # the package, if the entry imported it, is this tree's src/gribs
assert pkg is None or Path(pkg.__file__).resolve().parent == (root / "src" / "gribs").resolve(), pkg
print("ok")
"""


@pytest.mark.integration
def test_foreign_gribs_module_is_neither_reused_nor_replaced():
    """Runs in a fresh interpreter: after the split the entry imports the package `gribs`, so a foreign module
    placed in this process would be hit by the entry's own imports (ADR_A6_PACKAGE_SPLIT, unit 0b)."""
    import subprocess
    r = subprocess.run([sys.executable, "-B", "-P", "-s", "-c", _FOREIGN_CHECK, str(ROOT)],
                       capture_output=True, text=True, timeout=600)
    assert r.returncode == 0 and r.stdout.strip().endswith("ok"), (r.returncode, r.stdout[-2000:], r.stderr[-4000:])


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
