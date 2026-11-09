"""A-6 stage 1 guard (ADR_A6_PACKAGE_SPLIT): names moved into the package are patched only on their owner.

MOVED grows, in the same commit as each move, by exactly the names that unit moves; ALLOWLIST (names not yet
moved, with their exact number of patch / assignment sites in tests) shrinks by the same names and is empty at
the end of stage 1. tests/red/ is a frozen historical record (A-5.1-3) and is not checked.
"""
import ast
import importlib
from collections import Counter
from pathlib import Path
from typing import NamedTuple, Optional

import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
TESTS = ROOT / "tests"

# name -> (owner module fully qualified, kind, exact number of patch sites on the owner in tests); kind is
# "function_or_global" or "class". T0 fails when a site on the owner disappears or appears (report 116).
MOVED = {
    "validate": ("gribs.config", "function_or_global", 2),  # unit 2: test_gribs.py (2)
    # unit 3b (36 sites before the move; the 4 restoring assignments of test_gribs.py became monkeypatch, so 32)
    "_cache_rename_noreplace": ("gribs.thermo", "function_or_global", 2),  # _cache_v7_child.py, test_cache_atomic_v7.py
    "_cache_write_checkpoint": ("gribs.thermo", "function_or_global", 15),  # the two children (2), test_cache_atomic_v7.py (13)
    "gas_props": ("gribs.thermo", "function_or_global", 1),
    "initialize_thermochemistry": ("gribs.thermo", "function_or_global", 6),  # unit 10: +1 test_a6_public_api.py
    "is_extrapolated": ("gribs.thermo", "function_or_global", 1),
    "phase_fractions": ("gribs.thermo", "function_or_global", 1),
    "psi_and_deriv": ("gribs.thermo", "function_or_global", 2),
    "thermochemistry": ("gribs.thermo", "function_or_global", 2),
    "PressureTableBackend": ("gribs.thermo", "class", 1),
    "ReferenceConstantBackend": ("gribs.thermo", "class", 1),
    "_ThermoSwap": ("gribs.thermo", "class", 1),
    # unit 4 (12 sites before the move; the 2 restoring assignments of test_gribs.py became monkeypatch, so 10, plus
    # the 2 sites of test_a6_cross_module_reads.py)
    "burn_rate": ("gribs.propellant.burning", "function_or_global", 2),  # test_gribs.py (1), cross-module test (1)
    "choke_limit_pressure": ("gribs.nozzle", "function_or_global", 2),
    "equilibrium_pressure": ("gribs.chamber.equilibrium", "function_or_global", 1),
    "nozzle_state": ("gribs.nozzle", "function_or_global", 5),  # test_gribs.py (4), cross-module test (1)
    "NozzleTwoPhaseModel": ("gribs.nozzle", "class", 2),
    # unit 5 (6 sites): test_integration_boundary_restart.py (the three limits), test_gribs.py (run_model)
    "INTEGRATION_MAX_RHS_CALLS": ("gribs.simulation.integrator", "function_or_global", 1),
    "INTEGRATION_MAX_TABLE_RESTARTS": ("gribs.simulation.integrator", "function_or_global", 1),
    "INTEGRATION_MAX_WALL_S": ("gribs.simulation.integrator", "function_or_global", 1),
    "run_model": ("gribs.simulation.integrator", "function_or_global", 4),  # unit 10: +1 test_a6_public_api.py
    # unit 6 (9 sites before the move; the assign / restore pair of the _phase_grid try/finally test of test_gribs.py
    # became one monkeypatch, so 8)
    "_phase_grid": ("gribs.simulation.postprocess", "function_or_global", 1),
    "mass_energy_ledger": ("gribs.simulation.postprocess", "function_or_global", 2),
    "model_applicability_error": ("gribs.simulation.postprocess", "function_or_global", 1),
    "nozzle_transition_diagnostics": ("gribs.simulation.postprocess", "function_or_global", 2),
    "sample": ("gribs.simulation.postprocess", "function_or_global", 1),
    "summarize": ("gribs.simulation.postprocess", "function_or_global", 1),
    # unit 7
    "make_figure": ("gribs.plotting", "function_or_global", 1),
    "summary_text": ("gribs.reporting", "function_or_global", 1),
    "write_csv": ("gribs.reporting", "function_or_global", 1),
    # unit 8
    "self_tests_post_run": ("gribs.selftest", "function_or_global", 1),
    # unit 9
    "load_json_configuration": ("gribs.config.loader", "function_or_global", 2),
    "_cea_theoretical_diagnostics": ("gribs.cli", "function_or_global", 1),
    "_print_backend_banner": ("gribs.cli", "function_or_global", 3),
    "environment_info": ("gribs.cli", "function_or_global", 1),
}

# names not moved yet -> exact number of patch / assignment sites in tests (AST census, report 110 / 111). Empty since
# unit 9: every patched name is in MOVED, so any patch of a GRIBS name on the entry fails T0.
ALLOWLIST = {}


# package module -> unit that introduced it (ADR_A6_PACKAGE_SPLIT; report 112 order). A module may import only
# gribs modules of the same or an earlier unit. Every module under src/gribs must be listed. The values are layers
# (report 114 addendum): gribs._version is layer 1 (gribs.constants imports it); the package root "gribs" is
# layer 10 since unit 10 (public API), so a module that imports a public name through `from gribs import X` fails.
MODULE_UNIT = {
    "gribs": 10,
    "gribs._version": 1,
    "gribs.constants": 1,
    "gribs.errors": 1,
    "gribs.util": 1,
    "gribs.config": 2,
    "gribs.thermo": 3,
    "gribs.propellant": 4,
    "gribs.propellant.burning": 4,
    "gribs.nozzle": 4,
    "gribs.chamber": 4,
    "gribs.chamber.equilibrium": 4,
    "gribs.simulation": 5,
    "gribs.simulation.integrator": 5,
    "gribs.simulation.postprocess": 6,
    "gribs.plotting": 7,
    "gribs.reporting": 7,
    "gribs.selftest": 8,
    "gribs.config.loader": 9,
    "gribs.cli": 9,
    "gribs.api": 10,
    "gribs.simulation.result": 10,
}


#: Module state that only its accessors may touch (ADR_A6_PACKAGE_SPLIT decision 2): name -> (reader, writer). Unit 3a
#: put the accessors in the entry; unit 3b moved the state, the accessors and every reader to gribs.thermo. T10.
ACCESSOR_ONLY = {"_THERMO": ("_thermo_state", "_set_thermo")}
#: The one module that defines the state and its accessors (unit 3b). T10 fails if the entry or another package module
#: defines them, so that there is never a second copy of the state.
ACCESSOR_OWNER = {"_THERMO": "gribs.thermo"}


def _entry():
    bodies = sorted(ROOT.glob("GRIBS_v*.py"))
    if len(bodies) != 1:
        raise RuntimeError(f"expected exactly one GRIBS_v*.py, found {[b.name for b in bodies]}")
    return bodies[0]


def _test_files():
    for path in sorted(TESTS.rglob("*.py")):
        rel = path.relative_to(TESTS)
        if rel.parts[0] == "red" or path.resolve() == Path(__file__).resolve():
            continue
        yield path


def _chain(node):
    parts = []
    while isinstance(node, ast.Attribute):
        parts.append(node.attr)
        node = node.value
    if isinstance(node, ast.Name):
        parts.append(node.id)
        return ".".join(reversed(parts))
    return None


#: Every importlib registration in tests, keyed by (file, name argument as written) -> what the loaded file is:
#: "entry" (GRIBS_v*.py) or "other" (a tool or harness, not GRIBS code). Default-deny: an unlisted registration fails
#: T9, so a new registration name is classified before T0 relies on it (report 116).
REGISTRATIONS = {
    ("tests/_cache_v7_child.py", "'gribs_child'"): "entry",
    ("tests/_cache_v7_sync_child.py", "'gribs_sync_child'"): "entry",
    ("tests/_a6_state_guard_cases.py", "'gribs_state_guard_cases'"): "entry",
    ("tests/_a6_state_guard_late_import_cases.py", "'gribs_state_guard_late'"): "entry",
    ("tests/_a6_two_cases_child.py", "'gribs_two_cases'"): "entry",
    ("tests/reference/test_reference_constant_backend.py", "'gribs_a5_reference'"): "entry",
    ("tests/test_a5_1_g3b_falsification.py", "'gribs_entry'"): "entry",
    ("tests/test_a5_1_g4a_tool_throat_validity.py", "'a5_1_g4a_cea_comparison'"): "other",
    ("tests/test_a5_1_g4a_tool_throat_validity.py", "name"): "other",
    ("tests/test_a6_reexports.py", "'gribs_entry'"): "entry",
    ("tests/test_cache_atomic_v7.py", "'gribs_cache_v7'"): "entry",
    ("tests/test_cea_theoretical_diagnostics_output.py", "'gribs_entry'"): "entry",
    ("tests/test_cleanup_p3_compat.py", "'gribs_entry'"): "entry",
    ("tests/test_cleanup_semantic_labels.py", "'gribs_entry'"): "entry",
    ("tests/test_cli_backend_override.py", "'gribs_cli_override'"): "entry",
    ("tests/test_cli_backend_override_boundary.py", "'gribs_cli_override_boundary'"): "entry",
    ("tests/test_cli_backend_override_edge.py", "'gribs_cli_override_edge'"): "entry",
    ("tests/test_convergence_sweep.py", "module_name"): "other",
    ("tests/test_gribs.py", "'gribs_entry'"): "entry",
    ("tests/test_integration_boundary_restart.py", "'gribs_boundary_restart'"): "entry",
    ("tests/test_legacy_psi_recursion_regression.py", "'gribs_entry'"): "entry",
    ("tests/test_n6_homogeneous_equals_legacy_totals.py", "'gribs_entry'"): "entry",
    ("tests/test_nonfinite_config_inputs.py", "'gribs_nonfinite_inputs'"): "entry",
    ("tests/test_selftest_backend_modes.py", "'gribs_entry'"): "entry",
}
ENTRY_NAMES = frozenset(ast.literal_eval(k[1]) for k, v in REGISTRATIONS.items() if v == "entry")

#: Default-deny remainder: sites whose target object T0 cannot resolve to the package or to the entry (objects made in
#: the test, fixture parameters, other modules such as os). Key (file, root expression) -> exact number of sites.
#: A new unresolved site fails T0 until it is classified here; a moved name at an unresolved site always fails.
UNRESOLVED = {
    ('tests/_cache_v7_child.py', 'os'): 1,
    ('tests/_cache_v7_sync_child.py', 'os'): 1,
    ('tests/test_a5_1_g4a_tool_throat_validity.py', 'fake'): 2,
    ('tests/test_a5_1_g4a_tool_throat_validity.py', 'ns'): 2,
    ('tests/test_a6_thermo_resolve_executable.py', 'backend'): 1,
    ('tests/test_cache_atomic_v7.py', 'os'): 34,
    ('tests/test_cleanup_semantic_labels.py', 'backend'): 1,
    ('tests/test_gribs.py', 'backend'): 49,
    ('tests/test_gribs.py', 'c'): 2,
    ('tests/test_gribs.py', 'config'): 2,
    ('tests/test_gribs.py', 'legacy'): 1,
    ('tests/test_gribs.py', 'result'): 12,
    ('tests/test_integration_boundary_restart.py', 'p_at_6e6'): 1,
    ('tests/test_integration_boundary_restart.py', 'p_down_6e6'): 2,
    ('tests/test_nonfinite_config_inputs.py', 'config'): 6,
}


def _rel(path):
    return path.relative_to(ROOT).as_posix()


def _const_str(node):
    return node.value if isinstance(node, ast.Constant) and isinstance(node.value, str) else None


def _by_module_name(name):
    """A module given by its import / sys.modules name -> ("package", dotted) | ("entry", name) | None."""
    if name is None:
        return None
    if name == "gribs" or name.startswith("gribs."):
        return ("package", name)
    if name in ENTRY_NAMES:
        return ("entry", name)
    return None


def _module_expr(node, bound):
    """What a module-valued expression denotes: ("package", dotted) | ("entry", name) | None (unresolved)."""
    if isinstance(node, ast.Name):
        return bound.get(node.id)
    if isinstance(node, ast.Subscript) and _chain(node.value) == "sys.modules":
        return _by_module_name(_const_str(node.slice))
    if isinstance(node, ast.Call):
        last = (_chain(node.func) or "").split(".")[-1]
        if last in ("import_module", "__import__") and node.args:
            return _by_module_name(_const_str(node.args[0]))
        if last == "module_from_spec":
            return ("entry", None)  # any loaded file counts as the entry (conservative; tools included)
    return None


def _bindings(tree):
    """Names bound anywhere in the file to a module: `import gribs...`, `from gribs[.x] import y`, and assignments from
    sys.modules[...], importlib.import_module(...), __import__(...) or module_from_spec(...). The entry wins over the
    package when one name is bound to both."""
    out = {}

    def bind(name, what):
        if what is not None and (out.get(name, ("",))[0] != "entry"):
            out[name] = what

    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for a in node.names:
                if a.name == "gribs" or a.name.startswith("gribs."):
                    bind(a.asname or "gribs", ("package", a.name if a.asname else "gribs"))
        elif isinstance(node, ast.ImportFrom) and node.level == 0 and \
                (node.module == "gribs" or (node.module or "").startswith("gribs.")):
            for a in node.names:
                bind(a.asname or a.name, ("package", f"{node.module}.{a.name}"))
        elif isinstance(node, (ast.Assign, ast.AnnAssign)) and node.value is not None:
            what = _module_expr(node.value, {})
            for t in (node.targets if isinstance(node, ast.Assign) else [node.target]):
                if isinstance(t, ast.Name):
                    bind(t.id, what)
    return out


def _moved_target(dotted):
    """Last component of "pkg.mod.NAME[.attr]" that is a moved name -> (name, prefix before it)."""
    parts = dotted.split(".")
    hit = None
    for i, part in enumerate(parts):
        if part in MOVED:
            hit = (part, ".".join(parts[:i]))
    return hit or (None, None)


class Site(NamedTuple):
    path: Path
    line: int
    sym: Optional[str]
    form: str
    owner_ok: bool
    root: str


def _split_root(node):
    """expr.a.b -> (root expr, ["a", "b"])."""
    parts = []
    while isinstance(node, ast.Attribute):
        parts.append(node.attr)
        node = node.value
    return node, list(reversed(parts))


def _classify(path, line, root, parts, bound, kind):
    """One site: target object `root` reached through attribute names `parts` (the last one is the name that is set;
    None when it is not a literal). kind is "object" (setattr / setitem) or "assign" (store / delete)."""
    if isinstance(root, ast.Name) and root.id in ("self", "cls"):
        return None  # attributes of instances / classes defined in the test itself
    what = _module_expr(root, bound)
    names = [p for p in parts if p is not None]
    moved = [p for p in names if p in MOVED]
    if what is not None and what[0] == "package":
        dotted = ".".join([what[1]] + names)
        sym, prefix = _moved_target(dotted)
        if sym is not None:
            return Site(path, line, sym, f"package-{kind}", kind == "object" and prefix == MOVED[sym][0], "")
        return Site(path, line, names[-1] if names else "<dynamic>", f"package-{kind}-unclassified", False, "")
    if what is not None and what[0] == "entry":
        sym = moved[-1] if moved else (names[0] if names else "<dynamic>")
        return Site(path, line, sym, f"entry-{kind}", False, "")
    # default-deny: unresolved target
    sym = moved[-1] if moved else None
    return Site(path, line, sym, f"unresolved-{kind}", False, ast.unparse(root))


def _sites(paths=None):
    """Every patch / assignment site in tests (default-deny, report 116):
    - every call whose name ends in setattr / delattr (monkeypatch.* or builtins), object form or string form;
    - monkeypatch.setitem / delitem on X.__dict__ or vars(X);
    - every attribute store / delete (Assign, AugAssign, AnnAssign, Delete; any depth), and item stores into
      X.__dict__ / vars(X).
    A target that resolves to the package or to the entry is classified; any other target is "unresolved" and is
    counted per (file, root) in UNRESOLVED; self / cls are not sites."""
    for path in (_test_files() if paths is None else paths):
        tree = ast.parse(path.read_text(), filename=str(path))
        bound = _bindings(tree)
        for node in ast.walk(tree):
            if isinstance(node, ast.Call) and node.args:
                last = (_chain(node.func) or "").split(".")[-1]
                first = node.args[0]
                if last in ("setattr", "delattr"):
                    string_form = (last == "setattr" and len(node.args) == 2) or (last == "delattr" and len(node.args) == 1
                                                                                  and isinstance(node.func, ast.Attribute))
                    if string_form:
                        target = _const_str(first)
                        if target is None:
                            yield Site(path, node.lineno, None, "unresolved-string", False, ast.unparse(first))
                            continue
                        head, *rest = target.split(".")
                        what = _by_module_name(head)
                        sym, prefix = _moved_target(target)
                        if what and what[0] == "package":
                            if sym is not None:
                                yield Site(path, node.lineno, sym, "string", prefix == MOVED[sym][0], "")
                            else:
                                yield Site(path, node.lineno, target.split(".")[-1], "string-unclassified", False, "")
                        elif what and what[0] == "entry":
                            yield Site(path, node.lineno, sym or (rest[0] if rest else "<dynamic>"), "string-entry", False, "")
                        else:
                            yield Site(path, node.lineno, sym, "unresolved-string", False, repr(head))
                        continue
                    name = _const_str(node.args[1]) if len(node.args) > 1 else None
                    root, parts = _split_root(first)
                    site = _classify(path, node.lineno, root, parts + [name], bound, "object")
                    if site:
                        yield site
                elif last in ("setitem", "delitem") and isinstance(node.func, ast.Attribute) and len(node.args) > 1:
                    obj = _namespace_of(first)
                    if obj is not None:
                        root, parts = _split_root(obj)
                        site = _classify(path, node.lineno, root, parts + [_const_str(node.args[1])], bound, "object")
                        if site:
                            yield site
            elif isinstance(node, (ast.Assign, ast.AugAssign, ast.AnnAssign, ast.Delete)):
                targets = node.targets if isinstance(node, (ast.Assign, ast.Delete)) else [node.target]
                for t in targets:
                    for x in ast.walk(t):
                        if isinstance(x, ast.Attribute) and isinstance(x.ctx, (ast.Store, ast.Del)):
                            root, parts = _split_root(x)
                        elif isinstance(x, ast.Subscript) and isinstance(x.ctx, (ast.Store, ast.Del)) \
                                and _namespace_of(x.value) is not None:
                            root, parts = _split_root(_namespace_of(x.value))
                            parts = parts + [_const_str(x.slice)]
                        else:
                            continue
                        site = _classify(path, node.lineno, root, parts, bound, "assign")
                        if site:
                            yield site


def _namespace_of(node):
    """X for X.__dict__ or vars(X); otherwise None."""
    if isinstance(node, ast.Attribute) and node.attr == "__dict__":
        return node.value
    if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == "vars" and node.args:
        return node.args[0]
    return None


@pytest.mark.fast
def test_t0_every_patch_site_is_classified():
    """Default-deny completeness: each site is a moved name patched on its owner (exactly MOVED's count), a not-yet-
    moved name counted exactly in ALLOWLIST, or an unresolved target counted exactly in UNRESOLVED."""
    assert not set(MOVED) & set(ALLOWLIST), sorted(set(MOVED) & set(ALLOWLIST))
    pending, on_owner, unresolved = Counter(), Counter(), Counter()
    problems = []
    for s in _sites():
        if s.sym in MOVED:
            if s.owner_ok:
                on_owner[s.sym] += 1
            else:
                problems.append(f"{s.path}:{s.line}: moved name {s.sym} patched as {s.form}, not on {MOVED[s.sym][0]}")
        elif s.form.startswith("unresolved"):
            unresolved[(_rel(s.path), s.root)] += 1
        else:
            pending[s.sym] += 1
    assert not problems, "\n".join(problems)
    expected = {name: count for name, (_owner, _kind, count) in MOVED.items()}
    assert dict(on_owner) == expected, {"moved_sites_found": dict(on_owner), "expected": expected}
    assert dict(pending) == ALLOWLIST, {"unexpected_or_changed": {k: v for k, v in pending.items() if ALLOWLIST.get(k) != v},
                                        "missing": {k: v for k, v in ALLOWLIST.items() if pending.get(k) != v}}
    assert dict(unresolved) == UNRESOLVED, {"unresolved_new_or_changed": {k: v for k, v in unresolved.items() if UNRESOLVED.get(k) != v},
                                            "missing": {k: v for k, v in UNRESOLVED.items() if unresolved.get(k) != v}}


@pytest.mark.fast
def test_t9_every_registration_is_classified():
    """Every importlib.util.spec_from_file_location in tests is listed in REGISTRATIONS (default-deny)."""
    found = set()
    for path in _test_files():
        for node in ast.walk(ast.parse(path.read_text(), filename=str(path))):
            if isinstance(node, ast.Call) and (_chain(node.func) or "").split(".")[-1] == "spec_from_file_location" \
                    and node.args:
                found.add((_rel(path), ast.unparse(node.args[0])))
    assert found == set(REGISTRATIONS), {"unlisted": sorted(found - set(REGISTRATIONS)),
                                         "stale": sorted(set(REGISTRATIONS) - found)}


@pytest.mark.fast
def test_t1_owner_module_defines_each_moved_name():
    for name, (owner, _kind, _count) in MOVED.items():
        assert hasattr(importlib.import_module(owner), name), (name, owner)


@pytest.mark.fast
def test_t2_moved_names_are_patched_only_on_their_owner():
    """Whitelist: a moved name may be patched only as monkeypatch.setattr("<owner>.<name>[.attr]", value) or with a
    first argument that is the owner module bound by `import gribs...` in that file. The entry (any alias) never."""
    problems = [f"{s.path}:{s.line}: {s.sym} ({s.form})" for s in _sites() if s.sym in MOVED and not s.owner_ok]
    assert not problems, "\n".join(problems)


@pytest.mark.fast
def test_t3_package_modules_do_not_from_import_moved_names():
    """Decision (A): inside the package a moved name is never from-imported."""
    problems = []
    for path in sorted(SRC.rglob("*.py")):
        tree = ast.parse(path.read_text(), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                for alias in node.names:
                    if alias.name in MOVED or alias.name == "*":
                        problems.append(f"{path}:{node.lineno}: from {'.' * node.level}{node.module or ''} import {alias.name}")
    assert not problems, "\n".join(problems)


@pytest.mark.fast
def test_t4_entry_does_not_read_moved_names_by_bare_name():
    tree = ast.parse(_entry().read_text())
    problems = [f"line {n.lineno}: {n.id}" for n in ast.walk(tree)
                if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Load)
                and MOVED.get(n.id, (None, None))[1] == "function_or_global"]
    assert not problems, "\n".join(problems)


@pytest.mark.fast
def test_t6_package_does_not_import_the_entry():
    """src/gribs/** never imports or loads the entry GRIBS_v*.py (no cycle back into the shim)."""
    problems = []
    for path in sorted(SRC.rglob("*.py")):
        text = path.read_text()
        tree = ast.parse(text, filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                problems += [f"{path}:{node.lineno}: import {a.name}" for a in node.names if a.name.startswith("gribs_entry")]
            elif isinstance(node, ast.ImportFrom) and (node.module or "").startswith("gribs_entry"):
                problems.append(f"{path}:{node.lineno}: from {node.module}")
        # no source text of the package names the entry file at all (catches spec_from_file_location, read_text, glob,
        # string building, comments and docstrings alike)
        for lineno, line in enumerate(text.splitlines(), 1):
            if "GRIBS_v" in line:
                problems.append(f"{path}:{lineno}: source mentions GRIBS_v")
    assert not problems, "\n".join(problems)


@pytest.mark.fast
def test_t7_module_names_do_not_collide_with_moved_names():
    """Static (before any import): no module or package basename under src/gribs equals a moved name, so a patch path
    like "gribs.<module>.<name>" is never ambiguous (e.g. the thermochemistry module is gribs.thermo)."""
    basenames = set()
    for path in SRC.rglob("*.py"):
        basenames.add(path.parent.name if path.name == "__init__.py" else path.stem)
    assert not basenames & set(MOVED), sorted(basenames & set(MOVED))


def _module_name(path):
    rel = path.relative_to(SRC).with_suffix("")
    parts = rel.parts[:-1] if rel.name == "__init__" else rel.parts
    return ".".join(parts)


@pytest.mark.fast
def test_t8_package_modules_import_only_same_or_earlier_units():
    problems = []
    for path in sorted(SRC.rglob("*.py")):
        me = _module_name(path)
        if me not in MODULE_UNIT:
            problems.append(f"{path}: module {me} is not in MODULE_UNIT")
            continue
        package = me if path.name == "__init__.py" else me.rpartition(".")[0]
        for node in ast.walk(ast.parse(path.read_text(), filename=str(path))):
            targets = []
            if isinstance(node, ast.Import):
                targets = [a.name for a in node.names]
            elif isinstance(node, ast.ImportFrom):
                if node.level:
                    base_parts = package.split(".")[: len(package.split(".")) - (node.level - 1)]
                    base = ".".join(base_parts + ([node.module] if node.module else []))
                else:
                    base = node.module or ""
                targets = [f"{base}.{a.name}" if f"{base}.{a.name}" in MODULE_UNIT else base for a in node.names]
            for target in targets:
                if not (target == "gribs" or target.startswith("gribs.")):
                    continue
                unit = MODULE_UNIT.get(target)
                if unit is None:
                    problems.append(f"{path}:{node.lineno}: imports {target}, which is not in MODULE_UNIT")
                elif unit > MODULE_UNIT[me]:
                    problems.append(f"{path}:{node.lineno}: {me} (unit {MODULE_UNIT[me]}) imports {target} (unit {unit})")
    assert not problems, "\n".join(problems)


def _names_state(node, name):
    """node mentions the state `name` as an attribute, an imported name or a string "<name>" / "<...>.<name>"."""
    if isinstance(node, ast.Attribute) and node.attr == name:
        return True
    if isinstance(node, ast.alias) and name in (node.name, node.name.rsplit(".", 1)[-1]):
        return True
    s = _const_str(node)
    return s is not None and (s == name or s.endswith("." + name))


def _accessor_violations_in_tests(paths):
    """Tests reach accessor-only state through the accessors only: no attribute, import or string of the name."""
    out = []
    for path in paths:
        for node in ast.walk(ast.parse(path.read_text(), filename=str(path))):
            for name in ACCESSOR_ONLY:
                if _names_state(node, name) or (isinstance(node, ast.Name) and node.id == name):
                    out.append(f"{_rel(path) if path.is_relative_to(ROOT) else path.name}:{getattr(node, 'lineno', '?')}: {name}")
    return out


def _accessor_violations_in_source(path):
    """In GRIBS code, the bare name is used only by its module-level definition and inside its two accessors; `global`
    only in the writer; never as an attribute or a string (globals()["_THERMO"], module._THERMO)."""
    tree = ast.parse(path.read_text(), filename=str(path))
    out = []
    allowed = set()
    for node in tree.body:
        if isinstance(node, (ast.Assign, ast.AnnAssign)):
            for t in (node.targets if isinstance(node, ast.Assign) else [node.target]):
                if isinstance(t, ast.Name) and t.id in ACCESSOR_ONLY:
                    allowed.add(id(t))
        if isinstance(node, ast.FunctionDef):
            for name, (reader, writer) in ACCESSOR_ONLY.items():
                if node.name in (reader, writer):
                    for x in ast.walk(node):
                        if isinstance(x, ast.Name) and x.id == name:
                            allowed.add(id(x))
                        if isinstance(x, ast.Global) and name in x.names and node.name == writer:
                            allowed.add(id(x))
    for node in ast.walk(tree):
        for name in ACCESSOR_ONLY:
            bad = (isinstance(node, ast.Name) and node.id == name and id(node) not in allowed) \
                or (isinstance(node, ast.Global) and name in node.names and id(node) not in allowed) \
                or _names_state(node, name)
            if bad:
                out.append(f"{path.name}:{getattr(node, 'lineno', '?')}: {name}")
    return out


@pytest.mark.fast
def test_t10_accessor_only_state():
    """_THERMO is read and written only by _thermo_state() / _set_thermo() (ADR decision 2): in tests (red excluded) the
    name never appears; in the entry and the package only the definition and the two accessors use it."""
    problems = _accessor_violations_in_tests(_test_files())  # conftest.py and the helpers included
    for path in [_entry()] + sorted(SRC.rglob("*.py")):
        problems += _accessor_violations_in_source(path)
    assert not problems, "\n".join(problems)
    for name, (reader, writer) in ACCESSOR_ONLY.items():
        definers = {}
        for path in [_entry()] + sorted(SRC.rglob("*.py")):
            where = "entry" if path == _entry() else _module_name(path)
            for n in ast.parse(path.read_text()).body:
                targets = [n.name] if isinstance(n, (ast.FunctionDef, ast.ClassDef)) else \
                    [t.id for t in (n.targets if isinstance(n, ast.Assign) else [n.target] if isinstance(n, ast.AnnAssign) else [])
                     if isinstance(t, ast.Name)]
                for d in targets:
                    if d in (name, reader, writer):
                        definers.setdefault(d, set()).add(where)
        owner = ACCESSOR_OWNER[name]
        assert definers == {name: {owner}, reader: {owner}, writer: {owner}}, definers


_T10_TEST_CASES = [
    ("x = gribs._THERMO", 1), ("gribs._THERMO = None", 1), ('monkeypatch.setattr(gribs, "_THERMO", None)', 1),
    ('monkeypatch.setattr("gribs_entry._THERMO", None)', 1), ('getattr(gribs, "_THERMO")', 1),
    ('vars(gribs)["_THERMO"] = None', 1), ("from gribs_entry import _THERMO", 1),
    ("gribs._set_thermo(None)\nx = gribs._thermo_state()", 0), ('k = gribs._THERMO_ROOT_KEYS\nc = "E_OLD_THERMO_CACHE"', 0),
]
_T10_SOURCE_CASES = [
    ("_THERMO = None\ndef _thermo_state():\n    return _THERMO\n"
     "def _set_thermo(b):\n    global _THERMO\n    p = _THERMO\n    _THERMO = b\n    return p", 0),
    ("_THERMO = None\ndef thermochemistry():\n    return _THERMO", 1),
    ("_THERMO = None\ndef _thermo_state():\n    global _THERMO\n    return _THERMO", 1),
    ('_THERMO = None\ndef f():\n    globals()["_THERMO"] = 1', 1),
    ("import sys\n_THERMO = None\ndef f():\n    sys.modules[__name__]._THERMO = 1", 1),
    ('"""``saved = _THERMO; _THERMO = x`` in a docstring is text."""\n_THERMO = None', 0),
]


@pytest.mark.fast
@pytest.mark.parametrize("source, expected", _T10_TEST_CASES)
def test_t10_classifier_tests(tmp_path, source, expected):
    path = tmp_path / "test_case.py"
    path.write_text(source + "\n")
    assert len(_accessor_violations_in_tests([path])) == expected


@pytest.mark.fast
@pytest.mark.parametrize("source, expected", _T10_SOURCE_CASES)
def test_t10_classifier_source(tmp_path, source, expected):
    path = tmp_path / "module.py"
    path.write_text(source + "\n")
    assert len(_accessor_violations_in_source(path)) == expected


# T0 classifier self-test (reports 115, 116): the form of every kind of site, including the forms of the review
# (S17c-S17g: the entry reached through sys.modules, importlib.import_module or an alias of it)
_E = 'import importlib, sys\n'
_CLASSIFIER_CASES = [
    ('monkeypatch.setattr("gribs.config.validate", f)', [("validate", "string", True)]),
    ('monkeypatch.setattr("gribs.constants.G0", 1)', [("G0", "string-unclassified", False)]),
    ('import gribs.constants as gc\nmonkeypatch.setattr(gc, "G0", 1)', [("G0", "package-object-unclassified", False)]),
    ("import gribs.constants as gc\ngc.G0 = 1", [("G0", "package-assign-unclassified", False)]),
    ('from gribs import config\nmonkeypatch.setattr(config, "validate", f)', [("validate", "package-object", True)]),
    (_E + 'm = importlib.import_module("gribs.config")\nm.validate = f', [("validate", "package-assign", False)]),
    ('monkeypatch.setattr("gribs_entry.G0", 1)', [("G0", "string-entry", False)]),
    ('monkeypatch.setattr("gribs_cli_override.validate", f)', [("validate", "string-entry", False)]),
    ('import importlib.util\ng = importlib.util.module_from_spec(spec)\nmonkeypatch.setattr(g, "validate", f)',
     [("validate", "entry-object", False)]),
    # S17c-S17g
    (_E + 'setattr(sys.modules["gribs_entry"], "validate", f)', [("validate", "entry-object", False)]),
    (_E + 'setattr(importlib.import_module("gribs_entry"), "validate", f)', [("validate", "entry-object", False)]),
    (_E + 'm = importlib.import_module("gribs_entry")\nmonkeypatch.setattr(m, "validate", f)',
     [("validate", "entry-object", False)]),
    (_E + 'm = importlib.import_module("gribs_entry")\nm.validate = f', [("validate", "entry-assign", False)]),
    (_E + 'sys.modules["gribs_entry"].validate = f', [("validate", "entry-assign", False)]),
    # namespaces, deletes, dynamic names, unresolved targets
    (_E + 'monkeypatch.setitem(vars(sys.modules["gribs_entry"]), "validate", f)', [("validate", "entry-object", False)]),
    (_E + 'sys.modules["gribs_entry"].__dict__["G0"] = 1', [("G0", "entry-assign", False)]),
    (_E + 'del sys.modules["gribs_entry"].G0', [("G0", "entry-assign", False)]),
    ('monkeypatch.delattr("gribs_entry.G0")', [("G0", "string-entry", False)]),
    ('import importlib.util\ng = importlib.util.module_from_spec(spec)\nsetattr(g, name, 1)',
     [("<dynamic>", "entry-object", False)]),
    ("def test_x(mod):\n    mod.validate = f", [("validate", "unresolved-assign", False)]),
    ("def test_x(mod):\n    mod.G0 = 1", [(None, "unresolved-assign", False)]),
    ('import os\nmonkeypatch.setattr(os, "sep", "/")', [(None, "unresolved-object", False)]),
    ('monkeypatch.setattr("os.sep", "/")\nmonkeypatch.setattr("gribs_a3_convergence_harness_test.x", 1)',
     [(None, "unresolved-string", False), (None, "unresolved-string", False)]),
    ("monkeypatch.setattr(target, 1)", [(None, "unresolved-string", False)]),
    ("class A:\n    def f(self):\n        self.validate = 1", []),
]


@pytest.mark.fast
@pytest.mark.parametrize("source, expected", _CLASSIFIER_CASES)
def test_t0_classifier_counts_each_form(tmp_path, source, expected):
    path = tmp_path / "test_case.py"
    path.write_text(source + "\n")
    assert [(x.sym, x.form, x.owner_ok) for x in _sites([path])] == expected
