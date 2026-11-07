"""A-6 stage 1 guard (ADR_A6_PACKAGE_SPLIT): names moved into the package are patched only on their owner.

MOVED grows, in the same commit as each move, by exactly the names that unit moves; ALLOWLIST (names not yet
moved, with their exact number of patch / assignment sites in tests) shrinks by the same names and is empty at
the end of stage 1. tests/red/ is a frozen historical record (A-5.1-3) and is not checked.
"""
import ast
import importlib
from collections import Counter
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
TESTS = ROOT / "tests"

# name -> (owner module fully qualified, kind); kind is "function_or_global" or "class"
MOVED = {
}

# names not moved yet -> exact number of patch / assignment sites in tests (AST census, report 110 / 111)
ALLOWLIST = {
    "INTEGRATION_MAX_RHS_CALLS": 1,
    "INTEGRATION_MAX_TABLE_RESTARTS": 1,
    "INTEGRATION_MAX_WALL_S": 1,
    "NozzleTwoPhaseModel": 2,
    "PressureTableBackend": 1,
    "ReferenceConstantBackend": 1,
    "_THERMO": 7,
    "_ThermoSwap": 1,
    "_cache_rename_noreplace": 2,
    "_cache_write_checkpoint": 15,
    "_cea_theoretical_diagnostics": 1,
    "_phase_grid": 2,
    "_print_backend_banner": 3,
    "burn_rate": 1,
    "choke_limit_pressure": 2,
    "environment_info": 1,
    "equilibrium_pressure": 2,
    "gas_props": 2,
    "initialize_thermochemistry": 5,
    "is_extrapolated": 2,
    "load_json_configuration": 2,
    "make_figure": 1,
    "mass_energy_ledger": 2,
    "model_applicability_error": 1,
    "nozzle_state": 5,
    "nozzle_transition_diagnostics": 2,
    "phase_fractions": 2,
    "psi_and_deriv": 3,
    "run_model": 3,
    "sample": 1,
    "self_tests_post_run": 1,
    "summarize": 1,
    "summary_text": 1,
    "thermochemistry": 2,
    "validate": 2,
    "write_csv": 1,
}


# package module -> unit that introduced it (ADR_A6_PACKAGE_SPLIT; report 112 order). A module may import only
# gribs modules of the same or an earlier unit. Every module under src/gribs must be listed.
MODULE_UNIT = {
    "gribs": 0,
    "gribs.constants": 1,
    "gribs.errors": 1,
    "gribs.util": 1,
}


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


def _entry_aliases(tree):
    """Names bound to a loaded entry: x = importlib.util.module_from_spec(...) or x = sys.modules[...]."""
    out = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign):
            v = node.value
            src = (isinstance(v, ast.Call) and (_chain(v.func) or "").endswith("module_from_spec")) or \
                  (isinstance(v, ast.Subscript) and _chain(v.value) == "sys.modules")
            if src:
                out |= {t.id for t in node.targets if isinstance(t, ast.Name)}
    return out


def _package_roots(tree):
    """Names bound by `import gribs...` (the package), not by loading the entry -> the module they denote.
    `import gribs.thermo` binds gribs -> "gribs"; `import gribs.thermo as th` binds th -> "gribs.thermo"."""
    out = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for a in node.names:
                if a.name == "gribs" or a.name.startswith("gribs."):
                    out[a.asname] = a.name if a.asname else None
                    if not a.asname:
                        out["gribs"] = "gribs"
        elif isinstance(node, ast.ImportFrom) and node.level == 0 and \
                (node.module == "gribs" or (node.module or "").startswith("gribs.")):
            # from gribs import thermo as th -> th denotes "gribs.thermo"
            for a in node.names:
                out[a.asname or a.name] = f"{node.module}.{a.name}"
        elif isinstance(node, ast.Assign) and isinstance(node.value, ast.Call) and node.value.args \
                and (_chain(node.value.func) or "").split(".")[-1] == "import_module" \
                and isinstance(node.value.args[0], ast.Constant) and isinstance(node.value.args[0].value, str) \
                and (node.value.args[0].value == "gribs" or node.value.args[0].value.startswith("gribs.")):
            # mod = importlib.import_module("gribs.thermo") -> mod denotes "gribs.thermo"
            for t in node.targets:
                if isinstance(t, ast.Name):
                    out[t.id] = node.value.args[0].value
    out.pop(None, None)
    for alias in _entry_aliases(tree):
        out.pop(alias, None)
    return out


def _moved_target(dotted):
    """Last component of "pkg.mod.NAME[.attr]" that is a moved name -> (name, prefix before it)."""
    parts = dotted.split(".")
    hit = None
    for i, part in enumerate(parts):
        if part in MOVED:
            hit = (part, ".".join(parts[:i]))
    return hit or (None, None)


def _sites():
    """Every patch / assignment site in tests: (path, line, symbol, form, owner_ok)."""
    for path in _test_files():
        tree = ast.parse(path.read_text(), filename=str(path))
        aliases, pkgroots = _entry_aliases(tree), _package_roots(tree)
        for node in ast.walk(tree):
            if isinstance(node, ast.Call):
                name = _chain(node.func) or ""
                if name.split(".")[-1] != "setattr" or not node.args:
                    continue
                first = node.args[0]
                second = node.args[1] if len(node.args) > 1 else None
                if isinstance(first, ast.Constant) and isinstance(first.value, str):
                    sym, prefix = _moved_target(first.value)
                    if sym is not None:
                        yield path, node.lineno, sym, "string", prefix == MOVED[sym][0]
                    elif first.value == "gribs" or first.value.startswith("gribs."):
                        # a gribs target that names no moved name: counted, so T0 needs it in ALLOWLIST
                        yield path, node.lineno, first.value.split(".")[-1], "string-unclassified", False
                    continue
                chain = _chain(first)
                if chain is None or not isinstance(second, ast.Constant) or not isinstance(second.value, str):
                    continue
                root = chain.split(".")[0]
                if root in aliases:
                    # candidates: the last attribute of the chain (class attribute patches) and the name argument;
                    # the one that is a moved name wins, otherwise the chain attribute (or the argument)
                    parts = chain.split(".")
                    candidates = ([parts[-1]] if len(parts) > 1 else []) + [second.value]
                    sym = next((c for c in candidates if c in MOVED), candidates[0])
                    yield path, node.lineno, sym, "entry-object", False
                elif root in pkgroots:
                    parts = pkgroots[root].split(".") + chain.split(".")[1:] + [second.value]
                    sym, prefix = _moved_target(".".join(parts))
                    if sym is not None:
                        yield path, node.lineno, sym, "package-object", prefix == MOVED[sym][0]
                    else:
                        yield path, node.lineno, parts[-1], "package-object-unclassified", False
            elif isinstance(node, (ast.Assign, ast.AugAssign, ast.AnnAssign)):
                for t in (node.targets if isinstance(node, ast.Assign) else [node.target]):
                    if isinstance(t, ast.Attribute) and isinstance(t.value, ast.Name):
                        if t.value.id in aliases:
                            yield path, node.lineno, t.attr, "entry-assign", False
                        elif t.value.id in pkgroots:
                            sym, _prefix = _moved_target(f"{pkgroots[t.value.id]}.{t.attr}")
                            yield path, node.lineno, sym if sym is not None else t.attr, \
                                "package-assign" if sym is not None else "package-assign-unclassified", False


@pytest.mark.fast
def test_t0_every_patch_site_is_classified():
    """Each site is a moved name patched on its owner, or a not-yet-moved name counted exactly in ALLOWLIST."""
    assert not set(MOVED) & set(ALLOWLIST), sorted(set(MOVED) & set(ALLOWLIST))
    pending = Counter()
    problems = []
    for path, line, sym, form, owner_ok in _sites():
        if sym in MOVED:
            if not owner_ok:
                problems.append(f"{path}:{line}: moved name {sym} patched as {form}, not on {MOVED[sym][0]}")
        else:
            pending[sym] += 1
    assert not problems, "\n".join(problems)
    assert dict(pending) == ALLOWLIST, {"unexpected_or_changed": {k: v for k, v in pending.items() if ALLOWLIST.get(k) != v},
                                        "missing": {k: v for k, v in ALLOWLIST.items() if pending.get(k) != v}}


@pytest.mark.fast
def test_t1_owner_module_defines_each_moved_name():
    for name, (owner, _kind) in MOVED.items():
        assert hasattr(importlib.import_module(owner), name), (name, owner)


@pytest.mark.fast
def test_t2_moved_names_are_patched_only_on_their_owner():
    """Whitelist: a moved name may be patched only as monkeypatch.setattr("<owner>.<name>[.attr]", value) or with a
    first argument that is the owner module bound by `import gribs...` in that file. The entry (any alias) never."""
    problems = [f"{p}:{l}: {s} ({f})" for p, l, s, f, ok in _sites() if s in MOVED and not ok]
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
