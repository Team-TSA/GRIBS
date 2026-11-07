"""A-6 stage 1 guard (ADR_A6_PACKAGE_SPLIT): names moved into the package are patched only on their owner.

MOVED grows, in the same commit as each move, by exactly the names that unit moves (report 109 table,
column `unit`). It is empty in unit 0, so every check below is inert until names are added.
tests/red/ is a frozen historical record (A-5.1-3) and is not checked.
"""
import ast
import importlib
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"

# name -> (owner module fully qualified, kind); kind is "function_or_global" or "class"
MOVED = {
}


def _entry():
    bodies = sorted(ROOT.glob("GRIBS_v*.py"))
    if len(bodies) != 1:
        raise RuntimeError(f"expected exactly one GRIBS_v*.py, found {[b.name for b in bodies]}")
    return bodies[0]


def _test_files():
    tests = ROOT / "tests"
    for path in sorted(tests.rglob("*.py")):
        rel = path.relative_to(tests)
        if rel.parts[0] == "red" or path.resolve() == Path(__file__).resolve():
            continue
        yield path


def _moved_target(dotted):
    """For "pkg.mod.NAME" or "pkg.mod.Class.attr", return (NAME-or-Class, prefix) if NAME/Class is moved."""
    parts = dotted.split(".")
    for i, part in enumerate(parts):
        if part in MOVED:
            return part, ".".join(parts[:i])
    return None, None


@pytest.mark.fast
def test_t1_owner_module_defines_each_moved_name():
    for name, (owner, _kind) in MOVED.items():
        assert hasattr(importlib.import_module(owner), name), (name, owner)


@pytest.mark.fast
def test_t2_moved_names_are_patched_only_through_the_owner_dotted_path():
    """Whitelist: for a moved name the only accepted patch is monkeypatch.setattr("<owner>.<name>[.attr]", value).
    Any setattr with an object as first argument, and any assignment to <expr>.<name>, fails."""
    problems = []
    for path in _test_files():
        tree = ast.parse(path.read_text(), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.Call):
                func = node.func
                is_setattr = (isinstance(func, ast.Attribute) and func.attr == "setattr") or \
                    (isinstance(func, ast.Name) and func.id == "setattr")
                if not is_setattr or not node.args:
                    continue
                first = node.args[0]
                if isinstance(first, ast.Constant) and isinstance(first.value, str):
                    name, prefix = _moved_target(first.value)
                    if name is not None and prefix != MOVED[name][0]:
                        problems.append(f"{path}:{node.lineno}: {first.value!r} is not under {MOVED[name][0]}")
                    continue
                # object form: setattr(obj, "NAME", ...) or setattr(obj.Class, "attr", ...)
                second = node.args[1] if len(node.args) > 1 else None
                names = set()
                if isinstance(second, ast.Constant) and isinstance(second.value, str):
                    names.add(second.value)
                for sub in ast.walk(first):
                    if isinstance(sub, ast.Attribute):
                        names.add(sub.attr)
                    elif isinstance(sub, ast.Name):
                        names.add(sub.id)
                hit = sorted(names & set(MOVED))
                if hit:
                    problems.append(f"{path}:{node.lineno}: object-form setattr of moved name(s) {hit}")
            elif isinstance(node, (ast.Assign, ast.AugAssign, ast.AnnAssign)):
                targets = node.targets if isinstance(node, ast.Assign) else [node.target]
                for target in targets:
                    if isinstance(target, ast.Attribute) and target.attr in MOVED:
                        problems.append(f"{path}:{node.lineno}: assignment to .{target.attr} (moved name)")
    assert not problems, "\n".join(problems)


@pytest.mark.fast
def test_t3_package_modules_do_not_from_import_moved_names():
    """Decision (A): inside the package a moved name is never from-imported (its owner defines it; every
    other module calls it as <owner module>.<name>)."""
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
