"""The only version written to any output is PROGRAM_VERSION (ADR_A6_PACKAGE_SPLIT). The pyproject version is a
placeholder and importlib.metadata.version("gribs") is never used."""
import ast
import json
import subprocess
import sys
import tomllib
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _entry():
    bodies = sorted(ROOT.glob("GRIBS_v*.py"))
    if len(bodies) != 1:
        raise RuntimeError(f"expected exactly one GRIBS_v*.py, found {[b.name for b in bodies]}")
    return bodies[0]


def _program_version():
    tree = ast.parse((ROOT / "src" / "gribs" / "constants.py").read_text())
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == "PROGRAM_VERSION" for t in node.targets):
            return node.value.value
    raise AssertionError("PROGRAM_VERSION not found")


def _pyproject_version():
    return tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]["version"]


@pytest.mark.fast
def test_no_metadata_version_of_gribs_in_sources():
    problems = []
    for path in [_entry(), *sorted((ROOT / "src").rglob("*.py"))]:
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.Call) and node.args and isinstance(node.args[0], ast.Constant) \
                    and node.args[0].value == "gribs":
                func = node.func
                name = func.attr if isinstance(func, ast.Attribute) else getattr(func, "id", "")
                if name in ("version", "metadata", "distribution"):
                    problems.append(f"{path}:{node.lineno}")
    assert not problems, problems


@pytest.mark.integration
def test_cli_version_and_outputs_carry_program_version_only(tmp_path):
    pv, pp = _program_version(), _pyproject_version()
    assert pv != pp
    r = subprocess.run([sys.executable, "-B", "-P", "-s", str(_entry()), "--version"], capture_output=True, text=True, timeout=120)
    assert r.returncode == 0 and r.stdout.strip() == f"GRIBS {pv}", (r.stdout, r.stderr)
    doc = json.loads((ROOT / "cases" / "golden" / "g1_reference_constant.json").read_text())
    doc["output"]["directory"] = str(tmp_path / "out")
    cfg = tmp_path / "g1.json"; cfg.write_text(json.dumps(doc))
    r = subprocess.run([sys.executable, "-B", "-P", "-s", str(_entry()), "--config", str(cfg)], capture_output=True, text=True,
                       timeout=600, cwd=str(tmp_path))
    assert r.returncode == 0, r.stderr[-3000:]
    summary = json.loads((tmp_path / "out" / "summary.json").read_text())
    assert summary["program"]["version"] == pv
    for name in ("summary.json", "summary.txt"):
        assert pp not in (tmp_path / "out" / name).read_text(), name
