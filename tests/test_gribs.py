import importlib.util
import json
import re
import subprocess
import sys
from dataclasses import asdict
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
MODULE_FILES = sorted(ROOT.glob("GRIBS_v*.py"))

if len(MODULE_FILES) != 1:
    raise RuntimeError(
        f"Expected exactly one GRIBS_v*.py file, found {len(MODULE_FILES)}: "
        f"{[path.name for path in MODULE_FILES]}"
    )

MODULE_PATH = MODULE_FILES[0]
spec = importlib.util.spec_from_file_location("gribs", MODULE_PATH)

if spec is None or spec.loader is None:
    raise ImportError(f"Could not load GRIBS module from {MODULE_PATH}")

gribs = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = gribs
spec.loader.exec_module(gribs)


@pytest.mark.integration
def test_default_self_tests_pass():
    result = subprocess.run(
        [sys.executable, str(MODULE_PATH), "--selftest"],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, (
        f"Self-test failed with exit code {result.returncode}\n"
        f"stdout:\n{result.stdout}\n"
        f"stderr:\n{result.stderr}"
    )
    assert "-> all checks passed" in result.stdout


@pytest.mark.integration
def test_cd_nozzle_transition_classification(tmp_path):
    """C-D nozzle mass flow passes T4 while thrust jump remains diagnostic."""
    source = ROOT / "gribs_config.json"
    supplied = json.loads(source.read_text(encoding="utf-8"))
    supplied["nozzle"]["expansion_ratio"] = 4.0

    config_path = tmp_path / "gribs_config_eps4.json"
    config_path.write_text(
        json.dumps(supplied, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    result = subprocess.run(
        [
            sys.executable,
            str(MODULE_PATH),
            "--config",
            str(config_path),
            "--selftest",
        ],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, (
        "Ae/At = 4 self-test should pass after A-2 classification.\n"
        f"exit code: {result.returncode}\n"
        f"stdout:\n{result.stdout}\n"
        f"stderr:\n{result.stderr}"
    )

    diagnostic_lines = [
        line
        for line in result.stdout.splitlines()
        if "C-D nozzle transition mass-flow continuity" in line
    ]
    assert len(diagnostic_lines) == 1, (
        "Expected exactly one C-D nozzle transition diagnostic, "
        f"found {len(diagnostic_lines)}:\n"
        + "\n".join(diagnostic_lines)
    )

    diagnostic = diagnostic_lines[0]
    assert "[PASS]" in diagnostic
    assert "Ae/At = 4" in diagnostic
    assert "known jump-model limitation" in diagnostic

    match = re.search(
        r"mdot jump = ([0-9.eE+-]+), F jump = ([0-9.eE+-]+)",
        diagnostic,
    )
    assert match is not None, (
        f"Could not parse C-D nozzle transition diagnostics:\n{diagnostic}"
    )

    mdot_jump = float(match.group(1))
    thrust_jump = float(match.group(2))

    assert mdot_jump == pytest.approx(4.02e-08, rel=1.0e-2)
    assert thrust_jump == pytest.approx(1.04, rel=1.0e-2)
    assert "-> all checks passed" in result.stdout
    assert "[FAIL]" not in result.stdout


@pytest.mark.fast
def test_geometry_identity_finite_difference():
    c = gribs.Config()
    x = 0.5 * gribs.web_thickness(c)
    h = gribs.web_thickness(c) * 1e-6
    derivative = (
        gribs.geometry(x + h, c)[3]
        - gribs.geometry(x - h, c)[3]
    ) / (2 * h)
    area = gribs.geometry(x, c)[2]
    assert abs(derivative / area - 1.0) < 1e-7


@pytest.mark.fast
def test_example_configuration_keys_are_valid():
    path = ROOT / "examples" / "example_config.json"
    supplied = json.loads(path.read_text(encoding="utf-8"))
    defaults = asdict(gribs.Config())
    assert set(supplied).issubset(defaults)
