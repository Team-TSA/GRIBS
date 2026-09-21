import importlib.util
import subprocess
import sys
from dataclasses import asdict
from pathlib import Path


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


def test_example_configuration_keys_are_valid():
    import json

    path = ROOT / "examples" / "example_config.json"
    supplied = json.loads(path.read_text(encoding="utf-8"))
    defaults = asdict(gribs.Config())
    assert set(supplied).issubset(defaults)
