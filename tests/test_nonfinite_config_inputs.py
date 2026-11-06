"""Non-finite configuration inputs are rejected; finite behaviour is unchanged.

P1  the configuration file parser rejects the non-standard JSON constants
    NaN / Infinity / -Infinity (RFC 8259).
P2  configuration_from_document() rejects any non-finite float in the
    document, including disabled sections and values such as 1e999 that are
    valid JSON numbers but overflow to inf in Python.
V   validate() rejects a non-finite value in any float-valued Config field,
    also for Config objects built or changed directly in Python (run_model()
    calls validate() first).

Unchanged on purpose: finite values (including finite values in disabled
sections, which stay unchecked and map to 0.0), every range rule, the
existing handling of huge integers, and the CLI --backend override.
"""
import copy
import dataclasses
import decimal
import fractions
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
MAIN = ROOT / "GRIBS_v0.6.1-alpha.py"
GOLD = ROOT / "cases" / "golden"
spec = importlib.util.spec_from_file_location("gribs_nonfinite_inputs", MAIN)
gribs = importlib.util.module_from_spec(spec)
sys.modules["gribs_nonfinite_inputs"] = gribs
spec.loader.exec_module(gribs)

# Tracked configuration files (explicit list: untracked files are never read).
CONFIG_FILES = [
    "gribs_config.json",
    "examples/example_config.json",
    "cases/golden/g1_reference_constant.json",
    "cases/golden/g2_single_phase_cea.json",
    "cases/golden/g3_low_condensed.json",
    "cases/golden/g4_medium_condensed.json",
    "cases/golden/g5_high_condensed_ap_al.json",
    "cases/golden/g6_cd_nozzle_blowdown.json",
    "cases/reference/level1_end_faces_0.json",
    "cases/reference/level1_end_faces_1.json",
    "cases/reference/level1_end_faces_2.json",
    "tests/fixtures/cli_backend_override/v030_1e50615.json",
    "tests/fixtures/cli_backend_override/v043_b1c9bfd.json",
    "Specification/gribs_config_eps1.json",
    "Specification/gribs_config_wide.json",
    "verification/a2-pre/gribs_config_eps4.json",
    "verification/a4-ledger/final-production/default/config.json",
    "verification/a4-ledger/final-production/igniter-bdf/config.json",
    "verification/a4-ledger/final-production/igniter-lsoda/config.json",
    "verification/a4-ledger/final-production/igniter-radau/config.json",
    "verification/a5-1-2b/rtol-1e-10/config.json",
    "verification/a5-1-2b/rtol-1e-7/config.json",
    "verification/a5-1-2c/gribs_config_wide_prefixed_evidence.json",
    "verification/a5-1-2d/baseline-rtol-1e-10/config.json",
    "verification/a5-1-2d/reproduction-rtol-1e-7/config.json",
    "verification/a5-1-2d/unit-1-rtol-1e-10/config.json",
    "verification/a5-1-2d/unit-1-rtol-1e-7/config.json",
    "verification/a5-1-2d/unit-2-rtol-1e-10/config.json",
    "verification/a5-1-2d/unit-2-rtol-1e-7/config.json",
]
NONFINITE = [float("nan"), float("inf"), float("-inf")]
NONFINITE_IDS = ["nan", "inf", "-inf"]
NUMERIC_FIELDS = [f.name for f in dataclasses.fields(gribs.Config) if f.type in (float, "float", int, "int")]


def _doc(case="g2"):
    return json.loads(next(GOLD.glob(f"{case}_*.json")).read_text(encoding="utf-8"))


def _config(case="g3"):
    return gribs.configuration_from_document(_doc(case)).config


def _validate_cli(path, *extra):
    return subprocess.run([sys.executable, str(MAIN), "--config", str(path), "--validate-config", *extra],
                          capture_output=True, text=True, timeout=120, check=False)


@pytest.mark.integration
@pytest.mark.parametrize("rel", CONFIG_FILES)
def test_tracked_configuration_files_still_validate(rel):
    p = _validate_cli(ROOT / rel)
    assert p.returncode == 0, p.stderr[-300:]


@pytest.mark.fast
@pytest.mark.parametrize("field", NUMERIC_FIELDS)
@pytest.mark.parametrize("value", NONFINITE, ids=NONFINITE_IDS)
def test_validate_rejects_nonfinite_numeric_config_field(field, value):
    config = _config()
    setattr(config, field, value)
    with pytest.raises(ValueError, match=f"{field} must be a finite number"):
        gribs.validate(config)


@pytest.mark.fast
@pytest.mark.parametrize("case", ["g1", "g2", "g3", "g4", "g5", "g6"])
def test_golden_configs_validate_unchanged(case):
    gribs.validate(_config(case))


@pytest.mark.fast
def test_change_after_validation_is_caught_by_run_model():
    config = _config("g1")
    gribs.validate(config)
    config.V_g0 = float("nan")
    with pytest.raises(ValueError, match="V_g0 must be a finite number"):
        gribs.run_model(config)


@pytest.mark.integration
@pytest.mark.parametrize("token", ["NaN", "Infinity", "-Infinity", "1e999", "-1e999"])
@pytest.mark.parametrize("old", ['"alpha": 0.0', '"pressure_points": 81', '"packing_fraction": 1.0'],
                         ids=["disabled-alpha", "pressure_points", "packing_fraction"])
def test_nonfinite_file_values_are_configuration_errors(tmp_path, token, old):
    text = next(GOLD.glob("g2_*.json")).read_text(encoding="utf-8")
    assert old in text
    path = tmp_path / "config.json"
    path.write_text(text.replace(old, old.split(":")[0] + f": {token}", 1), encoding="utf-8")
    p = _validate_cli(path)
    assert p.returncode == 2
    assert ("non-standard JSON constant" in p.stderr) or ("non-finite number" in p.stderr)


@pytest.mark.fast
@pytest.mark.parametrize("value", NONFINITE, ids=NONFINITE_IDS)
def test_nonfinite_python_value_rejected_in_disabled_section(value):
    doc = _doc()
    doc["propellant"]["burn_law"]["erosive_burning"]["alpha"] = value
    with pytest.raises(gribs.ConfigurationError, match="non-finite number"):
        gribs.configuration_from_document(copy.deepcopy(doc))


@pytest.mark.fast
@pytest.mark.parametrize("value", [-1.0, 0.0, 1e300])
def test_disabled_section_finite_values_keep_previous_behaviour(value):
    doc = _doc()
    doc["propellant"]["burn_law"]["erosive_burning"]["alpha"] = value
    doc["nozzle"]["throat_erosion"]["rate_m_s_at_reference_pressure"] = value
    config = gribs.configuration_from_document(doc).config
    gribs.validate(config)
    assert config.ero_alpha == 0.0 and config.ero_throat_c == 0.0


@pytest.mark.fast
def test_huge_integers_keep_previous_handling():
    doc = _doc()
    doc["propellant"]["packing_fraction"] = 10 ** 5000
    with pytest.raises(OverflowError):
        gribs.configuration_from_document(copy.deepcopy(doc))
    config = _config()
    config.V_g0 = 10 ** 5000          # an int, not a float: V does not apply
    gribs.validate(config)


@pytest.mark.integration
@pytest.mark.parametrize("backend, expected", [("reference_constant", 0), ("cea_python", 0), ("cea_legacy_executable", 2)])
def test_backend_override_unchanged(backend, expected):
    p = _validate_cli(next(GOLD.glob("g2_*.json")), "--backend", backend)
    assert p.returncode == expected, p.stderr[-300:]


# --- non-float numeric types (NumPy scalars / 0-d arrays / Decimal) ----------
# Finite values are never converted: a finite long double beyond the float
# range, huge ints, Fractions and finite NumPy scalars keep their previous
# handling; only NaN / +-Infinity are rejected.

NONFLOAT_NONFINITE = {
    "float32-nan": np.float32("nan"), "float32-inf": np.float32("inf"), "float32--inf": np.float32("-inf"),
    "float16-nan": np.float16("nan"), "float16--inf": np.float16("-inf"),
    "longdouble-nan": np.longdouble("nan"), "longdouble-inf": np.longdouble("inf"),
    "array0d-nan": np.array(np.nan), "array0d--inf": np.array(-np.inf),
    "Decimal-NaN": decimal.Decimal("NaN"), "Decimal-sNaN": decimal.Decimal("sNaN"),
    "Decimal-Infinity": decimal.Decimal("Infinity"), "Decimal--Infinity": decimal.Decimal("-Infinity"),
}


@pytest.mark.fast
@pytest.mark.parametrize("field", ["p_a", "V_g0", "rtol", "cea_pressure_points"])
@pytest.mark.parametrize("key", list(NONFLOAT_NONFINITE))
def test_validate_rejects_nonfinite_nonfloat_types(field, key):
    config = _config()
    setattr(config, field, NONFLOAT_NONFINITE[key])
    with pytest.raises(ValueError, match=f"{field} must be a finite number"):
        gribs.validate(config)


@pytest.mark.fast
@pytest.mark.parametrize("value", [np.longdouble("1e400"), np.float32(1.5), np.float16(1.5), np.array(1.5),
                                   decimal.Decimal("1.5"), fractions.Fraction(3, 2), 10 ** 5000],
                         ids=["longdouble-1e400", "float32", "float16", "array0d", "Decimal", "Fraction", "huge-int"])
def test_validate_keeps_finite_nonfloat_values_unchanged(value):
    config = _config()
    config.V_g0 = value
    gribs.validate(config)          # accepted before and after: finiteness only, no float conversion
    assert config.V_g0 is value


@pytest.mark.fast
@pytest.mark.parametrize("path", [("propellant", "burn_law", "erosive_burning", "alpha"),
                                  ("nozzle", "throat_erosion", "rate_m_s_at_reference_pressure"),
                                  ("environment", "ambient_pressure_Pa"),
                                  ("propellant", "reactants", 0, "wt_percent")],
                         ids=["disabled-alpha", "disabled-throat-rate", "enabled-ambient", "list-wt-percent"])
@pytest.mark.parametrize("key", list(NONFLOAT_NONFINITE))
def test_document_rejects_nonfinite_nonfloat_types_in_every_section(path, key):
    doc = _doc()
    node = doc
    for k in path[:-1]:
        node = node[k]
    node[path[-1]] = NONFLOAT_NONFINITE[key]
    with pytest.raises(gribs.ConfigurationError, match="non-finite number"):
        gribs.configuration_from_document(doc)


@pytest.mark.fast
@pytest.mark.parametrize("value", [np.float32(-1.0), decimal.Decimal(-1), np.array(-1.0)],
                         ids=["float32", "Decimal", "array0d"])
def test_disabled_section_finite_nonfloat_values_keep_previous_behaviour(value):
    doc = _doc()
    doc["propellant"]["burn_law"]["erosive_burning"]["alpha"] = value
    config = gribs.configuration_from_document(doc).config
    gribs.validate(config)
    assert config.ero_alpha == 0.0


@pytest.mark.fast
def test_enabled_field_finite_nonfloat_value_keeps_previous_type_error():
    doc = _doc()
    doc["environment"]["ambient_pressure_Pa"] = np.float32(101325.0)
    with pytest.raises(gribs.ConfigurationError, match="expected a number"):
        gribs.configuration_from_document(doc)


@pytest.mark.fast
@pytest.mark.parametrize("key", ["float32-nan", "Decimal-Infinity", "array0d-nan"])
def test_initialize_thermochemistry_alone_rejects_nonfinite(tmp_path, key):
    config = _config("g1")
    config.p_a = NONFLOAT_NONFINITE[key]
    with pytest.raises(ValueError, match="p_a must be a finite number"):
        gribs.initialize_thermochemistry(config, tmp_path / "cache")


@pytest.mark.fast
def test_initialize_before_validate_order_still_works(tmp_path):
    config = _config("g1")
    gribs.initialize_thermochemistry(config, tmp_path / "cache")   # finiteness only, no range rules
    gribs.validate(config)
    assert gribs.run_model(config)["stop_burn"] in ("burnout", "ambient", "time-limit")
