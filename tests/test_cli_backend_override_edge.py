"""Edge-case regression tests for the CLI --backend fix candidate v4 (round 12, work copy only).

Pins, for an override applied before validation:
  * the file must still name a valid backend of its own (``_backend_choice``: never assumed);
  * reference_constant values missing from the file come from the documented Config
    defaults (absent section or missing keys), as before the fix (Phase A-5 contract);
  * migrated (v0.3.0 / v0.4.3 / v0.4.4) documents: original_configuration is the migrated
    document without the override, the migration notes are unchanged and followed by the
    override note, and resolved_configuration carries the override.
Fixtures v030_1e50615.json / v043_b1c9bfd.json are gribs_config.json from git history
(commits 1e50615 and b1c9bfd).  No CEA run, no fcea2.
"""
import copy
import importlib.util
import json
import sys
from dataclasses import asdict
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = sorted(ROOT.glob("GRIBS_v*.py"))[-1]
spec = importlib.util.spec_from_file_location("gribs_cli_override_edge", MODULE_PATH)
gribs = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = gribs
spec.loader.exec_module(gribs)

pytestmark = pytest.mark.fast

G2 = ROOT / "cases" / "golden" / "g2_single_phase_cea.json"
FIX = Path(__file__).resolve().parent / "fixtures" / "cli_backend_override"
MIGRATED = [FIX / "v030_1e50615.json", FIX / "v043_b1c9bfd.json",
            ROOT / "Specification" / "gribs_config_wide.json"]
OVERRIDE_NOTE = "thermochemistry.backend overridden on the command line to {!r}."
DEFAULTS = (3000.0, 1.2, 25.0)


def _g2(tmp_path, mutate, name="cfg.json"):
    doc = json.loads(G2.read_text(encoding="utf-8"))
    doc["output"]["directory"] = str(tmp_path / "out")
    doc["thermochemistry"].pop("reference_constant", None)
    mutate(doc)
    path = tmp_path / name
    path.write_text(json.dumps(doc), encoding="utf-8")
    return path, doc


def _ref_values(config):
    return (config.reference_temperature_K, config.reference_gamma,
            config.reference_molecular_weight_kg_kmol)


def _first_error_line(path, override=None):
    with pytest.raises(gribs.ConfigurationError) as info:
        if override is None:
            gribs.load_json_configuration(path)
        else:
            gribs.load_json_configuration(path, backend_override=override)
    return str(info.value).splitlines()[0]


# 1. reference_constant section empty / partial / invalid under --backend reference_constant
@pytest.mark.parametrize("section, expected", [
    ({}, DEFAULTS),
    ({"temperature_K": 2800.0}, (2800.0, 1.2, 25.0)),
    ({"gamma": 1.3, "molecular_weight_kg_kmol": 22.0}, (3000.0, 1.3, 22.0)),
    ({"temperature_K": 2800.0, "gamma": 1.3, "molecular_weight_kg_kmol": 22.0}, (2800.0, 1.3, 22.0)),
])
def test_reference_override_fills_missing_keys_from_documented_defaults(tmp_path, section, expected):
    path, doc = _g2(tmp_path, lambda d: d["thermochemistry"].__setitem__("reference_constant", copy.deepcopy(section)))
    loaded = gribs.load_json_configuration(path, backend_override="reference_constant")
    assert loaded.config.thermo_backend == "reference_constant"
    assert _ref_values(loaded.config) == expected
    assert loaded.document == doc                                   # recorded section stays as written
    assert loaded.document["thermochemistry"]["reference_constant"] == section


@pytest.mark.parametrize("section, message", [
    ({"gamma": 0.9}, "thermochemistry.reference_constant.gamma: must be > 1"),
    ({"temperature_K": 2800.0, "foo": 1}, "thermochemistry.reference_constant: unknown key(s)"),
    (None, "thermochemistry.reference_constant: expected a JSON object"),
    ([], "thermochemistry.reference_constant: expected a JSON object"),
])
def test_reference_override_still_rejects_invalid_sections(tmp_path, capsys, section, message):
    path, _ = _g2(tmp_path, lambda d: d["thermochemistry"].__setitem__("reference_constant", section))
    assert gribs.main(["--config", str(path), "--backend", "reference_constant", "--validate-config"]) == 2
    assert message in capsys.readouterr().err


# 2. the file's own backend: missing / invalid / empty is rejected even with --backend
_MISSING = object()


@pytest.mark.parametrize("value", [_MISSING, "foo", "", "   ", "legacy_fit", 123, None, "cea2"],
                         ids=["missing", "unknown", "empty", "blank", "removed", "int", "null", "v03-alias"])
@pytest.mark.parametrize("override", ["cea_python", "reference_constant", "cea_legacy_executable"])
def test_file_backend_is_still_validated_under_override(tmp_path, capsys, value, override):
    def mutate(d):
        if value is _MISSING:
            d["thermochemistry"].pop("backend")
        else:
            d["thermochemistry"]["backend"] = value
    path, _ = _g2(tmp_path, mutate)
    without = _first_error_line(path)
    assert without.startswith("thermochemistry.backend")
    assert _first_error_line(path, override) == without           # same error as without override
    assert gribs.main(["--config", str(path), "--backend", override, "--validate-config"]) == 2
    assert without in capsys.readouterr().err


# 3. migrated documents: original / resolved configuration and notes
@pytest.mark.parametrize("fixture", MIGRATED, ids=lambda p: p.stem)
@pytest.mark.parametrize("override", ["reference_constant", "cea_python"])
def test_migrated_document_original_resolved_and_notes(fixture, override):
    plain = gribs.load_json_configuration(fixture)
    over = gribs.load_json_configuration(fixture, backend_override=override)
    assert plain.document["schema_version"] == gribs.SCHEMA_VERSION      # migrated, not the file
    assert over.document == plain.document                              # original = migrated doc, no override
    assert over.config.original_configuration == plain.document
    assert over.migration_notes == plain.migration_notes + [OVERRIDE_NOTE.format(override)]
    assert plain.migration_notes                                         # a migration did happen
    assert over.config.thermo_backend == override
    a, b = asdict(plain.config), asdict(over.config)
    changed = {k for k in a if a[k] != b[k] and not (a[k] != a[k] and b[k] != b[k])}
    allowed = {"thermo_backend", "cea_cache_enabled", "cea_rebuild_cache", "cea_cache_dir"}
    assert changed <= allowed, changed - allowed
    if override == "reference_constant":
        assert _ref_values(over.config) == DEFAULTS                      # A-5 defaults (no section in fixtures)


def test_migrated_v030_legacy_override_keeps_single_phase_and_validates(capsys):
    fixture = FIX / "v030_1e50615.json"
    plain = gribs.load_json_configuration(fixture)
    assert plain.config.two_phase_mode == "single_phase_legacy"
    assert gribs.main(["--config", str(fixture), "--backend", "cea_legacy_executable", "--validate-config"]) == 0


@pytest.mark.parametrize("fixture", [FIX / "v043_b1c9bfd.json", ROOT / "Specification" / "gribs_config_wide.json"],
                         ids=lambda p: p.stem)
def test_migrated_homogeneous_document_rejects_legacy_override(fixture, capsys):
    assert gribs.load_json_configuration(fixture).config.two_phase_mode == "homogeneous_equilibrium"
    assert gribs.main(["--config", str(fixture), "--backend", "cea_legacy_executable", "--validate-config"]) == 2
    assert "cannot be used with backend 'cea_legacy_executable'" in capsys.readouterr().err
