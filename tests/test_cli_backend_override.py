"""Regression-test candidate v3 (round 10; v2 was round 9, work copy only — NOT in the repo).

`--backend` must be validated as the *effective* configuration: the CLI
override has to go through the same document-level rules as a backend written
in the file (legacy/homogeneous exclusion, backend-conditional sections).
All tests are fcea2-free: they stop at --validate-config or at the loader.
"""
import copy
import importlib.util
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = sorted(ROOT.glob("GRIBS_v*.py"))[-1]
spec = importlib.util.spec_from_file_location("gribs_cli_override", MODULE_PATH)
gribs = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = gribs
spec.loader.exec_module(gribs)

G2 = ROOT / "cases" / "golden" / "g2_single_phase_cea.json"


def _write(tmp_path, mutate, name="cfg.json"):
    doc = json.loads(G2.read_text(encoding="utf-8"))
    doc["output"]["directory"] = str(tmp_path / "out")
    mutate(doc)
    path = tmp_path / name
    path.write_text(json.dumps(doc), encoding="utf-8")
    return path, doc


def _drop(*sections):
    def mutate(doc):
        for s in sections:
            doc["thermochemistry"].pop(s, None)
    return mutate


def _validate(path, *extra):
    return gribs.main(["--config", str(path), "--validate-config", *extra])


pytestmark = pytest.mark.fast


def test_override_to_legacy_rejects_homogeneous_mode(tmp_path, capsys):
    path, doc = _write(tmp_path, _drop())
    assert doc["two_phase_model"]["mode"] == "homogeneous_equilibrium"
    assert _validate(path) == 0                       # control: file as written
    capsys.readouterr()
    assert _validate(path, "--backend", "cea_legacy_executable") == 2
    err = capsys.readouterr().err
    assert "cannot be used with backend 'cea_legacy_executable'" in err


def test_override_to_reference_without_section_uses_documented_a5_defaults(tmp_path):
    """v3: keeps the Phase A-5 contract (A5_DESIGN.md): documented Config defaults."""
    path, doc = _write(tmp_path, _drop("reference_constant"))
    loaded = gribs.load_json_configuration(path, backend_override="reference_constant")
    c = loaded.config
    assert (c.reference_temperature_K, c.reference_gamma, c.reference_molecular_weight_kg_kmol) == (3000.0, 1.2, 25.0)
    assert "reference_constant" not in loaded.document["thermochemistry"]   # recorded file unchanged
    assert loaded.config.original_configuration == doc


def test_override_to_cea_python_requires_cea_python_section(tmp_path, capsys):
    ref = json.loads((ROOT / "cases" / "golden" / "g1_reference_constant.json").read_text(encoding="utf-8"))
    ref["output"]["directory"] = str(tmp_path / "out"); ref["thermochemistry"].pop("cea_python")
    path = tmp_path / "ref.json"; path.write_text(json.dumps(ref), encoding="utf-8")
    assert _validate(path, "--backend", "cea_python") == 2
    assert "thermochemistry.cea_python" in capsys.readouterr().err


def test_override_to_legacy_requires_legacy_section(tmp_path, capsys):
    def mutate(doc):
        _drop("cea_legacy_executable")(doc)
        doc["two_phase_model"]["mode"] = "single_phase_legacy"
    path, _ = _write(tmp_path, mutate)
    assert _validate(path, "--backend", "cea_legacy_executable") == 2
    assert "cea_legacy_executable" in capsys.readouterr().err


def test_valid_override_is_accepted_and_reported(tmp_path, capsys):
    def mutate(doc):
        doc["thermochemistry"]["reference_constant"] = {
            "temperature_K": 2500.0, "gamma": 1.25, "molecular_weight_kg_kmol": 24.0}
    path, _ = _write(tmp_path, mutate)
    assert _validate(path, "--backend", "reference_constant") == 0
    captured = capsys.readouterr()
    assert "backend reference_constant" in captured.out
    assert "overridden on the command line to 'reference_constant'" in captured.err


def test_override_keeps_recorded_document_and_sets_effective_backend(tmp_path):
    def mutate(doc):
        doc["thermochemistry"]["reference_constant"] = {
            "temperature_K": 2500.0, "gamma": 1.25, "molecular_weight_kg_kmol": 24.0}
    path, doc = _write(tmp_path, mutate)
    loaded = gribs.load_json_configuration(path, backend_override="reference_constant")
    assert loaded.config.thermo_backend == "reference_constant"
    assert loaded.config.reference_temperature_K == 2500.0
    assert loaded.document["thermochemistry"]["backend"] == doc["thermochemistry"]["backend"]
    assert loaded.config.original_configuration == doc
    assert any("overridden on the command line" in n for n in loaded.migration_notes)


def test_no_override_is_unchanged(tmp_path):
    path, doc = _write(tmp_path, _drop())
    a = gribs.load_json_configuration(path)
    b = gribs.load_json_configuration(path, backend_override=None)
    assert a.config == b.config and a.document == b.document == doc
    assert a.migration_notes == b.migration_notes == []


@pytest.mark.parametrize("case", sorted((ROOT / "cases" / "golden").glob("g[0-9]_*.json")))
def test_golden_configs_with_their_own_backend_still_validate(case, tmp_path):
    doc = json.loads(case.read_text(encoding="utf-8"))
    backend = doc["thermochemistry"]["backend"]
    doc["output"]["directory"] = str(tmp_path / "out")
    path = tmp_path / case.name
    path.write_text(json.dumps(doc), encoding="utf-8")
    assert _validate(path) == 0
    assert _validate(path, "--backend", backend) == 0
