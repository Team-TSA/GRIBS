"""Boundary regression tests for the CLI --backend fix candidate v5 (round 13, work copy only).

Rule pinned here: an override never makes a file acceptable that is rejected without it,
and the reported error is the same.  Covers a file written for reference_constant with an
absent / empty / partial section (with the same and with another --backend), a file whose
own backend section is missing, and a file whose backend and two_phase mode are inconsistent.
No CEA run, no fcea2.
"""
import copy
import importlib.util
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = sorted(ROOT.glob("GRIBS_v*.py"))[-1]
spec = importlib.util.spec_from_file_location("gribs_cli_override_boundary", MODULE_PATH)
gribs = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = gribs
spec.loader.exec_module(gribs)

pytestmark = pytest.mark.fast

GOLDEN = ROOT / "cases" / "golden"
REF_FULL = {"temperature_K": 2800.0, "gamma": 1.3, "molecular_weight_kg_kmol": 22.0}
_DROP = object()


def _doc(base, mode=None, **thermo):
    doc = json.loads((GOLDEN / base).read_text(encoding="utf-8"))
    for key, value in thermo.items():
        if value is _DROP:
            doc["thermochemistry"].pop(key, None)
        else:
            doc["thermochemistry"][key] = copy.deepcopy(value)
    if mode is not None:
        doc["two_phase_model"]["mode"] = mode
    return doc


def _write(tmp_path, doc):
    doc = copy.deepcopy(doc)
    doc["output"]["directory"] = str(tmp_path / "out")
    path = tmp_path / "cfg.json"
    path.write_text(json.dumps(doc), encoding="utf-8")
    return path


def _error(path, override=None):
    with pytest.raises(gribs.ConfigurationError) as info:
        if override is None:
            gribs.load_json_configuration(path)
        else:
            gribs.load_json_configuration(path, backend_override=override)
    return str(info.value).splitlines()[0]


INVALID_FILES = {
    "ref_file_section_absent": (_doc("g1_reference_constant.json", reference_constant=_DROP),
                                ("reference_constant", "cea_python")),
    "ref_file_section_empty": (_doc("g1_reference_constant.json", reference_constant={}),
                               ("reference_constant", "cea_python")),
    "ref_file_section_partial_T": (_doc("g1_reference_constant.json", reference_constant={"temperature_K": 2800.0}),
                                   ("reference_constant", "cea_python")),
    "ref_file_section_partial_gM": (_doc("g1_reference_constant.json",
                                         reference_constant={"gamma": 1.3, "molecular_weight_kg_kmol": 22.0}),
                                    ("reference_constant", "cea_python")),
    "legacy_file_without_legacy_section": (_doc("g2_single_phase_cea.json", mode="single_phase_legacy",
                                                backend="cea_legacy_executable", cea_legacy_executable=_DROP),
                                           ("cea_python", "reference_constant")),
    "py_file_without_py_section": (_doc("g2_single_phase_cea.json", mode="single_phase_legacy", backend="cea_python",
                                        cea_python=_DROP, reference_constant=REF_FULL),
                                   ("reference_constant", "cea_legacy_executable")),
    "ref_file_without_ref_section_homogeneous": (_doc("g2_single_phase_cea.json", backend="reference_constant"),
                                                 ("cea_python",)),
    "legacy_file_with_homogeneous_mode": (_doc("g2_single_phase_cea.json", backend="cea_legacy_executable",
                                               reference_constant=REF_FULL),
                                          ("cea_python", "reference_constant", "cea_legacy_executable")),
}
PARAMS = [(name, ov) for name, (_, ovs) in INVALID_FILES.items() for ov in ovs]


@pytest.mark.parametrize("name, override", PARAMS, ids=[f"{n}->{o}" for n, o in PARAMS])
def test_override_never_accepts_a_file_rejected_without_it(tmp_path, capsys, name, override):
    path = _write(tmp_path, INVALID_FILES[name][0])
    without = _error(path)
    assert _error(path, override) == without
    assert gribs.main(["--config", str(path), "--validate-config"]) == 2
    assert gribs.main(["--config", str(path), "--backend", override, "--validate-config"]) == 2
    err = capsys.readouterr().err
    assert err.count(without) == 2


@pytest.mark.parametrize("override", [None, "reference_constant"])
def test_complete_reference_file_uses_its_own_values(tmp_path, override):
    path = _write(tmp_path, _doc("g1_reference_constant.json", reference_constant=REF_FULL))
    loaded = (gribs.load_json_configuration(path) if override is None
              else gribs.load_json_configuration(path, backend_override=override))
    c = loaded.config
    assert c.thermo_backend == "reference_constant"
    assert (c.reference_temperature_K, c.reference_gamma, c.reference_molecular_weight_kg_kmol) == (2800.0, 1.3, 22.0)


def test_same_backend_override_changes_nothing_but_the_note(tmp_path):
    path = _write(tmp_path, _doc("g1_reference_constant.json", reference_constant=REF_FULL))
    plain = gribs.load_json_configuration(path)
    over = gribs.load_json_configuration(path, backend_override="reference_constant")
    assert over.config == plain.config and over.document == plain.document
    assert over.migration_notes == plain.migration_notes + [
        "thermochemistry.backend overridden on the command line to 'reference_constant'."]


def test_complete_reference_file_overridden_to_cea_python(tmp_path):
    path = _write(tmp_path, _doc("g1_reference_constant.json", reference_constant=REF_FULL))
    loaded = gribs.load_json_configuration(path, backend_override="cea_python")
    assert loaded.config.thermo_backend == "cea_python"
    assert loaded.document["thermochemistry"]["backend"] == "reference_constant"
