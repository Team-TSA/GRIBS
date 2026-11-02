"""Cleanup 2026-11-01: semantic labels and evidence-honest status tokens.

Covers K-02/K-03/K-04/K-05/K-07/K-08/K-09/K-10 and inventory item 30.  No
physics, threshold, default, cache schema or cache key is changed.
"""
import copy
import importlib.util
import json
import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = sorted(ROOT.glob("GRIBS_v*.py"))[0]
if "gribs" in sys.modules:
    gribs = sys.modules["gribs"]
else:
    spec = importlib.util.spec_from_file_location("gribs", MODULE_PATH)
    gribs = importlib.util.module_from_spec(spec)
    sys.modules["gribs"] = gribs
    spec.loader.exec_module(gribs)

SOURCE = MODULE_PATH.read_text(encoding="utf-8")


@pytest.mark.fast
@pytest.mark.parametrize("active,mode,yg,expected", [
    (False, gribs.TWO_PHASE_SINGLE_LEGACY, 1.0, "regression_comparison_only"),
    (True, gribs.TWO_PHASE_SINGLE_LEGACY, 0.5, "regression_comparison_only"),
    (False, gribs.TWO_PHASE_HOMOGENEOUS, 1.0, "in_scope_not_validated"),
    (True, gribs.TWO_PHASE_HOMOGENEOUS, 0.995, "in_scope_not_validated"),
    (True, gribs.TWO_PHASE_HOMOGENEOUS, 0.99, "approximate_with_limitations"),
    (True, gribs.TWO_PHASE_HOMOGENEOUS, 0.78, "approximate_with_limitations"),
    (True, gribs.TWO_PHASE_HOMOGENEOUS, 0.65, "out_of_scope"),
    (True, gribs.TWO_PHASE_HOMOGENEOUS, 0.60, "out_of_scope"),
])
def test_model_status_never_claims_validated_and_thresholds_unchanged(active, mode, yg, expected):
    assert gribs.MODEL_STATUS_VALIDATED_YG_THRESHOLD == 0.99
    assert gribs.MODEL_STATUS_OUT_OF_SCOPE_YG_THRESHOLD == 0.65
    status = gribs.model_status_from_yg(active, mode, yg)
    assert status == expected
    assert status != "validated_in_scope"
    assert status in gribs.MODEL_STATUS_PREVIOUS_TOKENS


@pytest.mark.fast
def test_semantic_status_is_not_pass_from_schema_number_alone():
    assert '"pass" if getattr(backend, \'cache_schema\', 0) >= 4' not in SOURCE
    assert gribs.SEMANTIC_STATUS_PREVIOUS_TOKENS == {
        "schema_ok_semantic_probe_not_linked": "pass", "fail": "fail"}


@pytest.mark.fast
def test_wrong_gas_phase_molecular_weight_claims_removed():
    for claim in (
        "M = float(solution.M)          # gas-phase molecular weight",
        "non-physical gas-phase molecular weight",
        "must use the GAS-PHASE molecular weight",
        "use the GAS-PHASE molecular weight",
        "require the GAS-PHASE molecular weight",
        "EqSolution.density equals the gas-phase ideal",
        "must use the gas-phase molecular weight",
    ):
        assert claim not in SOURCE, claim
    # rejection behaviour of 'total_MW' is unchanged
    assert "is no \"\n                \"longer accepted (v0.5.0-alpha)" in SOURCE or \
        "longer accepted (v0.5.0-alpha)" in SOURCE


@pytest.mark.fast
def test_total_mw_still_rejected():
    c = gribs.Config(cea_py_molecular_weight="total_MW")
    backend = gribs.CEAPythonBackend.__new__(gribs.CEAPythonBackend)
    backend.c = c
    with pytest.raises(gribs.ConfigurationError, match="M_CEA"):
        backend._check_option_consistency()


@pytest.mark.integration
def test_cea_backend_labels_and_deprecated_aliases(tmp_path):
    cea = pytest.importorskip("cea")
    doc = json.loads((ROOT / "cases" / "golden" / "g4_medium_condensed.json").read_text())
    doc = copy.deepcopy(doc)
    doc["output"]["directory"] = str(tmp_path / "out")
    doc["thermochemistry"]["cea_python"]["cache_directory"] = str(tmp_path / "cache")
    conf = gribs.configuration_from_document(doc).config
    backend = gribs.initialize_thermochemistry(conf, tmp_path / "cache")
    pv = backend.property_validation
    # K-07: new correct keys and deprecated aliases carry identical values
    for side in ("min", "max"):
        new = pv[f"equilibrium_molecular_weight_{side}_kg_kmol"]
        old = pv[f"gas_phase_molecular_weight_{side}_kg_kmol"]
        assert new == old
        assert pv["deprecated_keys"][f"gas_phase_molecular_weight_{side}_kg_kmol"] == \
            f"equilibrium_molecular_weight_{side}_kg_kmol"
    diag = backend.table["diagnostics"]
    # K-08
    assert diag["density_note"].startswith("EqSolution.density = p*M_CEA/(Ru*T)")
    # K-09: CEA 3.3.4 EqSolution has no sonic_velocity attribute
    if not hasattr(cea.EqSolution, "sonic_velocity"):
        assert diag["sonic_velocity_status"] == "UNAVAILABLE"
        assert all(math.isnan(v) for v in diag["sonic_velocity_m_s"])
    # K-10
    assert "kmol-gas/kg-total" in backend.table["n_gas_per_kg_total_definition"]
