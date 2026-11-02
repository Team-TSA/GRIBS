"""Cleanup P3 v2 (proposal, not reflected): cache / output / status compatibility.

Covers: per-backend semantic status + basis + previous token (B, C),
robust sonic-column classification (D), fact-only cache provenance (A),
table ranges only when a property table exists (E), and the outputs that
carry them: summary.json and --dump-thermo (F).
"""
import copy, importlib.util, json, math, sys
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = sorted(ROOT.glob("GRIBS_v*.py"))[0]
if "gribs" in sys.modules:
    gribs = sys.modules["gribs"]
else:
    spec = importlib.util.spec_from_file_location("gribs", MODULE_PATH)
    gribs = importlib.util.module_from_spec(spec); sys.modules["gribs"] = gribs
    spec.loader.exec_module(gribs)

P1_DENSITY_NOTE = ("EqSolution.density = p*M_CEA/(Ru*T) "
                   "[kg-total/m^3]: total mixture mass per gas-phase "
                   "volume (condensed-phase volume neglected); equals "
                   "the gas-phase ideal density p/(R_gas*T) divided "
                   "by Yg. It equals the two-phase mixture density "
                   "only under the negligible condensed-volume "
                   "assumption. Diagnostic only")


def _doc(tmp_path, case, rebuild=True, cache_enabled=True):
    doc = copy.deepcopy(json.loads((ROOT / "cases" / "golden" / case).read_text()))
    doc["output"]["directory"] = str(tmp_path / "out")
    cp = doc["thermochemistry"]["cea_python"]
    cp["cache_directory"] = str(tmp_path / "cache")
    cp["rebuild_cache"] = rebuild
    cp["cache_enabled"] = cache_enabled
    return doc


def _ref_config():
    doc = json.loads((ROOT / "cases" / "reference" / "level1_end_faces_0.json").read_text())
    return gribs.configuration_from_document(doc).config


def _run(tmp_path, doc, extra=()):
    cfg = tmp_path / "c.json"; cfg.write_text(json.dumps(doc))
    return gribs.main(["--config", str(cfg), *extra])


# ---------------------------------------------------------------- B / C
@pytest.mark.fast
def test_semantic_status_per_backend_kind():
    c = _ref_config()
    ref = gribs.ReferenceConstantBackend(c, c.reactants_spec_list, c.packing_fraction)
    syn = gribs._SyntheticTestBackend(c)
    assert syn.cache_schema == 4 and syn.uses_cea is False   # inherited schema
    tok, basis, prev = gribs.semantic_status_for_backend(ref)
    assert (tok, prev) == ("not_applicable_no_cea_table", "fail")
    assert "does not use CEA" in basis and "schema" not in basis
    tok, basis, prev = gribs.semantic_status_for_backend(syn)
    assert (tok, prev) == ("not_applicable_no_cea_table", "pass")
    assert "does not use CEA" in basis

    class Cea4: uses_cea = True; cache_schema = 4; name = "cea_python"
    class Cea3: uses_cea = True; cache_schema = 3; name = "cea_python"
    tok, basis, prev = gribs.semantic_status_for_backend(Cea4())
    assert (tok, prev) == ("schema_ok_semantic_probe_not_linked", "pass")
    assert "schema number 4 only" in basis
    tok, basis, prev = gribs.semantic_status_for_backend(Cea3())
    assert (tok, prev) == ("fail", "fail")
    assert "< 4" in basis


@pytest.mark.fast
def test_model_status_basis_matches_branch():
    leg = gribs.model_status_basis(False, gribs.TWO_PHASE_SINGLE_LEGACY)
    off = gribs.model_status_basis(False, gribs.TWO_PHASE_HOMOGENEOUS)
    thr = gribs.model_status_basis(True, gribs.TWO_PHASE_HOMOGENEOUS)
    assert "not applied" in leg and "regression" in leg
    assert "not applied" in off and "inactive" in off
    assert "0.99 / 0.65" in thr


# ---------------------------------------------------------------- D
N = 5
@pytest.mark.fast
@pytest.mark.parametrize("column,expected", [
    (None, "UNAVAILABLE"),
    ([], "UNAVAILABLE"),
    ([float("nan")] * N, "UNAVAILABLE"),
    ([None] * N, "UNAVAILABLE"),
    ([900.0] * N, "AVAILABLE"),
    ([900.0] * (N - 1), "INVALID"),            # length mismatch
    ([900.0, float("nan")] + [900.0] * (N - 2), "INVALID"),   # partial
    ([900.0, None] + [900.0] * (N - 2), "INVALID"),           # partial null
    ([0.0] * N, "INVALID"),                    # non-positive
    ([-5.0] * N, "INVALID"),
    ([float("inf")] * N, "INVALID"),
    (["900"] * N, "INVALID"),                  # non-numeric
    ([True] * N, "INVALID"),
    ("900", "INVALID"),                        # wrong type
    ({"a": 1}, "INVALID"),
])
def test_sonic_column_status(column, expected):
    status, reason = gribs._sonic_column_status(column, N)
    assert status == expected, reason
    assert reason


@pytest.mark.fast
def test_density_note_constant_is_p1_text():
    # new caches keep exactly the P1 diagnostic text
    assert gribs.CEA_DENSITY_NOTE == P1_DENSITY_NOTE


# ---------------------------------------------------------------- A / D integration
def _built_cache(tmp_path):
    doc = _doc(tmp_path, "g4_medium_condensed.json", rebuild=True)
    conf = gribs.configuration_from_document(doc).config
    b = gribs.initialize_thermochemistry(conf, tmp_path / "cache")
    return b


def _reuse(tmp_path):
    doc = _doc(tmp_path, "g4_medium_condensed.json", rebuild=False)
    conf = gribs.configuration_from_document(doc).config
    return gribs.initialize_thermochemistry(conf, tmp_path / "cache")


@pytest.mark.integration
def test_new_cache_and_reuse_of_same_writer(tmp_path):
    pytest.importorskip("cea")
    b = _built_cache(tmp_path)
    sd = b.metadata()["semantic_definitions"]
    assert sd["table_source"] == "generated_this_run"
    assert sd["table_fields_missing"] == []
    assert sd["table_density_note_matches_current"] is True
    assert sd["sonic_velocity_status"] == "UNAVAILABLE"      # CEA 3.3.4 EqSolution
    b2 = _reuse(tmp_path)
    sd2 = b2.metadata()["semantic_definitions"]
    assert b2.cache_used and sd2["table_source"] == "cache"
    assert sd2["table_fields_missing"] == [] and sd2["table_density_note_matches_current"]


@pytest.mark.integration
@pytest.mark.parametrize("mutation,expected_status", [
    ("drop_p1_fields", "UNAVAILABLE"),
    ("sonic_missing", "UNAVAILABLE"),
    ("sonic_null", "UNAVAILABLE"),
    ("sonic_partial", "INVALID"),
    ("sonic_short", "INVALID"),
    ("sonic_nonpositive", "INVALID"),
    ("sonic_huge_int", "INVALID"),        # JSON int beyond the float range
])
def test_reused_cache_with_altered_fields(tmp_path, mutation, expected_status):
    pytest.importorskip("cea")
    b = _built_cache(tmp_path)
    path = b.cache_path
    t = json.loads(path.read_text())
    d = t["diagnostics"]; n = len(t["pressure_Pa"])
    if mutation == "drop_p1_fields":
        d.pop("sonic_velocity_status"); d.pop("sonic_velocity_note")
        d["density_note"] = "EqSolution.density equals the gas-phase ideal density p/(Rg*T)"
        t.pop("n_gas_per_kg_total_definition")
    elif mutation == "sonic_missing":
        d.pop("sonic_velocity_m_s")
    elif mutation == "sonic_null":
        d["sonic_velocity_m_s"] = [None] * n
    elif mutation == "sonic_partial":
        d["sonic_velocity_m_s"] = [1000.0] * (n - 1) + [None]
    elif mutation == "sonic_short":
        d["sonic_velocity_m_s"] = [1000.0] * (n - 1)
    elif mutation == "sonic_nonpositive":
        d["sonic_velocity_m_s"] = [0.0] * n
    elif mutation == "sonic_huge_int":
        d["sonic_velocity_m_s"] = [10**400] * n
    path.write_text(json.dumps(t)); before = path.read_bytes(); mtime_before = path.stat().st_mtime_ns
    b2 = _reuse(tmp_path)
    assert b2.cache_used and b2.key == b.key
    sd = b2.metadata()["semantic_definitions"]
    assert sd["sonic_velocity_status"] == expected_status
    assert sd["density_definition"] == gribs.CEA_DENSITY_NOTE
    if mutation == "drop_p1_fields":
        assert set(sd["table_fields_missing"]) == {
            "sonic_velocity_status", "sonic_velocity_note", "n_gas_per_kg_total_definition"}
        assert sd["table_density_note_matches_current"] is False
        assert "creation" not in json.dumps(sd) and "before_2026" not in json.dumps(sd)
    assert path.read_bytes() == before             # cache file never rewritten
    assert path.stat().st_mtime_ns == mtime_before


@pytest.mark.integration
def test_run_with_invalid_sonic_cache_completes_and_reports(tmp_path):
    pytest.importorskip("cea")
    b = _built_cache(tmp_path)
    t = json.loads(b.cache_path.read_text())
    t["diagnostics"]["sonic_velocity_m_s"] = [1000.0] * 3
    b.cache_path.write_text(json.dumps(t))
    assert _run(tmp_path, _doc(tmp_path, "g4_medium_condensed.json", rebuild=False)) == 0
    s = json.loads((tmp_path / "out" / "summary.json").read_text())
    txt = json.dumps(s)
    assert '"sonic_velocity_status": "INVALID"' in txt


# ---------------------------------------------------------------- E / F outputs
@pytest.mark.integration
def test_cea_run_summary_and_dump(tmp_path):
    pytest.importorskip("cea")
    doc = _doc(tmp_path, "g4_medium_condensed.json", rebuild=True)
    assert _run(tmp_path, doc) == 0
    s = json.loads((tmp_path / "out" / "summary.json").read_text())
    v = s["results"]["validation"]
    assert v["semantic_status"] == "schema_ok_semantic_probe_not_linked"
    assert v["status_previous_tokens"]["semantic_status"] == "pass"
    assert v["validated_pressure_range_Pa"] == v["table_pressure_range_Pa"]
    assert v["validated_Yg_range"] == v["table_Yg_range"]
    assert "not an externally validated range" in v["validated_range_basis"]
    assert "G4a" not in v["validated_range_basis"]
    assert "0.99 / 0.65" in v["model_status_basis"]
    assert '"semantic_definitions"' in json.dumps(s)
    txt = (tmp_path / "out" / doc["output"]["summary_text_filename"]).read_text()
    assert "table pressure range [not validated]" in txt
    dump = tmp_path / "dump.json"
    assert _run(tmp_path, _doc(tmp_path, "g4_medium_condensed.json", rebuild=False),
                ["--dump-thermo", str(dump)]) == 0
    prov = json.loads(dump.read_text())["provenance"]
    assert prov["semantic_definitions"]["density_definition"] == gribs.CEA_DENSITY_NOTE


@pytest.mark.integration
def test_cache_disabled_run(tmp_path):
    pytest.importorskip("cea")
    doc = _doc(tmp_path, "g4_medium_condensed.json", rebuild=False, cache_enabled=False)
    assert _run(tmp_path, doc) == 0
    assert not (tmp_path / "cache").exists() or not any((tmp_path / "cache").iterdir())
    s = json.loads((tmp_path / "out" / "summary.json").read_text())
    assert '"table_source": "generated_this_run"' in json.dumps(s)


@pytest.mark.integration
def test_reference_constant_run(tmp_path):
    doc = _doc(tmp_path, "g1_reference_constant.json", rebuild=True)
    assert _run(tmp_path, doc) == 0
    v = json.loads((tmp_path / "out" / "summary.json").read_text())["results"]["validation"]
    assert v["semantic_status"] == "not_applicable_no_cea_table"
    assert v["status_previous_tokens"]["semantic_status"] == "fail"
    assert "does not use CEA" in v["semantic_status_basis"]
    assert v["table_pressure_range_Pa"] is None and v["table_Yg_range"] is None
    assert "no chamber-property table" in v["validated_range_basis"]
    # existing keys are kept unchanged for current readers
    assert "validated_pressure_range_Pa" in v and "validated_Yg_range" in v
    assert v["model_status"] == "regression_comparison_only"
    assert "regression" in v["model_status_basis"]


# ---------------------------------------------------------------- regression (2026-11-02 review)
class _FloatConversionRaises(float):
    def __float__(self):
        raise RuntimeError("no float")


class _ReprRaises:
    def __repr__(self):
        raise RuntimeError("no repr")


@pytest.mark.fast
@pytest.mark.parametrize("column", [
    [10**400] * N,                       # OverflowError in float()
    [_FloatConversionRaises(900.0)] * N,  # __float__ raises
    [_ReprRaises()] * N,                  # non-numeric and its repr raises
    [10**5000] * N,                       # beyond the int string-digit limit
], ids=["int_1e400", "float_conversion_raises", "repr_raises", "int_1e5000"])
def test_sonic_column_status_never_raises(column):
    status, reason = gribs._sonic_column_status(column, N)
    assert status == "INVALID"
    assert isinstance(reason, str) and len(reason) < 300


class _Stub:
    def __init__(self, **kw):
        self.pmin, self.pmax = 1.0e5, 8.0e6
        self.__dict__.update(kw)


@pytest.mark.fast
def test_reference_range_basis_states_the_analytic_values():
    pv = {"gas_mass_fraction_min": 1.0, "gas_mass_fraction_max": 1.0}
    f = gribs._validation_range_fields(_Stub(), pv, False)
    assert f["validated_Yg_range"] == [1.0, 1.0]
    assert f["table_Yg_range"] is None and f["table_pressure_range_Pa"] is None
    assert "analytic phase-fraction range [1, 1]" in f["validated_range_basis"]
    assert "defaults" not in f["validated_range_basis"]


@pytest.mark.fast
@pytest.mark.parametrize("pv, validated, table", [
    ({"gas_mass_fraction_min": 0.9}, [0.9, 1.0], None),                 # max missing
    ({"gas_mass_fraction_max": 0.95}, [0.0, 0.95], None),               # min missing
    ({}, [0.0, 1.0], None),
    ({"gas_mass_fraction_min": 0.8, "gas_mass_fraction_max": 1.0}, [0.8, 1.0], [0.8, 1.0]),
], ids=["max_missing", "min_missing", "empty", "complete"])
def test_range_fields_with_a_table_and_incomplete_yg(pv, validated, table):
    """validated_Yg_range keeps the former expression; table_Yg_range is built
    only from a complete, finite range and never raises."""
    f = gribs._validation_range_fields(_Stub(), pv, True)
    assert f["validated_Yg_range"] == validated
    assert f["table_Yg_range"] == table
    assert f["table_pressure_range_Pa"] == [1.0e5, 8.0e6]
    assert "chamber-property table range" in f["validated_range_basis"]


@pytest.mark.fast
def test_range_fields_without_table_and_without_yg_use_defaults_wording():
    f = gribs._validation_range_fields(_Stub(), {}, False)
    assert f["validated_Yg_range"] == [0.0, 1.0]
    assert "defaults 0 / 1" in f["validated_range_basis"]


@pytest.mark.fast
@pytest.mark.parametrize("schema", [None, "4", True, float("nan"), _ReprRaises()],
                         ids=["none", "string", "bool", "nan", "repr_raises"])
def test_unknown_cache_schema_is_fail_and_never_raises(schema):
    cea_backend = _Stub(name="cea_python", uses_cea=True, cache_schema=schema)
    tok, basis, prev = gribs.semantic_status_for_backend(cea_backend)
    assert (tok, prev) == ("fail", "fail")
    assert "unknown cache schema" in basis
    ref = _Stub(name="reference_constant", uses_cea=False, cache_schema=schema)
    tok, basis, prev = gribs.semantic_status_for_backend(ref)
    assert tok == "not_applicable_no_cea_table" and prev == "fail"


@pytest.mark.fast
@pytest.mark.parametrize("schema, tok, prev", [(4, "schema_ok_semantic_probe_not_linked", "pass"),
                                               (4.0, "schema_ok_semantic_probe_not_linked", "pass"),
                                               (3, "fail", "fail")])
def test_numeric_cache_schema_keeps_the_former_rule(schema, tok, prev):
    t, basis, p = gribs.semantic_status_for_backend(_Stub(name="cea_python", uses_cea=True, cache_schema=schema))
    assert (t, p) == (tok, prev) and "unknown" not in basis


@pytest.mark.integration
def test_reference_run_states_its_analytic_yg_range(tmp_path):
    for mode in ("single_phase_legacy", "homogeneous_equilibrium"):
        d = tmp_path / mode
        doc = _doc(d, "g1_reference_constant.json", rebuild=True)
        doc["two_phase_model"]["mode"] = mode
        d.mkdir(parents=True, exist_ok=True)
        assert _run(d, doc) == 0
        v = json.loads((d / "out" / "summary.json").read_text())["results"]["validation"]
        assert v["validated_Yg_range"] == [1.0, 1.0]
        assert v["table_Yg_range"] is None
        assert "analytic phase-fraction range [1, 1]" in v["validated_range_basis"]
        assert "defaults" not in v["validated_range_basis"]
