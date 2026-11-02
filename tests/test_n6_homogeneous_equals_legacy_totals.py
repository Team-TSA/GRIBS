"""10/09 N-6 (L1 identity): after A-5.1-4 the homogeneous model equals
single_phase_legacy in total quantities.

Each mode generates its own CEA table (the cache key includes the two-phase
mode, so the two cache files have different keys).  The test checks the claim
that both runs use the same table data: the two stored cache files must be
identical apart from ``cache_key``.  Requires P2 (CEA backend +
single_phase_legacy recursion fix)."""
import copy, importlib.util, json, sys
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

TOTALS = ("p_max_Pa", "impulse_total_Ns", "Isp_s", "cstar_eff_m_s", "CF_eff",
          "F_max_N", "burn_time_s")


def _run(tmp_path, case, mode):
    doc = copy.deepcopy(json.loads((ROOT / "cases" / "golden" / case).read_text()))
    d = tmp_path / mode
    doc["output"]["directory"] = str(d / "out")
    doc["two_phase_model"]["mode"] = mode
    doc["thermochemistry"]["cea_python"]["cache_directory"] = str(d / "cache")
    d.mkdir(parents=True)
    cfg = d / "c.json"; cfg.write_text(json.dumps(doc))
    assert gribs.main(["--config", str(cfg)]) == 0
    return json.loads((d / "out" / "summary.json").read_text())["results"]


def _stored_table(tmp_path, mode):
    files = sorted((tmp_path / mode / "cache").glob("thermo_cea_python_*.json"))
    assert len(files) == 1, files
    table = json.loads(files[0].read_text())
    key = table.pop("cache_key")
    return key, json.dumps(table, sort_keys=True)   # NaN-safe comparison of the stored data


@pytest.mark.integration
@pytest.mark.parametrize("case", ["g4_medium_condensed.json", "g5_high_condensed_ap_al.json"])
def test_homogeneous_equals_legacy_in_totals(tmp_path, case):
    pytest.importorskip("cea")
    hom = _run(tmp_path, case, "homogeneous_equilibrium")
    leg = _run(tmp_path, case, "single_phase_legacy")
    # both modes really take their own paths
    assert hom["two_phase_model"]["mode"] == "homogeneous_equilibrium"
    assert hom["two_phase_model"]["active"] is True
    assert leg["two_phase_model"]["mode"] == "single_phase_legacy"
    assert leg["two_phase_model"]["active"] is False
    assert hom["gas_mass_fraction_min"] < 0.9                  # condensed phase present
    assert hom["maximum_condensed_mass_fraction_during_run"] > 0.1
    assert leg["maximum_condensed_mass_fraction_during_run"] == 0.0   # Yg := 1 path
    for k in TOTALS:
        assert hom[k] == pytest.approx(leg[k], rel=1e-9), k
    # the two runs used separately generated tables with identical stored data
    key_hom, data_hom = _stored_table(tmp_path, "homogeneous_equilibrium")
    key_leg, data_leg = _stored_table(tmp_path, "single_phase_legacy")
    assert key_hom != key_leg
    assert data_hom == data_leg
