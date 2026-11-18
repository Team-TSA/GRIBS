"""Wide-table variants of golden g2 and g5 (cases/wide_table/): runs that stay inside the property table.

Golden g2 (table 0.1-8 MPa, MEOP about 44.5 MPa) is a clamped-table regression reference with about 84 % of its
history above the table limit.  The variants change only ``pressure_range.maximum_Pa`` so that the whole run is
inside the table.  They are NOT frozen references: only structural facts are asserted, no recorded values.
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
BODY = sorted(ROOT.glob("GRIBS_v*.py"))[0]
spec = importlib.util.spec_from_file_location("gribs_wide_table", BODY)
g = importlib.util.module_from_spec(spec)
sys.modules["gribs_wide_table"] = g
spec.loader.exec_module(g)

pytestmark = pytest.mark.integration
WIDE = ROOT / "cases" / "wide_table"


def run_case(name, tmp_path):
    doc = json.loads((WIDE / f"{name}.json").read_text())
    doc["output"]["directory"] = str(tmp_path / "out")
    doc["thermochemistry"]["cea_python"]["cache_directory"] = str(tmp_path / "cache")
    c = g.configuration_from_document(doc).config
    g.initialize_thermochemistry(c, tmp_path / "cache")
    res = g.run_model(c)
    h = g.sample(res, c)
    return c, res, h, g.summarize(res, h, c)


@pytest.mark.parametrize("name, condensed", [
    ("g2_wide_table_single_phase_cea", False),
    ("g5_wide_table_high_condensed_ap_al", True),
])
def test_wide_table_case_stays_inside_the_table(name, condensed, tmp_path):
    c, res, h, summ = run_case(name, tmp_path)
    codes = [w["code"] for w in summ["warnings"]]
    assert not np.any(h["extrap"])
    assert summ["property_extrapolation"] is False and summ["property_extrapolation_fraction"] == 0.0
    assert "W_PROPERTY_TABLE_EXTRAPOLATED" not in codes
    assert "W_PROPERTY_TABLE_RANGE_CROSSED" not in codes
    assert res["integration_budget"]["table_restarts_burn"] == []
    assert summ["post_run_consistency"]["passed"] is True and summ["ledger"]["passed"] is True
    assert float(np.max(h["p0"])) < c.p_fit_max
    if condensed:
        assert summ["gas_mass_fraction_min"] < 0.75
    else:
        assert float(np.max(h["p0"])) > 8.0e6      # the case golden g2 clamps
