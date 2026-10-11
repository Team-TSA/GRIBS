"""A-8 golden freeze regression tests.

Each golden case is run from its config and the key performance metrics are
compared against the frozen golden_values.json.  Any physics change that
alters these values will cause a test failure, requiring an explicit golden
update with a documented reason.

Golden cases:
  g1  reference_constant   single-phase constant-property (analytic)
  g2  single_phase_cea     AP/HTPB, no aluminum, Yg = 1.0
  g3  low_condensed        NH4NO3/AL/C 88.7/0.7/10.6, Yg ~ 0.998
  g4  medium_condensed     AP/Al/HTPB 70/12/18, Yg ~ 0.78
  g5  high_condensed_ap_al AP/Al/HTPB 68/18/14, Yg ~ 0.70
  g6  cd_nozzle_blowdown   C-D nozzle (eps=2.5) + blowdown

model_status for g5 is ``approximate_with_limitations`` by design.
"""

import json
import subprocess
import sys
from pathlib import Path

import pytest

_GOLDEN_DIR = Path(__file__).resolve().parent.parent / "cases" / "golden"
_GOLDEN_VALUES = _GOLDEN_DIR / "golden_values.json"
_MAIN = Path(__file__).resolve().parent.parent / "GRIBS_v0.5.0-alpha.py"

# Relative tolerances for each metric group
_RTOL_TIGHT = 1e-9    # mass balance, timing
_RTOL_LOOSE = 1e-6    # performance (c*, Isp, CF)

# Metrics and their tolerances
_METRICS = {
    "p_max_Pa":                              _RTOL_LOOSE,
    "burn_time_s":                           _RTOL_TIGHT,
    "total_time_s":                          _RTOL_TIGHT,
    "impulse_total_Ns":                      _RTOL_LOOSE,
    "Isp_s":                                 _RTOL_LOOSE,
    "cstar_eff_m_s":                         _RTOL_LOOSE,
    "CF_eff":                                _RTOL_LOOSE,
    "F_max_N":                               _RTOL_LOOSE,
    "gas_mass_fraction_min":                 _RTOL_TIGHT,
    "gas_mass_fraction_max":                 _RTOL_TIGHT,
    "maximum_condensed_mass_fraction_during_run": _RTOL_TIGHT,
    "propellant_mass_balance_error":         _RTOL_TIGHT,
    "total_mass_consistency_error":          _RTOL_TIGHT,
}


def _golden_cases():
    """Yield (case_id, config_path) pairs from the golden directory."""
    if not _GOLDEN_VALUES.exists():
        pytest.skip("golden_values.json not found")
    data = json.loads(_GOLDEN_VALUES.read_text())
    for case_id in sorted(data):
        matches = list(_GOLDEN_DIR.glob(f"{case_id}_*.json"))
        if matches:
            yield pytest.param(case_id, matches[0], id=case_id)


@pytest.mark.parametrize("case_id,config_path", list(_golden_cases()))
def test_golden_regression(case_id, config_path):
    """Run a golden case and compare all metrics against frozen values."""
    golden = json.loads(_GOLDEN_VALUES.read_text())[case_id]

    # Run the simulation
    result = subprocess.run(
        [sys.executable, str(_MAIN), "--config", str(config_path)],
        capture_output=True, text=True, timeout=120,
    )
    assert result.returncode == 0, (
        f"Golden case {case_id} failed (exit {result.returncode}):\n"
        f"{result.stderr[-500:]}"
    )

    # Load the output summary
    repo_root = Path(__file__).resolve().parent.parent
    out_dir = repo_root / f"results_golden_{case_id}"
    summary_path = out_dir / "summary.json"
    assert summary_path.exists(), f"No summary.json for {case_id}"
    actual = json.loads(summary_path.read_text())["results"]

    # Compare each metric
    failures = []
    for metric, rtol in _METRICS.items():
        expected = golden.get(metric)
        got = actual.get(metric)
        if expected is None or got is None:
            continue
        if expected == 0.0 and got == 0.0:
            continue
        if expected == 0.0:
            rel_err = abs(got)
        else:
            rel_err = abs(got / expected - 1.0)
        if rel_err > rtol:
            failures.append(
                f"  {metric}: expected {expected:.6e}, got {got:.6e}, "
                f"rel_err = {rel_err:.2e} (tol {rtol:.0e})"
            )

    assert not failures, (
        f"Golden case {case_id} regression detected:\n" + "\n".join(failures)
    )
