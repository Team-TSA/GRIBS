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
  g7  low_condensed_graphite      as g3 with graphite C(gr) instead of gaseous C
  g8  cd_nozzle_blowdown_graphite as g6 with graphite C(gr) instead of gaseous C

g3 and g6 use the CEA species "C" (gaseous atomic carbon); they are kept
unchanged as numerical regression references and are not representative of a
real propellant.  g7 and g8 were added in v0.6.1-alpha (GitHub issue #7) with
the carbon density 2200 kg/m3 of g3 and g6 (the default configuration uses
2260).  Like g1-g6 they were frozen before the release gate G5, with the
maintainer's approval, as regression references that are not validated.

Each case is checked on 15 fields of golden_values.json: the 13 numeric
metrics below and the two-phase model mode and activity (exact match).

model_status for g5 is ``approximate_with_limitations`` by design.
"""

import json
import subprocess
import sys
from pathlib import Path

import pytest

_GOLDEN_DIR = Path(__file__).resolve().parent.parent / "cases" / "golden"
_GOLDEN_VALUES = _GOLDEN_DIR / "golden_values.json"
_MAIN = Path(__file__).resolve().parent.parent / "GRIBS_v0.6.2-alpha.py"

# Relative tolerances for each metric group
_RTOL_TIGHT = 1e-9    # mass balance, timing
_RTOL_LOOSE = 1e-6    # performance (c*, Isp, CF)

# Metrics and their tolerances
# Residual metrics (propellant mass balance, total mass consistency, maximum condensed mass fraction) are numerical
# residuals at round-off level.  A last-digit difference of the CEA table (another CPU, for example) or a different
# integration path changes them by tens of percent while the physical results move by <= ~1e-9 (reports 72, 76, 78),
# so a relative comparison with the frozen value is not meaningful.  They are judged by absolute bounds equal to the
# program's own verification thresholds:
#   propellant_mass_balance_error  <= max(MASS_LEDGER_PROPELLANT_TOLERANCE_FLOOR = 1e-9,
#                                         MASS_LEDGER_PROPELLANT_TOLERANCE_RTL_FACTOR = 100 x solver rtol)
#                                     (the mass-ledger propellant gate; the golden cases use rtol 1e-10 -> 1e-8)
#   total_mass_consistency_error   <= POST_RUN_MASS_TOLERANCE = 1e-6
#                                     (tolerance of the post-run mass-consistency check; this field is the retained
#                                      instantaneous form, the check itself uses the invariant form with the same bound)
#   maximum_condensed_mass_fraction_during_run: |got - expected| <= max(1e-9 x |expected|, 1e-12)
#                                     (1e-12 = _PHASE_FRACTION_ROUNDING_TOL; meaningful also for single-phase ~0)
# The physical and timing metrics keep the relative comparison with the frozen golden values.
_PROPELLANT_GATE_FLOOR = 1.0e-9
_PROPELLANT_GATE_RTOL_FACTOR = 100.0
_POST_RUN_MASS_TOLERANCE = 1.0e-6
_PHASE_FRACTION_ROUNDING_TOL = 1.0e-12
_RESIDUAL_METRICS = ("propellant_mass_balance_error", "total_mass_consistency_error",
                     "maximum_condensed_mass_fraction_during_run")


def _residual_check(metric, got, expected, solver_rtol):
    """Return (checked value, bound) for a residual metric."""
    if metric == "propellant_mass_balance_error":
        return got, max(_PROPELLANT_GATE_FLOOR, _PROPELLANT_GATE_RTOL_FACTOR * solver_rtol)
    if metric == "total_mass_consistency_error":
        return got, _POST_RUN_MASS_TOLERANCE
    return abs(got - expected), max(1.0e-9 * abs(expected), _PHASE_FRACTION_ROUNDING_TOL)


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
def test_golden_regression(case_id, config_path, tmp_path):
    """Run a golden case and compare all metrics against frozen values.

    The run output goes to pytest's tmp_path.  For the cea_python cases
    (g2-g8) the CEA cache directory is set explicitly to
    tmp_path/results_golden_<id>/thermo_cache (the same place the default
    rule would choose, made independent of it); the reference_constant case
    (g1) has no cache.  The frozen inputs (cases/golden/*.json,
    golden_values.json) and the tracked results_golden_* reference
    directories are never written.  The values of g1-g6 are unchanged;
    g7 and g8 were added in v0.6.1-alpha.
    """
    golden = json.loads(_GOLDEN_VALUES.read_text())[case_id]

    doc = json.loads(Path(config_path).read_text())
    out_dir = tmp_path / f"results_golden_{case_id}"
    doc["output"]["directory"] = str(out_dir)
    if "cea_python" in doc.get("thermochemistry", {}):
        doc["thermochemistry"]["cea_python"]["cache_directory"] = str(out_dir / "thermo_cache")
    run_config = tmp_path / Path(config_path).name
    run_config.write_text(json.dumps(doc))

    # Run the simulation
    result = subprocess.run(
        [sys.executable, str(_MAIN), "--config", str(run_config)],
        capture_output=True, text=True, timeout=120,
    )
    assert result.returncode == 0, (
        f"Golden case {case_id} failed (exit {result.returncode}):\n"
        f"{result.stderr[-500:]}"
    )

    # Load the output summary
    summary_path = out_dir / "summary.json"
    assert summary_path.exists(), f"No summary.json for {case_id}"
    actual = json.loads(summary_path.read_text())["results"]

    # Compare each metric
    failures = []
    for key, got in (("two_phase_mode", actual["two_phase_model"]["mode"]),
                     ("two_phase_active", actual["two_phase_model"]["active"])):
        if got != golden[key]:
            failures.append(f"  {key}: expected {golden[key]!r}, got {got!r}")
    for metric, rtol in _METRICS.items():
        expected = golden.get(metric)
        got = actual.get(metric)
        if expected is None or got is None:
            continue
        if metric in _RESIDUAL_METRICS:
            value, bound = _residual_check(metric, got, expected, float(doc["solver"]["relative_tolerance"]))
            if not value <= bound:
                failures.append(f"  {metric}: got {got:.6e} (checked {value:.3e}) > bound {bound:.3e}")
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
