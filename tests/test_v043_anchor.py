"""External regression anchor: current GRIBS == published v0.4.3-alpha (total quantities).

The v0.4.3-alpha chamber EOS and nozzle flow are already on the total-product-mass basis
(p0*Vg = m*Theta, Theta = (Ru/M_CEA)*T0, mdot = Cd*At*p0*sqrt(gamma/(R T0))*f(gamma)), so for a condensed
AP/Al composition (Yg about 0.72) the corrected v0.5.0 model reproduces it exactly (independent review
2026-10-09, N-3; measured relative difference 0.0 for all total quantities).

The v0.4.3-alpha program is taken from the git tag v0.4.3-alpha of this repository (``git show``); the tests
are skipped when git, the tag or the ``cea`` package is not available.  The configuration is a v0.4.3-alpha
document; the current program migrates it automatically (the migration warning goes to stderr).
epsilon = 1 is used because v0.4.3-alpha stops its own checks for epsilon > 1 (known, fixed in v0.5.0).
"""
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
BODY = sorted(ROOT.glob("GRIBS_v*.py"))[0]
TAG = "v0.4.3-alpha"
TAG_FILE = "GRIBS_v0.4.3-alpha.py"

pytestmark = pytest.mark.integration

KEYS = ("burn_time_s", "p_max_Pa", "F_max_N", "impulse_burn_Ns", "impulse_total_Ns",
        "impulse_blowdown_Ns", "Isp_s", "cstar_eff_m_s", "CF_eff")


def v043_document(out_dir, cache_dir):
    return {
        "schema_version": "0.4.3-alpha",
        "propellant": {
            "packing_fraction": 0.95,
            "reactants": [
                {"name": "NH4CLO4(I)", "wt_percent": 70.0, "temperature_K": 298.15, "density_kg_m3": 1950.0},
                {"name": "AL(cr)", "wt_percent": 16.0, "temperature_K": 298.15, "density_kg_m3": 2700.0},
                {"name": "C7H8(L)", "wt_percent": 14.0, "temperature_K": 298.15, "density_kg_m3": 920.0}],
            "burn_law": {"coefficient_m_s": 3.0e-3, "pressure_exponent": 0.40,
                         "reference_pressure_Pa": 1.0e6, "temperature_sensitivity_1_K": 0.0,
                         "grain_temperature_K": 294.0, "reference_temperature_K": 294.0,
                         "erosive_burning": {"enabled": False, "alpha": 0.0, "beta": 53.0}}},
        "grain": {"initial_bore_radius_m": 0.008, "outer_radius_m": 0.030, "initial_length_m": 0.150,
                  "initial_free_volume_m3": 200e-6, "burning_end_faces": 1.0},
        "nozzle": {"initial_throat_radius_m": 0.004, "expansion_ratio": 1.0, "discharge_coefficient": 0.98,
                   "thrust_efficiency": 0.98, "flow_separation": {"enabled": True, "pressure_ratio": 0.4},
                   "throat_erosion": {"enabled": False, "rate_m_s_at_reference_pressure": 0.0,
                                      "pressure_exponent": 1.0}},
        "environment": {"ambient_pressure_Pa": 101325.0, "initial_chamber_pressure_Pa": 5.0e5},
        "igniter": {"mass_flow_kg_s": 0.02, "duration_s": 0.1},
        "thermochemistry": {
            "backend": "cea_python", "temperature_efficiency": 1.0, "outside_range_policy": "clamp",
            "pressure_range": {"minimum_Pa": 1.0e5, "maximum_Pa": 3.0e7}, "pressure_points": 81,
            "cea_python": {"ions": False, "transport": False, "trace": 1e-10, "products_from_reactants": True,
                           "product_species": [], "omit_species": [], "insert_species": [],
                           "molecular_weight": "gas_phase_M", "smooth_truncation": False,
                           "truncation_width": -1.0, "rocket_diagnostics": True, "cache_enabled": True,
                           "rebuild_cache": True, "cache_directory": str(cache_dir)}},
        "solver": {"method": "LSODA", "relative_tolerance": 1e-10, "pressure_absolute_tolerance_Pa": 0.001,
                   "burn_depth_absolute_tolerance_m": 1e-13, "burning_time_limit_s": 600.0,
                   "maximum_burning_step_s": 0.02, "unchoked_policy": "switch",
                   "blowdown": {"enabled": True, "time_limit_s": 5.0, "maximum_step_s": 0.0005}},
        "output": {"directory": str(out_dir), "figure_filename": "ballistics.png",
                   "history_filename": "time_history.csv", "summary_text_filename": "summary.txt",
                   "summary_json_filename": "summary.json"},
    }


@pytest.fixture(scope="module")
def v043_program(tmp_path_factory):
    pytest.importorskip("cea")
    try:
        p = subprocess.run(["git", "-C", str(ROOT), "show", f"{TAG}:{TAG_FILE}"],
                           capture_output=True, timeout=60)
    except (OSError, subprocess.TimeoutExpired):
        pytest.skip("git is not available")
    if p.returncode != 0:
        pytest.skip(f"git tag {TAG} is not available in this checkout")
    path = tmp_path_factory.mktemp("v043") / TAG_FILE
    path.write_bytes(p.stdout)
    return path


def run(program, work):
    work.mkdir(parents=True, exist_ok=True)
    cfg = work / "config.json"
    cfg.write_text(json.dumps(v043_document(work / "out", work / "cache")))
    p = subprocess.run([sys.executable, "-B", str(program), "--config", str(cfg)], capture_output=True,
                       text=True, timeout=300, cwd=str(work))
    assert p.returncode == 0, p.stderr[-500:]
    return json.loads((work / "out" / "summary.json").read_text())


def test_current_program_matches_published_v043_total_quantities(v043_program, tmp_path):
    old = run(v043_program, tmp_path / "v043")
    new = run(BODY, tmp_path / "current")
    for key in KEYS:
        a, b = old["results"][key], new["results"][key]
        assert b == pytest.approx(a, rel=1e-9), (key, a, b)
    # the anchor must be non-trivial: condensed products (Yg about 0.72) and a MEOP that the old
    # Yg-multiplied EOS / extra 1/Yg flow division would have lowered by about 40 %
    assert new["results"]["gas_mass_fraction_min"] < 0.75
    assert 9.0e6 < new["results"]["p_max_Pa"] < 9.6e6
