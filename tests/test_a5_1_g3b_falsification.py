"""A-5.1 G3-b: falsification tests for the CEA mass-basis correction.

These tests are NOT a recovered red-before record of the historical bug.  They
show that independent CEA oracles detect deliberately re-injected versions of
the two errors corrected in A-5.1-4, applied to the current implementation:

  M1  chamber Psi with a second Yg factor   (Psi = Yg * R_eff * T)
  M2  nozzle total flow with an extra /Yg   (mdot_total = mdot_ref / Yg)

Oracles (computed from the official ``cea`` package, not from GRIBS tables):
  O1  Psi(p) == p / EqSolution.density(p)         (CEA 3.3.4, eta_T0 = 1)
  O2  GRIBS effective c* / CEA RocketSolution c*  (quasi-steady run)

Both mutations are invisible when Yg == 1, so a condensed-phase case is used:
the golden g4 input, Yg ~ 0.78.  Despite its file name (g4_medium_condensed),
g4 falls in the HIGH-condensation band of the 2026-10-07 definition (medium is
Yg 0.85-0.95); the file name is kept because it is a frozen golden input.
A single-phase case is checked to show that it cannot discriminate (design
note, not a pass criterion for M1/M2).
"""
import copy
import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np
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

cea = pytest.importorskip("cea")

G4 = ROOT / "cases" / "golden" / "g4_medium_condensed.json"
G2 = ROOT / "cases" / "golden" / "g2_single_phase_cea.json"

PSI_ORACLE_TOL = 1.0e-6     # interpolation at table nodes, eta_T0 = 1
CSTAR_ORACLE_TOL = 0.01     # quasi-steady transient vs CEA infinite-area c*
MUTATION_MIN_SHIFT = 0.10   # a re-injected Yg error must move the oracle >10 %


def _document(path, tmp_path):
    doc = json.loads(path.read_text(encoding="utf-8"))
    doc = copy.deepcopy(doc)
    doc["output"]["directory"] = str(tmp_path / "out")
    th = doc["thermochemistry"]
    th["cea_python"]["cache_directory"] = str(tmp_path / "cache")
    th["cea_python"]["rebuild_cache"] = True
    return doc


def _backend(path, tmp_path):
    conf = gribs.configuration_from_document(_document(path, tmp_path)).config
    assert conf.eta_T0 == 1.0
    return conf, gribs.initialize_thermochemistry(conf, tmp_path / "cache")


def _cea_density(conf, pressures_pa):
    names = [r.name for r in conf.reactants_spec_list]
    w = np.array([r.wt_percent for r in conf.reactants_spec_list], float)
    t = np.array([r.temperature_K for r in conf.reactants_spec_list], float)
    rmix = cea.Mixture(names)
    pmix = cea.Mixture(names, products_from_reactants=True, omit=[], ions=False)
    solver = cea.EqSolver(pmix, reactants=rmix, transport=False, ions=False,
                          trace=1e-10, insert=[])
    sol = cea.EqSolution(solver)
    hR = float(rmix.calc_property(cea.ENTHALPY, w, t)) / float(cea.R)
    out = []
    for p in pressures_pa:
        solver.solve(sol, cea.HP, hR, p / 1e5, w)
        assert sol.converged
        out.append(float(sol.density))
    return np.array(out)


def _psi_oracle_error(backend, conf):
    p_nodes = np.asarray(backend.table["pressure_Pa"], float)[::10]
    rho = _cea_density(conf, p_nodes)
    psi = np.array([backend.psi_and_derivative(float(p))[0] for p in p_nodes])
    return float(np.max(np.abs(psi / (p_nodes / rho) - 1.0)))


def _inject_m1(monkeypatch):
    original = gribs.PressureTableBackend.psi_and_derivative

    def psi_double_yg(self, p):
        psi, dpsi = original(self, p)
        yg, _ = self.phase_fractions(p)
        return yg * psi, yg * dpsi
    monkeypatch.setattr(gribs.PressureTableBackend, "psi_and_derivative", psi_double_yg)


def _inject_m2(monkeypatch):
    def split_extra_div_yg(self, mdot_total_reference, Yg_chamber):
        if self.mode == gribs.TWO_PHASE_SINGLE_LEGACY:
            return mdot_total_reference, mdot_total_reference, 0.0
        total = mdot_total_reference / Yg_chamber
        return total, Yg_chamber * total, total - Yg_chamber * total
    monkeypatch.setattr(gribs.NozzleTwoPhaseModel, "split", split_extra_div_yg)


def _run_cstar_ratio(path, tmp_path):
    doc = _document(path, tmp_path)
    cfg = tmp_path / "config.json"
    cfg.write_text(json.dumps(doc), encoding="utf-8")
    rc = gribs.main(["--config", str(cfg)])
    assert rc == 0
    summary = json.loads((tmp_path / "out" / "summary.json").read_text())
    res = summary["results"]
    cea_cstar = summary["results"]["cea_theoretical_diagnostics"]["cea_theoretical_cstar_m_s"] \
        if "cea_theoretical_diagnostics" in res else _find(summary, "cea_theoretical_cstar_m_s")
    return res["cstar_eff_m_s"] / cea_cstar, res["gas_mass_fraction_min"]


def _find(obj, key):
    if isinstance(obj, dict):
        if key in obj:
            return obj[key]
        for v in obj.values():
            r = _find(v, key)
            if r is not None:
                return r
    return None


@pytest.mark.integration
def test_o1_psi_matches_cea_density_current_implementation(tmp_path):
    conf, backend = _backend(G4, tmp_path)
    assert _psi_oracle_error(backend, conf) < PSI_ORACLE_TOL


@pytest.mark.integration
def test_o1_detects_reinjected_double_yg_in_psi(tmp_path, monkeypatch):
    conf, backend = _backend(G4, tmp_path)
    _inject_m1(monkeypatch)
    assert _psi_oracle_error(backend, conf) > MUTATION_MIN_SHIFT


@pytest.mark.integration
def test_o2_cstar_matches_cea_current_implementation(tmp_path):
    ratio, yg_min = _run_cstar_ratio(G4, tmp_path)
    assert yg_min < 0.9
    assert abs(ratio - 1.0) < CSTAR_ORACLE_TOL


@pytest.mark.integration
def test_o2_detects_reinjected_extra_div_yg_in_nozzle(tmp_path, monkeypatch):
    _inject_m2(monkeypatch)
    ratio, yg_min = _run_cstar_ratio(G4, tmp_path)
    assert abs(ratio - 1.0) > MUTATION_MIN_SHIFT
    # c*_mutated / c*_CEA should be ~Yg (flow overstated by 1/Yg)
    assert ratio == pytest.approx(yg_min, rel=0.05)


@pytest.mark.integration
def test_single_phase_case_cannot_discriminate_mutations(tmp_path, monkeypatch):
    """Design note: with Yg == 1 both mutations are identities."""
    _inject_m2(monkeypatch)
    _inject_m1(monkeypatch)
    ratio, yg_min = _run_cstar_ratio(G2, tmp_path)
    assert yg_min == pytest.approx(1.0, abs=1e-12)
    assert abs(ratio - 1.0) < CSTAR_ORACLE_TOL


@pytest.mark.integration
def test_internal_checks_detect_reinjected_double_yg_in_full_run(tmp_path, monkeypatch, capsys):
    """M1 in a full run: GRIBS' own post-run/ledger gates (not independent)."""
    _inject_m1(monkeypatch)
    doc = _document(G4, tmp_path)
    cfg = tmp_path / "config.json"
    cfg.write_text(json.dumps(doc), encoding="utf-8")
    rc = gribs.main(["--config", str(cfg)])
    err = capsys.readouterr().err
    assert rc != 0
    assert "W_POST_RUN_CONSISTENCY" in err
    assert "W_MASS_LEDGER_FAILED" in err
