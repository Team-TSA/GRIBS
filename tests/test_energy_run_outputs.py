"""A-6 stage 2, C12: the energy-mode output, the command line and the API; the legacy model unchanged (PRE_REGISTRATION_C12).

L-01..L-03: the legacy model is not changed by energy mode.  O-01..O-04: an energy configuration runs through the command line and
the API, with a dedicated minimal output (D-05 (a)) that says which legacy outputs it does not produce.  The case is the golden g1
with the analytic reference_constant backend, so the chemistry is synthetic (energy mode is verified on synthetic fixtures only);
a smoke run with the CEA provider of a golden propellant is in test_energy_cea_tp.py.  No monkeypatch, no assignment, no parametrize.
"""
import json
import subprocess
import sys
import types
from dataclasses import asdict, fields
from pathlib import Path

import pytest

import gribs
import gribs.config
import gribs.config.loader
import gribs.energy.wiring as ep_wiring
import gribs.simulation.energy_postprocess as ep
from gribs.config import Config
from gribs.errors import ConfigurationError, ThermochemistryError

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
G1 = ROOT / "cases" / "golden" / "g1_reference_constant.json"
FORMULAS = {"NH4CLO4(I)": {"N": 1, "H": 4, "Cl": 1, "O": 4}, "C4H6,butadiene": {"C": 4, "H": 6}}
ENERGY = {
    "initial_fill": {"temperature_K": 300.0, "composition": "main_propellant_products", "element_tolerance": 1.0e-6},
    "wall": {"adiabatic": True},
    "burn_surface": {"form": "a", "grain_far_temperature_K": 293.0,
                     "solid_enthalpy": {"reference_J_kg": -2.0e6, "reference_temperature_K": 293.0, "cp_J_kgK": 1300.0}},
    "condensed": {"volume": "external_density", "phases": [], "transitions": []},
    "recovery": {"phase_fraction_tolerance": 1.0e-6, "root_tolerance": 1.0e-10, "plateau_edge_tolerance": 1.0e-6,
                 "temperature_search_K": [200.0, 8000.0], "start_temperatures_K": [600.0, 3000.0]},
    "propellant_formulas": FORMULAS,
}


def case_file(tmp_path, name, energy=True):
    doc = json.loads(G1.read_text(encoding="utf-8"))
    doc["output"]["directory"] = str(tmp_path / f"{name}_out")
    doc["solver"].update(relative_tolerance=1.0e-8, blowdown={"enabled": True, "time_limit_s": 0.05, "maximum_step_s": 0.005})
    if energy:
        doc["two_phase_model"]["mode"] = "homogeneous_equilibrium"
        doc["energy_model"] = ENERGY
    path = tmp_path / f"{name}.json"
    path.write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")
    return path


def g5_energy_document(output_directory):
    """The golden g5 (backend cea_python, no reference_constant values) with an energy_model section; the cases of C12c only
    build the inputs, so no CEA call is made."""
    doc = json.loads((ROOT / "cases" / "golden" / "g5_high_condensed_ap_al.json").read_text(encoding="utf-8"))
    doc["output"]["directory"] = str(output_directory)
    doc["two_phase_model"]["mode"] = "homogeneous_equilibrium"
    doc["energy_model"] = {**ENERGY, "initial_fill": {**ENERGY["initial_fill"], "temperature_K": 3000.0},
                           "recovery": {**ENERGY["recovery"], "temperature_search_K": [500.0, 5000.0], "start_temperatures_K": [3000.0]},
                           "propellant_formulas": {"NH4CLO4(I)": FORMULAS["NH4CLO4(I)"], "AL(cr)": {"Al": 1},
                                                   "C4H6,butadiene": FORMULAS["C4H6,butadiene"]}}
    return doc


def entry():
    (body,) = sorted(ROOT.glob("GRIBS_v*.py"))
    return body


def run_cli(path, cwd):
    return subprocess.run([sys.executable, "-B", "-P", "-s", str(entry()), "--config", str(path)], capture_output=True, text=True,
                          timeout=600, cwd=str(cwd))


def child(code, *args):
    r = subprocess.run([sys.executable, "-B", "-P", "-s", "-c", code, str(SRC), *map(str, args)], capture_output=True, text=True,
                       timeout=600)
    assert r.returncode == 0, r.stderr[-3000:]
    return json.loads(r.stdout.strip().splitlines()[-1])


@pytest.fixture(scope="module")
def energy_cli(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("energy_case")
    path = case_file(tmp, "e1")
    proc = run_cli(path, tmp)
    return path, tmp / "e1_out", proc


@pytest.mark.integration
class TestLegacyUnchanged:
    def test_l01_a_legacy_run_does_not_import_any_energy_module(self, tmp_path):
        path = case_file(tmp_path, "l1", energy=False)
        code = ("import json, sys; sys.path.insert(0, sys.argv[1]); import gribs; res = gribs.run_simulation(sys.argv[2]); "
                "print(json.dumps([m for m in sys.modules if m.startswith('gribs.energy') or 'energy_' in m]))")
        assert child(code, path) == []

    def test_l02_a_legacy_command_line_run_writes_the_legacy_files_and_no_energy_keys(self, tmp_path):
        path = case_file(tmp_path, "l2", energy=False)
        proc = run_cli(path, tmp_path)
        assert proc.returncode == 0, proc.stderr[-2000:]
        out = tmp_path / "l2_out"
        summary = json.loads((out / "summary.json").read_text(encoding="utf-8"))
        assert not {"model_mode", "energy_model"} & set(summary)
        assert set(summary["output_files"]) == {"figure", "history_csv", "summary_text", "summary_json"}
        assert set(summary["resolved_configuration"]) == {f.name for f in fields(Config)}
        assert (out / "time_history.csv").is_file() and (out / "ballistics.png").is_file()
        assert not (out / ep.HISTORY_FILENAME).exists()

    def test_l03_the_condensed_volume_of_the_legacy_model_stays_neglected(self):
        doc = json.loads(G1.read_text(encoding="utf-8"))
        c = gribs.config.loader.configuration_from_document(doc).config
        assert c.two_phase_condensed_volume == "neglected"
        doc["two_phase_model"]["condensed_volume"] = "external_density"             # the only choice is "neglected"
        with pytest.raises(ConfigurationError, match="condensed_volume"):
            gribs.config.loader.configuration_from_document(doc)


@pytest.mark.integration
class TestEnergyRun:
    def test_o01_the_command_line_and_the_api_give_the_same_values(self, energy_cli):
        path, out, proc = energy_cli
        assert proc.returncode == 0, proc.stderr[-2000:]
        written = json.loads((out / "summary.json").read_text(encoding="utf-8"))
        res = gribs.run_simulation(path)
        assert res.run_passed is True and res.summary["model_mode"] == "energy"
        for key in ("configuration_source", "program", "derived_inputs", "two_phase_model", "thermochemistry",
                    "original_configuration", "resolved_configuration"):
            assert json.dumps(getattr(res, key), sort_keys=True) == json.dumps(written[key], sort_keys=True), key
        assert json.dumps(res.summary["energy_model"], sort_keys=True) == json.dumps(written["energy_model"], sort_keys=True)
        results = {k: v for k, v in res.summary.items() if k not in ("model_mode", "energy_model")}
        assert json.dumps(results, sort_keys=True) == json.dumps(written["results"], sort_keys=True)
        assert ep.history_csv_text(res.history) == (out / ep.HISTORY_FILENAME).read_text(encoding="utf-8")

    def test_o02_the_summary_and_the_command_line_say_what_energy_mode_does_not_produce(self, energy_cli):
        _path, out, proc = energy_cli
        s = json.loads((out / "summary.json").read_text(encoding="utf-8"))
        assert s["model_mode"] == "energy" and "energy" not in s["resolved_configuration"]
        em = s["energy_model"]
        assert em["burn_surface_interpretation"].startswith("A:")
        definition = em["H0_out_definition"]                  # H-18 closed (report 237): H0_out = h_cv is the stagnation enthalpy
        assert "H0_out = h_cv" in definition and "h_inlet + v_inlet^2/2 = h_cv" in definition and "steady adiabatic flow" in definition
        assert "H-05" in definition and "H-18" not in definition and "departs from the model definition" not in definition
        assert any(item.startswith("H-05") for item in em["held_items_adopted"])
        assert not any("H-18" in item for item in em["held_items_adopted"])
        assert "provisional (held item H-09)" in em["closure_tolerance_status"]
        assert "synthetic" in em["scope"] and em["gamma_s"]["min"] > 1.0 and em["gamma_s"]["max"] >= em["gamma_s"]["min"]
        assert em["provenance"]["provider"] == "ReferenceConstantProvider" and em["provenance"]["element_mass_fractions"]
        not_produced = " ".join(em["legacy_outputs_not_produced"])
        assert "figure" in not_produced and "diagnostics" in not_produced and "equilibrium pressure diagnostic" in not_produced
        assert "legacy figure and the legacy diagnostics are not produced" in proc.stdout     # the command line says so
        assert s["output_files"]["figure"] is None and not (out / "ballistics.png").exists()
        assert not (out / "time_history.csv").exists() and (out / ep.HISTORY_FILENAME).is_file()
        assert s["results"]["closure"]["passed"] is True and s["results"]["closure"]["tolerance_normalised"] == ep.CLOSURE_TOLERANCE

    def test_o03_the_diagnostic_the_backend_failure_and_the_legacy_flags_are_rejected_explicitly(self, tmp_path, capsys):
        c = gribs.config.loader.load_json_configuration(case_file(tmp_path, "o3")).config
        sub = type("ConfigWithEnergy", (Config,), {"energy": c.energy, "reactants_spec_list": c.reactants_spec_list})
        changed = sub(**{**asdict(c), "equilibrium_pressure_diagnostic": "sampled"})
        loaded = types.SimpleNamespace(config=changed, migration_notes=[], document={}, output_dir=tmp_path / "o3_out",
                                       output_names={})
        outcome = ep.run_energy_case(loaded)
        assert isinstance(outcome.error, ConfigurationError) and "equilibrium_pressure_diagnostic" in str(outcome.error)
        assert outcome.write_and_report(tmp_path / "o3_out", {}, {}) == 2 and not (tmp_path / "o3_out").exists()
        plain = ep.EnergyOutcome(error=ThermochemistryError("the cea package is missing"))        # as the legacy command line: 3
        assert plain.write_and_report(tmp_path / "t_out", {}, {}) == 3 and not (tmp_path / "t_out").exists()
        assert "Thermochemistry backend error" in capsys.readouterr().err

        class Coded(ThermochemistryError):
            code = "E_TEST_CODE"
        assert ep.EnergyOutcome(error=Coded("x")).write_and_report(tmp_path / "t_out", {}, {}) == 3
        assert "[E_TEST_CODE]" in capsys.readouterr().err
        assert ep.EnergyOutcome(error=RuntimeError("x")).write_and_report(tmp_path / "t_out", {}, {}) == 4   # other failures
        path = case_file(tmp_path, "o3b")
        for flag, extra in (("--selftest", []), ("--dump-thermo", [str(tmp_path / "t.json")]), ("--no-cache", []),
                            ("--allow-verification-failure", [])):
            proc = subprocess.run([sys.executable, "-B", "-P", "-s", str(entry()), "--config", str(path), flag, *extra],
                                  capture_output=True, text=True, timeout=600, cwd=str(tmp_path))
            assert proc.returncode == 2, (flag, proc.stderr[-500:])
            assert f"{flag} is not available for a configuration with an energy_model section" in proc.stderr, flag
            assert not (tmp_path / "o3b_out").exists() and not (tmp_path / "t.json").exists(), flag
        ok = subprocess.run([sys.executable, "-B", "-P", "-s", str(entry()), "--config", str(path), "--validate-config",
                             "--backend", "reference_constant"], capture_output=True, text=True, timeout=600, cwd=str(tmp_path))
        assert ok.returncode == 0, ok.stderr[-500:]                                            # --backend is not rejected

        # C12c (B-29 (d)): energy mode has no default values; the reference_constant values must be stated in the file
        g5 = tmp_path / "g5_energy.json"
        g5.write_text(json.dumps(g5_energy_document(tmp_path / "g5_out")), encoding="utf-8")
        proc = subprocess.run([sys.executable, "-B", "-P", "-s", str(entry()), "--config", str(g5), "--backend", "reference_constant"],
                              capture_output=True, text=True, timeout=600, cwd=str(tmp_path))
        assert proc.returncode == 2 and "energy mode has no default values" in proc.stderr, proc.stderr[-600:]
        assert all(k in proc.stderr for k in ("temperature_K", "gamma", "molecular_weight_kg_kmol"))
        assert not (tmp_path / "g5_out").exists()
        with pytest.raises(ConfigurationError, match="energy mode has no default values"):
            gribs.run_simulation(g5, backend_override="reference_constant")
        partial = g5_energy_document(tmp_path / "g5p_out")
        partial["thermochemistry"]["reference_constant"] = {"temperature_K": 3000.0, "gamma": 1.2}      # one value not stated
        (tmp_path / "g5_partial.json").write_text(json.dumps(partial), encoding="utf-8")
        with pytest.raises(ConfigurationError, match=r"\['molecular_weight_kg_kmol'\]"):
            ep_wiring.build_energy_inputs(gribs.config.loader.load_json_configuration(
                tmp_path / "g5_partial.json", backend_override="reference_constant").config)
        stated = g5_energy_document(tmp_path / "g5s_out")
        stated["thermochemistry"]["reference_constant"] = {"temperature_K": 2800.0, "gamma": 1.25, "molecular_weight_kg_kmol": 22.0}
        (tmp_path / "g5_stated.json").write_text(json.dumps(stated), encoding="utf-8")
        inputs = ep_wiring.build_energy_inputs(gribs.config.loader.load_json_configuration(
            tmp_path / "g5_stated.json", backend_override="reference_constant").config)       # all three stated: accepted
        assert (inputs.provider.gamma, inputs.provider.M, inputs.provider.T_ref) == (1.25, 22.0, 2800.0)

    def test_o04_the_time_history_has_the_columns_with_units_and_the_closure_residuals(self, energy_cli):
        _path, out, _proc = energy_cli
        lines = (out / ep.HISTORY_FILENAME).read_text(encoding="utf-8").splitlines()
        header = lines[0].split(",")
        assert tuple(header) == ep.HISTORY_COLUMNS
        suffixes = ("_s", "_kg", "_J", "_Pa", "_K", "_m", "_m3", "_m_s", "_kg_s", "_N", "_Ns", "_normalised", "gamma_s")
        assert all(h.endswith(suffixes) for h in header)
        rows = [[float(v) for v in line.split(",")] for line in lines[1:]]
        col = {h: [r[i] for r in rows] for i, h in enumerate(header)}
        s = json.loads((out / "summary.json").read_text(encoding="utf-8"))["results"]
        assert len(rows) == s["nodes"] > 20 and col["time_s"][0] == 0.0 and col["time_s"] == sorted(col["time_s"])
        assert max(col["time_s"]) > s["burning_end_time_s"]                                  # the blowdown nodes follow the burn
        assert max(abs(v) for h in header if h.startswith("R_") and h.endswith("_normalised") for v in col[h]) <= ep.CLOSURE_TOLERANCE
        assert all(g > 1.0 for g in col["gamma_s"]) and all(p > 0.0 for p in col["pressure_Pa"])
