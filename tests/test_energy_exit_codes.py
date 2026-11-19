"""A-6 stage 2, C14: the exit codes of the failures outside the model, and the wall-clock limit (PRE_REGISTRATION_C14 and REV1, X-01..X-12 without X-07).

A TP cache that cannot be written is exit code 3 (E_THERMO_CACHE_WRITE_FAILED, wording "TP cache"); an unexpected exception is
code 4 with its traceback (E_ENERGY_UNEXPECTED_ERROR, the API raises it unchanged); an output that cannot be written is code 4
(E_ENERGY_OUTPUT_WRITE_FAILED), not the closure code 1; the run is limited in wall time over both phases.  The failing
conditions are real ones where they can be made (a path below a regular file, a directory without write permission, a counting
clock passed to the run); an unexpected exception is made by a stand-in `loaded` argument whose `.config` raises.  Nothing in
this file patches or assigns anything (guard table 186 section 8, row 3; PRE_REGISTRATION_C14_REV1).
"""
import inspect
import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path

import pytest

import gribs
import gribs.config.loader
import gribs.simulation.energy_integrator as ei
import gribs.simulation.energy_postprocess as ep
from gribs.energy.cea_tp import TPCache
from gribs.energy.provider import TPState
from gribs.errors import ConfigurationError, ThermochemistryCacheError, ThermochemistryError
from gribs.simulation.integrator import INTEGRATION_MAX_WALL_S, IntegrationBudgetError

ROOT = Path(__file__).resolve().parents[1]
G1 = ROOT / "cases" / "golden" / "g1_reference_constant.json"
G5 = ROOT / "cases" / "golden" / "g5_high_condensed_ap_al.json"
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
STATE = TPState(u_J_per_kg=-1.0e6, n_gas_per_kg=0.03, Yg=0.7, condensed_mass_fraction={"AL2O3(L)": 0.3}, gamma_s=1.2)
ID = {"test": "c14"}


def case_file(tmp_path, name, blowdown=True):
    doc = json.loads(G1.read_text(encoding="utf-8"))
    doc["output"]["directory"] = str(tmp_path / f"{name}_out")
    doc["solver"].update(relative_tolerance=1.0e-8, blowdown={"enabled": blowdown, "time_limit_s": 0.05, "maximum_step_s": 0.005})
    doc["two_phase_model"]["mode"] = "homogeneous_equilibrium"
    doc["energy_model"] = ENERGY
    path = tmp_path / f"{name}.json"
    path.write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")
    return path


def loaded_case(tmp_path, name, blowdown=True):
    return gribs.config.loader.load_json_configuration(case_file(tmp_path, name, blowdown))


def below_a_file(tmp_path, name="afile"):
    """A path whose parent is a regular file: nothing can be created there (a real condition, not an injected one)."""
    blocker = tmp_path / name
    blocker.write_text("not a directory\n", encoding="utf-8")
    return blocker / "below"


def read_only_directory(tmp_path, name="ro"):
    d = tmp_path / name
    d.mkdir()
    d.chmod(0o500)
    return d


needs_permissions = pytest.mark.skipif(os.geteuid() == 0, reason="not executable: uid 0 ignores directory permissions")


class CountingClock:
    """A clock that returns the number of times it was read before (0, 1, 2, ...)."""

    def __init__(self):
        self.reads = 0

    def __call__(self):
        value = float(self.reads)
        self.reads += 1
        return value


class RaisingLoaded:
    """A stand-in for the `loaded` argument of run_energy_case: reading `.config` raises `error` (the first thing the run does)."""

    def __init__(self, error):
        self._error = error

    @property
    def config(self):
        raise self._error


@pytest.mark.fast
def test_x01_the_tp_cache_creates_its_directory_before_writing(tmp_path):
    cache = TPCache(tmp_path / "a" / "b" / "c")
    cache.store(ID, 3000.0, 2.0e6, STATE)
    assert cache.stores == 1 and len(list((tmp_path / "a" / "b" / "c").iterdir())) == 1
    assert cache.load(ID, 3000.0, 2.0e6) == STATE


@pytest.mark.fast
def test_x02_a_tp_cache_directory_that_cannot_be_created_is_a_cache_write_error_with_the_tp_wording(tmp_path):
    cache = TPCache(below_a_file(tmp_path))
    with pytest.raises(ThermochemistryCacheError) as info:
        cache.store(ID, 3000.0, 2.0e6, STATE)
    e = info.value
    assert e.code == "E_THERMO_CACHE_WRITE_FAILED" and e.context["phase"] == "create_directory" and e.context["cache"] == "tp"
    assert "TP cache" in str(e) and "Chamber-property" not in str(e) and cache.stores == 0


@needs_permissions
@pytest.mark.fast
def test_x03_a_tp_cache_directory_without_write_permission_is_a_cache_write_error(tmp_path):
    d = read_only_directory(tmp_path)
    try:
        with pytest.raises(ThermochemistryCacheError) as info:
            TPCache(d).store(ID, 3000.0, 2.0e6, STATE)
    finally:
        d.chmod(0o700)
    e = info.value
    assert e.code == "E_THERMO_CACHE_WRITE_FAILED" and e.context["phase"] != "create_directory" and "TP cache" in str(e)


@pytest.mark.integration
def test_x04_an_energy_run_whose_tp_cache_cannot_be_written_ends_with_exit_code_3(tmp_path):
    pytest.importorskip("cea")
    doc = json.loads(G5.read_text(encoding="utf-8"))
    out = tmp_path / "g5_out"
    doc["output"]["directory"] = str(out)
    doc["solver"].update(burning_time_limit_s=0.002, relative_tolerance=1.0e-7)
    doc["thermochemistry"]["cea_python"]["cache_directory"] = str(below_a_file(tmp_path))
    doc["two_phase_model"]["mode"] = "homogeneous_equilibrium"
    doc["energy_model"] = {
        "initial_fill": {"temperature_K": 3000.0, "composition": "main_propellant_products", "element_tolerance": 1.0e-6},
        "wall": {"adiabatic": True},
        "burn_surface": {"form": "a", "grain_far_temperature_K": 293.0,
                         "solid_enthalpy": {"reference_J_kg": -1.4e6, "reference_temperature_K": 293.0, "cp_J_kgK": 1300.0}},
        "condensed": {"volume": "external_density", "latent_heat_check_tolerance_rel": 0.5,
                      "transitions": [{"solid_species": "AL2O3(a)", "liquid_species": "AL2O3(L)", "melting_temperature_K": 2327.0}],
                      "phases": [{"species": "AL2O3(a)", "phase": "solid", "density_kg_m3": 3990.0, "cp_J_kgK": 1300.0,
                                  "enthalpy_reference_J_kg": -1.6e7, "reference_temperature_K": 298.15},
                                 {"species": "AL2O3(L)", "phase": "liquid", "density_kg_m3": 3000.0, "cp_J_kgK": 1300.0,
                                  "enthalpy_reference_J_kg": -1.5e7, "reference_temperature_K": 298.15}]},
        "recovery": {"phase_fraction_tolerance": 1.0e-6, "root_tolerance": 1.0e-10, "plateau_edge_tolerance": 1.0e-6,
                     "temperature_search_K": [500.0, 5000.0], "start_temperatures_K": [3000.0]},
        "propellant_formulas": {"NH4CLO4(I)": {"N": 1, "H": 4, "Cl": 1, "O": 4}, "AL(cr)": {"Al": 1},
                                "C4H6,butadiene": {"C": 4, "H": 6}},
    }
    path = tmp_path / "g5_energy.json"
    path.write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")
    (entry,) = sorted(ROOT.glob("GRIBS_v*.py"))
    proc = subprocess.run([sys.executable, "-B", "-P", "-s", str(entry), "--config", str(path)], capture_output=True, text=True,
                          timeout=600, cwd=str(tmp_path))
    assert proc.returncode == 3, proc.stderr
    assert "E_THERMO_CACHE_WRITE_FAILED" in proc.stderr and "TP cache" in proc.stderr
    assert "E_ENERGY_" not in proc.stderr and not (out / ep.HISTORY_FILENAME).exists()


@pytest.mark.fast
def test_x05_an_exception_that_the_run_does_not_know_is_returned_as_unexpected():
    error = ZeroDivisionError("a stand-in: an exception that no code of the run raises on purpose")
    outcome = ep.run_energy_case(RaisingLoaded(error))
    assert outcome.unexpected is True and outcome.error is error


@pytest.mark.fast
def test_x06_an_unexpected_exception_prints_the_traceback_and_ends_with_exit_code_4(tmp_path, capsys):
    outcome = ep.run_energy_case(RaisingLoaded(ZeroDivisionError("a stand-in")))
    out = tmp_path / "x06_out"
    code = outcome.write_and_report(out, {}, {})
    err = capsys.readouterr().err
    assert code == 4
    assert "Traceback" in err and "ZeroDivisionError" in err and "E_ENERGY_UNEXPECTED_ERROR" in err
    assert not out.exists()


@pytest.mark.integration
def test_x08_an_output_directory_that_cannot_be_created_ends_with_exit_code_4_not_1(tmp_path, capsys):
    loaded = loaded_case(tmp_path, "x08")
    outcome = ep.run_energy_case(loaded)
    assert outcome.error is None and outcome.run_passed is True
    code = outcome.write_and_report(below_a_file(tmp_path), loaded.output_names, {})
    assert code == 4 and "E_ENERGY_OUTPUT_WRITE_FAILED" in capsys.readouterr().err


@needs_permissions
@pytest.mark.integration
def test_x09_an_output_directory_without_write_permission_ends_with_exit_code_4(tmp_path, capsys):
    loaded = loaded_case(tmp_path, "x09")
    outcome = ep.run_energy_case(loaded)
    d = read_only_directory(tmp_path)
    try:
        code = outcome.write_and_report(d, loaded.output_names, {})
        written = list(d.iterdir())
    finally:
        d.chmod(0o700)
    assert code == 4 and "E_ENERGY_OUTPUT_WRITE_FAILED" in capsys.readouterr().err and written == []


@pytest.mark.integration
def test_x10_the_wall_clock_limit_counts_both_phases_together(tmp_path):
    full = loaded_case(tmp_path, "x10_full").config
    burn_only = loaded_case(tmp_path, "x10_burn", blowdown=False).config
    inf = ei.WallBudget(limit_s=math.inf, clock=CountingClock(), check_every=1)
    c_all = ei.run_model_energy(full, wall_budget=inf)["integration_budget"]["rhs_calls"]
    c_burn = ei.run_model_energy(burn_only, wall_budget=inf)["integration_budget"]["rhs_calls"]
    assert c_all - c_burn > 10                                     # the blowdown phase is long enough to tell the two apart
    clock = CountingClock()
    with pytest.raises(IntegrationBudgetError) as info:
        ei.run_model_energy(full, wall_budget=ei.WallBudget(limit_s=float(c_burn + 5), clock=clock, check_every=1))
    assert "E_INTEGRATION_BUDGET_EXCEEDED" in str(info.value) and "wall time" in str(info.value)
    assert clock.reads == 1 + (c_burn + 6)                         # one start reading, then the (c_burn + 6)-th evaluation


@pytest.mark.fast
def test_x11_the_default_wall_budget_is_the_legacy_one():
    b = ei.WallBudget.default()
    assert b.limit_s == INTEGRATION_MAX_WALL_S == 300.0 and b.clock is time.monotonic and b.check_every == 1024
    assert inspect.signature(ei.run_model_energy).parameters["wall_budget"].default is None


@pytest.mark.fast
def test_x12_the_known_errors_keep_their_exit_codes_and_a_keyboard_interrupt_is_not_caught(tmp_path):
    for error, expected in ((ConfigurationError("x"), 2), (ThermochemistryError("x"), 3), (RuntimeError("x"), 4)):
        outcome = ep.EnergyOutcome(error=error)
        assert outcome.unexpected is False
        assert outcome.write_and_report(tmp_path / "unused", {}, {}) == expected
    with pytest.raises(KeyboardInterrupt):
        ep.run_energy_case(RaisingLoaded(KeyboardInterrupt()))
