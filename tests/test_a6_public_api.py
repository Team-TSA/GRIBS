"""A-6 unit 10: the public API of the package (U-39 (a); ADR_A6_PACKAGE_SPLIT, report 114 addendum).

gribs exposes run_simulation, Result and __version__; `import gribs` does not import matplotlib (nor numpy / scipy),
the public names are the objects of their owner modules, and run_simulation returns the values that the command
line writes to summary.json and time_history.csv for the same configuration file.
"""
import ast
import importlib
import json
import subprocess
import sys
import types
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
PUBLIC = ["run_simulation", "Result", "__version__"]
#: Result field -> key of summary.json that the command line writes from the same value
SUMMARY_KEYS = {"configuration_source": "configuration_source", "configuration_migration": "configuration_migration",
                "original_configuration": "original_configuration",
                "resolved_configuration": "resolved_configuration", "derived_inputs": "derived_inputs",
                "two_phase_model": "two_phase_model", "thermochemistry": "thermochemistry", "summary": "results",
                "program": "program"}


def _entry():
    bodies = sorted(ROOT.glob("GRIBS_v*.py"))
    if len(bodies) != 1:
        raise RuntimeError(f"expected exactly one GRIBS_v*.py, found {[b.name for b in bodies]}")
    return bodies[0]


def _child(code, *args, timeout=120):
    """Run `code` in a fresh interpreter with src/ first on sys.path; its last stdout line is JSON."""
    r = subprocess.run([sys.executable, "-B", "-P", "-s", "-c", code, str(SRC), *map(str, args)],
                       capture_output=True, text=True, timeout=timeout)
    assert r.returncode == 0, r.stderr[-3000:]
    return json.loads(r.stdout.strip().splitlines()[-1])


def _case(tmp_path, source, name):
    """A copy of the configuration `source` whose output directory is inside tmp_path."""
    doc = json.loads(Path(source).read_text(encoding="utf-8"))
    doc["output"]["directory"] = str(tmp_path / f"{name}_out")
    path = tmp_path / f"{name}.json"
    path.write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")
    return path


def _text(value):
    return json.dumps(value, indent=2, ensure_ascii=False)


_LOADED = ("import json, sys; sys.path.insert(0, sys.argv[1]); import gribs; "
           "loaded = lambda: sorted(m for m in sys.modules "
           "if m.split('.')[0] in ('matplotlib', 'numpy', 'scipy', 'gribs')); ")


def test_import_gribs_does_not_import_matplotlib():
    out = _child(_LOADED + "print(json.dumps(loaded()))")
    assert out == ["gribs", "gribs._version"], out


def test_resolving_the_public_names_does_not_import_matplotlib():
    """Reaching run_simulation imports gribs.api and the simulation modules, but not gribs.cli / gribs.plotting."""
    out = _child(_LOADED + "gribs.run_simulation; gribs.Result; print(json.dumps(loaded()))")
    assert not [m for m in out if m.split(".")[0] == "matplotlib"], out
    assert "gribs.api" in out and "gribs.simulation.result" in out
    assert not {"gribs.cli", "gribs.plotting", "gribs.selftest"} & set(out), out


@pytest.mark.fast
def test_public_names_are_the_objects_of_their_owners():
    gribs = importlib.import_module("gribs")
    api = importlib.import_module("gribs.api")
    result = importlib.import_module("gribs.simulation.result")
    version = importlib.import_module("gribs._version")
    constants = importlib.import_module("gribs.constants")
    assert gribs.__all__ == PUBLIC
    assert gribs.run_simulation is api.run_simulation
    assert gribs.Result is result.Result and gribs.Result.__module__ == "gribs.simulation.result"
    assert gribs.__version__ is version.PROGRAM_VERSION is constants.PROGRAM_VERSION
    assert set(PUBLIC) <= set(dir(gribs))
    for name in ("no_such_name", "main", "run_model"):  # the command line and internals are not public names
        with pytest.raises(AttributeError):
            getattr(gribs, name)


@pytest.mark.fast
def test_api_blocks_are_the_text_of_main():
    """The blocks of gribs.api between the "same text as gribs.cli.main" markers are found, unchanged, in main."""
    cli_text = (SRC / "gribs" / "cli.py").read_text()
    main = next(n for n in ast.parse(cli_text).body if isinstance(n, ast.FunctionDef) and n.name == "main")
    main_text = ast.get_source_segment(cli_text, main, padded=True)
    blocks, block = [], None
    for line in (SRC / "gribs" / "api.py").read_text().splitlines(keepends=True):
        if line.strip() == "# --- begin: same text as gribs.cli.main ---":
            block = []
        elif line.strip() == "# --- end: same text as gribs.cli.main ---":
            blocks.append("".join(block)); block = None
        elif block is not None:
            block.append(line)
    assert len(blocks) == 4 and all(b.strip() for b in blocks), blocks      # C12: + the energy-mode branch
    for b in blocks:
        assert b in main_text, b


@pytest.mark.fast
def test_run_simulation_raises_applicability_error_before_the_solver(tmp_path, monkeypatch):
    """As main (exit 1, nothing written), run_simulation stops before the solver, with the code and context of main's
    record; the backend state of gribs.thermo is restored to the one installed before the call."""
    api = importlib.import_module("gribs.api"); thermo = importlib.import_module("gribs.thermo")
    config = _case(tmp_path, ROOT / "gribs_config.json", "out_of_scope")
    backend = types.SimpleNamespace(  # built with its values, so the test has no attribute assignment (guard T0)
        supports_two_phase=True, property_validation={"gas_mass_fraction_min": 0.65, "gas_mass_fraction_max": 1.0},
        pmin=1.0e5, pmax=8.0e6)

    def install(_config, _cache_dir):
        thermo._set_thermo(backend)
        return backend

    monkeypatch.setattr("gribs.thermo.initialize_thermochemistry", install)
    monkeypatch.setattr("gribs.simulation.integrator.run_model",
                        lambda *_args: pytest.fail("out-of-scope run reached the solver"))
    marker = object()
    thermo._set_thermo(marker)  # restored by the conftest state guard
    with pytest.raises(api.ModelApplicabilityError) as info:
        api.run_simulation(config)
    assert info.value.code == "E_OUT_OF_VALIDATED_RANGE"
    assert info.value.context["simulation_started"] is False
    assert thermo._thermo_state() is marker
    assert not (tmp_path / "out_of_scope_out").exists()


_RUN = ("import json, sys; from pathlib import Path; sys.path.insert(0, sys.argv[1]); import gribs; "
        "import gribs.reporting, gribs.thermo; ")


@pytest.mark.integration
def test_run_simulation_returns_the_values_the_command_line_writes(tmp_path):
    config = _case(tmp_path, ROOT / "cases" / "golden" / "g1_reference_constant.json", "g1")
    r = subprocess.run([sys.executable, "-B", "-P", "-s", str(_entry()), "--config", str(config)],
                       capture_output=True, text=True, timeout=600, cwd=str(tmp_path))
    assert r.returncode == 0, r.stderr[-3000:]
    written = json.loads((tmp_path / "g1_out" / "summary.json").read_text(encoding="utf-8"))
    code = _RUN + (
        "res = gribs.run_simulation(sys.argv[2]); "
        "gribs.reporting.write_csv(res.history, Path(sys.argv[3])); "
        f"out = {{f: getattr(res, f) for f in {sorted(SUMMARY_KEYS)!r}}}; "
        "out.update(run_passed=res.run_passed, warnings_is_summary=res.warnings is res.summary['warnings'], "
        "type=type(res).__module__ + '.' + type(res).__name__, state=gribs.thermo._thermo_state() is None, "
        "matplotlib=sorted(m for m in sys.modules if m.split('.')[0] == 'matplotlib')); "
        "print(json.dumps(out))")
    out = _child(code, config, tmp_path / "api_history.csv", timeout=600)
    for field, key in SUMMARY_KEYS.items():
        assert _text(out[field]) == _text(written[key]), field
    assert (tmp_path / "api_history.csv").read_bytes() == (tmp_path / "g1_out" / "time_history.csv").read_bytes()
    assert out["run_passed"] is True and out["warnings_is_summary"] is True
    assert out["type"] == "gribs.simulation.result.Result"
    assert out["state"] is True and out["matplotlib"] == []


@pytest.mark.integration
def test_run_simulation_runs_two_cases_in_one_process(tmp_path):
    """Case A, case B, case A again in one process: the two runs of A give the same values, B differs, and the
    backend state is empty after each call."""
    a = _case(tmp_path, ROOT / "cases" / "golden" / "g1_reference_constant.json", "a")
    b = _case(tmp_path, ROOT / "cases" / "reference" / "level1_end_faces_1.json", "b")
    code = _RUN + (
        "out = []\n"
        "for path in sys.argv[2:]:\n"
        "    res = gribs.run_simulation(path)\n"
        "    out.append([json.dumps(res.summary), json.dumps({k: list(map(str, v)) for k, v in res.history.items()}), "
        "gribs.thermo._thermo_state() is None])\n"
        "print(json.dumps(out))")
    (a1, h1, s1), (b1, _hb, s2), (a2, h2, s3) = _child(code, a, b, a, timeout=600)
    assert s1 and s2 and s3
    assert a1 == a2 and h1 == h2
    assert a1 != b1
