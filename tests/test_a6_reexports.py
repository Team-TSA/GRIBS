"""A-6: every top-level name of the moved low-level modules is re-exported by the entry, and the entry and the
package hold the same objects (ADR_A6_PACKAGE_SPLIT, unit 1)."""
import ast
import importlib
import importlib.util
import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
MODULES = ("gribs.constants", "gribs.errors", "gribs.util", "gribs.config", "gribs.thermo",
           "gribs.propellant.burning", "gribs.nozzle", "gribs.chamber.equilibrium", "gribs.simulation.integrator",
           "gribs.simulation.postprocess", "gribs.plotting", "gribs.reporting", "gribs.selftest", "gribs.config.loader",
           "gribs.cli")


def _entry_path():
    bodies = sorted(ROOT.glob("GRIBS_v*.py"))
    if len(bodies) != 1:
        raise RuntimeError(f"expected exactly one GRIBS_v*.py, found {[b.name for b in bodies]}")
    return bodies[0]


def _module_file(module):
    """Source file of a package module: <name>.py, or <name>/__init__.py for a package (gribs.config since unit 9)."""
    base = ROOT / "src" / Path(*module.split("."))
    return base / "__init__.py" if base.is_dir() else base.with_suffix(".py")


def _top_level_names(module):
    tree = ast.parse(_module_file(module).read_text())
    out = set()
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.ClassDef)):
            out.add(node.name)
        elif isinstance(node, (ast.Assign, ast.AnnAssign)):
            for t in (node.targets if isinstance(node, ast.Assign) else [node.target]):
                for x in (t.elts if isinstance(t, ast.Tuple) else [t]):  # a tuple target binds each element
                    if isinstance(x, ast.Name):
                        out.add(x.id)
    return out


def _state_names(module):
    """Module state: names the module rebinds with `global`. They are never re-exported (a from-import would copy the
    value once and go stale) and are reached through their accessors only (ADR decision 2, guard T10)."""
    tree = ast.parse(_module_file(module).read_text())
    return {name for node in ast.walk(tree) if isinstance(node, ast.Global) for name in node.names}


def _entry_imports():
    out = {}
    for node in ast.parse(_entry_path().read_text()).body:
        if isinstance(node, ast.ImportFrom) and node.level == 0 and node.module in MODULES:
            out.setdefault(node.module, set()).update(a.asname or a.name for a in node.names)
    return out


def _load_entry():
    if "gribs_entry" in sys.modules:
        return sys.modules["gribs_entry"]
    spec = importlib.util.spec_from_file_location("gribs_entry", _entry_path())
    module = importlib.util.module_from_spec(spec); sys.modules["gribs_entry"] = module
    spec.loader.exec_module(module)
    return module


@pytest.mark.fast
def test_entry_reexports_every_top_level_name():
    imports = _entry_imports()
    missing = {m: sorted(_top_level_names(m) - _state_names(m) - imports.get(m, set())) for m in MODULES}
    assert not any(missing.values()), missing
    assert _state_names("gribs.thermo"), "gribs.thermo keeps the backend state"
    extra = {m: sorted(imports.get(m, set()) & _state_names(m)) for m in MODULES}
    assert not any(extra.values()), extra


@pytest.mark.fast
def test_entry_and_package_share_the_exception_classes():
    entry = _load_entry(); errors = importlib.import_module("gribs.errors")
    names = sorted(_top_level_names("gribs.errors"))
    assert len(names) == 7, names
    for name in names:
        assert getattr(entry, name) is getattr(errors, name), name


@pytest.mark.fast
def test_entry_and_package_agree_on_representative_values():
    entry = _load_entry(); constants = importlib.import_module("gribs.constants"); util = importlib.import_module("gribs.util")
    for name in ("PROGRAM_NAME", "PROGRAM_VERSION", "SCHEMA_VERSION", "KNOWN_BACKENDS", "R_UNIVERSAL",
                 "MODEL_STATUS_VALIDATED_YG_THRESHOLD", "MODEL_STATUS_OUT_OF_SCOPE_YG_THRESHOLD",
                 "HIGH_CONDENSATE_G4B_UNVALIDATED_YG_THRESHOLD", "POST_RUN_MASS_TOLERANCE", "G0"):
        assert getattr(entry, name) == getattr(constants, name), name
        assert getattr(entry, name) is getattr(constants, name), name
    assert entry._report_value_text is util._report_value_text


@pytest.mark.fast
def test_entry_and_package_share_the_configuration_objects():
    entry = _load_entry(); config = importlib.import_module("gribs.config")
    for name in sorted(_top_level_names("gribs.config")):
        assert getattr(entry, name) is getattr(config, name), name
    assert entry.Config.__module__ == "gribs.config"


@pytest.mark.fast
def test_entry_and_package_share_the_thermochemistry_objects_and_state():
    """Unit 3b: the entry re-exports the objects of gribs.thermo; the backend state exists once, in gribs.thermo."""
    entry = _load_entry(); thermo = importlib.import_module("gribs.thermo")
    for name in sorted(_top_level_names("gribs.thermo") - _state_names("gribs.thermo")):
        assert getattr(entry, name) is getattr(thermo, name), name
    for name in _state_names("gribs.thermo"):
        assert not hasattr(entry, name), name  # one state, in the package
    assert entry.PressureTableBackend.__module__ == "gribs.thermo"
    marker = object()
    entry._set_thermo(marker)  # restored by the conftest state guard
    assert thermo._thermo_state() is marker and entry.thermochemistry() is marker


@pytest.mark.fast
def test_entry_and_package_share_the_burn_rate_nozzle_and_equilibrium_objects():
    """Unit 4: the entry re-exports the objects of gribs.propellant.burning, gribs.nozzle and gribs.chamber.equilibrium."""
    entry = _load_entry()
    for module in ("gribs.propellant.burning", "gribs.nozzle", "gribs.chamber.equilibrium"):
        mod = importlib.import_module(module)
        names = sorted(_top_level_names(module))
        assert names, module
        for name in names:
            assert getattr(entry, name) is getattr(mod, name), (module, name)
    assert entry.NozzleTwoPhaseModel.__module__ == "gribs.nozzle"


@pytest.mark.fast
def test_top_level_names_count_tuple_targets():
    """The state indices are bound by one tuple assignment; the completeness check must see each of them."""
    names = _top_level_names("gribs.simulation.integrator")
    assert {"IP", "IX", "IRT", "IMO", "IIMP", "IMG", "IMIGN"} <= names, sorted(names)


@pytest.mark.fast
def test_entry_and_package_share_the_integrator_objects():
    """Unit 5: the entry re-exports the objects of gribs.simulation.integrator. sample() (gribs.simulation.postprocess)
    and the self-tests (gribs.selftest) import the state indices and igniter_mdot from gribs.simulation.integrator, so
    the entry, the integrator and its readers share the same objects."""
    entry = _load_entry(); integrator = importlib.import_module("gribs.simulation.integrator")
    for name in sorted(_top_level_names("gribs.simulation.integrator")):
        assert getattr(entry, name) is getattr(integrator, name), name
    assert (entry.IP, entry.IX, entry.IRT, entry.IMO, entry.IIMP, entry.IMG, entry.IMIGN) == tuple(range(7))
    assert entry.igniter_mdot is integrator.igniter_mdot
    assert entry.IntegrationBudgetError.__module__ == "gribs.simulation.integrator"


@pytest.mark.fast
def test_entry_and_package_share_the_postprocess_objects():
    """Unit 6: the entry re-exports the objects of gribs.simulation.postprocess; the post-processing reads the state
    indices and igniter_mdot of gribs.simulation.integrator (U-57, option A), so there is one definition of each."""
    entry = _load_entry(); postprocess = importlib.import_module("gribs.simulation.postprocess")
    integrator = importlib.import_module("gribs.simulation.integrator")
    names = sorted(_top_level_names("gribs.simulation.postprocess"))
    assert {"_phase_grid", "sample", "summarize", "mass_energy_ledger", "model_applicability_error",
            "nozzle_transition_diagnostics"} <= set(names), names
    for name in names:
        assert getattr(entry, name) is getattr(postprocess, name), name
    for name in ("IP", "IX", "IRT", "IMO", "IIMP", "IMG", "IMIGN", "igniter_mdot"):
        assert getattr(postprocess, name) is getattr(integrator, name), name
    assert entry.PhaseFractionValidationError.__module__ == "gribs.simulation.postprocess"


@pytest.mark.fast
def test_entry_and_package_share_the_plotting_and_reporting_objects():
    """Unit 7: the entry re-exports the objects of gribs.plotting and gribs.reporting; the figure reads the CEA
    diagnostics status of gribs.reporting (U-35 placement), so there is one definition of each."""
    entry = _load_entry()
    plotting = importlib.import_module("gribs.plotting"); reporting = importlib.import_module("gribs.reporting")
    for module, patched in ((plotting, {"make_figure"}), (reporting, {"summary_text", "write_csv"})):
        names = sorted(_top_level_names(module.__name__))
        assert patched <= set(names), names
        for name in names:
            assert getattr(entry, name) is getattr(module, name), name
    assert plotting.cea_diagnostics_status is reporting.cea_diagnostics_status
    assert entry.make_figure.__module__ == "gribs.plotting"
    assert entry.cea_diagnostics_status.__module__ == "gribs.reporting"


def test_plotting_sets_the_agg_backend_when_imported_alone():
    """gribs.plotting sets the Agg backend before pyplot, as the entry does, so the figure does not depend on whether
    the entry or the package imports pyplot first (report 127). MPLBACKEND asks for another backend, so the test
    fails if the module leaves the choice to matplotlib (which picks Agg by itself on a machine without a display)."""
    code = ("import sys; sys.path.insert(0, sys.argv[1]); import gribs.plotting, matplotlib; "
            "print(matplotlib.get_backend().lower())")
    env = dict(os.environ, MPLBACKEND="svg")
    out = subprocess.run([sys.executable, "-c", code, str(ROOT / "src")], capture_output=True, text=True, env=env,
                         check=True).stdout.split()
    assert out[-1] == "agg", out


@pytest.mark.fast
def test_entry_and_package_share_the_selftest_objects():
    """Unit 8: the entry re-exports the objects of gribs.selftest. The self-tests read the non-moved names of their
    owners (from-imported) and the moved classes of gribs.thermo as gribs.thermo.<name>, so each exists once."""
    entry = _load_entry(); selftest = importlib.import_module("gribs.selftest")
    thermo = importlib.import_module("gribs.thermo"); integrator = importlib.import_module("gribs.simulation.integrator")
    postprocess = importlib.import_module("gribs.simulation.postprocess")
    names = sorted(_top_level_names("gribs.selftest"))
    assert {"self_tests", "self_tests_post_run", "_SyntheticTestBackend"} <= set(names), names
    for name in names:
        assert getattr(entry, name) is getattr(selftest, name), name
    assert entry.self_tests_post_run.__module__ == "gribs.selftest"
    assert selftest._SyntheticTestBackend.__bases__ == (thermo.PressureTableBackend,)
    for name, owner in (("IP", integrator), ("make_rhs", integrator), ("theta_and_deriv", thermo),
                        ("post_run_consistency_report", postprocess)):
        assert getattr(selftest, name) is getattr(owner, name), name


@pytest.mark.fast
def test_entry_and_package_share_the_configuration_loader_and_cli_objects():
    """Unit 9: the entry re-exports the objects of gribs.config.loader and gribs.cli; the CLI reads the patched
    load_json_configuration through its owner (ADR decision 1), so it does not hold a copy of it."""
    entry = _load_entry()
    loader = importlib.import_module("gribs.config.loader"); cli = importlib.import_module("gribs.cli")
    for module, patched in ((loader, {"load_json_configuration"}),
                            (cli, {"_print_backend_banner", "environment_info", "_cea_theoretical_diagnostics"})):
        names = sorted(_top_level_names(module.__name__))
        assert patched <= set(names), names
        for name in names:
            assert getattr(entry, name) is getattr(module, name), name
    assert not hasattr(cli, "load_json_configuration")
    assert entry.main.__module__ == "gribs.cli" and entry.Configuration.__module__ == "gribs.config.loader"
    assert importlib.import_module("gribs.config").validate is entry.validate  # gribs.config is a package now


@pytest.mark.fast
def test_configuration_paths_are_resolved_next_to_the_entry():
    """The two non-verbatim path lines of unit 9: a relative configuration path and the default --config are taken in
    the entry's directory, as before the move (where they read the entry's own __file__)."""
    loader = importlib.import_module("gribs.config.loader"); cli = importlib.import_module("gribs.cli")
    here = _entry_path().resolve().parent
    assert loader._resolve_from_script("results") == here / "results"
    assert cli.build_parser().get_default("config") == here / "gribs_config.json"
    assert cli.build_parser().parse_args([]).config == here / "gribs_config.json"


@pytest.mark.fast
def test_program_directory_of_the_moved_backend_is_the_entry_directory():
    """The relative fcea2 executable path is still looked up next to the entry (report 119, the one non-verbatim line)."""
    thermo = importlib.import_module("gribs.thermo")
    assert thermo._PROGRAM_DIR == _entry_path().resolve().parent
    assert Path(thermo.__file__).resolve() == _entry_path().resolve().parent / "src" / "gribs" / "thermo.py"
