"""A-6: every top-level name of the moved low-level modules is re-exported by the entry, and the entry and the
package hold the same objects (ADR_A6_PACKAGE_SPLIT, unit 1)."""
import ast
import importlib
import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
MODULES = ("gribs.constants", "gribs.errors", "gribs.util", "gribs.config", "gribs.thermo",
           "gribs.propellant.burning", "gribs.nozzle", "gribs.chamber.equilibrium")


def _entry_path():
    bodies = sorted(ROOT.glob("GRIBS_v*.py"))
    if len(bodies) != 1:
        raise RuntimeError(f"expected exactly one GRIBS_v*.py, found {[b.name for b in bodies]}")
    return bodies[0]


def _top_level_names(module):
    tree = ast.parse((ROOT / "src" / Path(*module.split("."))).with_suffix(".py").read_text())
    out = set()
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.ClassDef)):
            out.add(node.name)
        elif isinstance(node, (ast.Assign, ast.AnnAssign)):
            for t in (node.targets if isinstance(node, ast.Assign) else [node.target]):
                if isinstance(t, ast.Name):
                    out.add(t.id)
    return out


def _state_names(module):
    """Module state: names the module rebinds with `global`. They are never re-exported (a from-import would copy the
    value once and go stale) and are reached through their accessors only (ADR decision 2, guard T10)."""
    tree = ast.parse((ROOT / "src" / Path(*module.split("."))).with_suffix(".py").read_text())
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
def test_program_directory_of_the_moved_backend_is_the_entry_directory():
    """The relative fcea2 executable path is still looked up next to the entry (report 119, the one non-verbatim line)."""
    thermo = importlib.import_module("gribs.thermo")
    assert thermo._PROGRAM_DIR == _entry_path().resolve().parent
    assert Path(thermo.__file__).resolve() == _entry_path().resolve().parent / "src" / "gribs" / "thermo.py"
