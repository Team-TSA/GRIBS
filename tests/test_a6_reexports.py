"""A-6: every top-level name of the moved low-level modules is re-exported by the entry, and the entry and the
package hold the same objects (ADR_A6_PACKAGE_SPLIT, unit 1)."""
import ast
import importlib
import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
MODULES = ("gribs.constants", "gribs.errors", "gribs.util", "gribs.config")


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
    missing = {m: sorted(_top_level_names(m) - imports.get(m, set())) for m in MODULES}
    assert not any(missing.values()), missing


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
