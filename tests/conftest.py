"""Test-wide fixtures (A-6, ADR_A6_PACKAGE_SPLIT).

From unit 3a an autouse fixture saves the thermochemistry backend (_thermo_state()) of every loaded copy of the entry
before each test and restores it (_set_thermo()) after the test, so that a backend installed by one test is never seen
by the next one. Tests install a backend with _set_thermo() or initialize_thermochemistry() and rely on this guard to
undo it; they never assign _THERMO (guard T10). From unit 3b the state lives in the package and is shared by every copy
of the entry; the guard then saves and restores the package state. tests/_a6_two_cases_child.py runs two cases in one
process WITHOUT this guard (the "two cases in one process" acceptance item).
"""
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _entry_path():
    bodies = sorted(ROOT.glob("GRIBS_v*.py"))
    if len(bodies) != 1:
        raise RuntimeError(f"expected exactly one GRIBS_v*.py, found {[b.name for b in bodies]}")
    return bodies[0].resolve()


def _loaded_entries():
    """Every module in sys.modules whose file is the entry and that has finished defining the accessors."""
    entry = _entry_path()
    found = []
    for module in list(sys.modules.values()):
        path = getattr(module, "__file__", None)
        if path and Path(path).resolve() == entry and callable(getattr(module, "_set_thermo", None)):
            found.append(module)
    return found


@pytest.fixture(autouse=True)
def _a6_thermo_state_guard():
    saved = [(module, module._thermo_state()) for module in _loaded_entries()]
    yield
    for module, backend in saved:
        module._set_thermo(backend)
