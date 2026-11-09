"""Test-wide fixtures (A-6, ADR_A6_PACKAGE_SPLIT).

From unit 3b the thermochemistry backend lives in the package module gribs.thermo and is one state shared by every
loaded copy of the entry. An autouse fixture saves it (_thermo_state()) before each test and restores it
(_set_thermo()) after the test, so that a backend installed by one test is never seen by the next one. If gribs.thermo
is first imported during a test, its state is reset afterwards to None, the state of a freshly imported module.
Tests install a backend with _set_thermo() or initialize_thermochemistry() and rely on this guard to undo it; they
never name _THERMO (guard T10). tests/_a6_two_cases_child.py runs two cases in one process WITHOUT this guard (the "two
cases in one process" acceptance item); tests/test_a6_thermo_state_guard.py checks that this guard restores the state.
"""
import sys

import pytest

#: owner of the backend state and of its accessors (ADR decision 2)
THERMO_MODULE = "gribs.thermo"


@pytest.fixture(autouse=True)
def _a6_thermo_state_guard():
    before = sys.modules.get(THERMO_MODULE)
    saved = before._thermo_state() if before is not None else None
    yield
    after = sys.modules.get(THERMO_MODULE)
    if before is not None:
        before._set_thermo(saved)
    if after is not None and after is not before:
        after._set_thermo(None)
