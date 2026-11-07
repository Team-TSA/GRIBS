"""Test-wide fixtures (A-6, ADR_A6_PACKAGE_SPLIT).

From unit 4a an autouse fixture saves and restores the thermochemistry backend state (_thermo_state())
around each test, because every loaded copy of the entry then shares one package state. Until then the
fixture does nothing. A subprocess test that does not restore state is kept for the "two cases in one
process" acceptance item.
"""
import pytest


@pytest.fixture(autouse=True)
def _a6_thermo_state_guard():
    # unit 0: no-op; the save / restore is added in unit 4a.
    yield
