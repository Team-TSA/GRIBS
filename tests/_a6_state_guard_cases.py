"""Test helper (not collected; run by tests/test_a6_thermo_state_guard.py in a child pytest): the first test installs
a backend and does not undo it; the second expects the state from before the first. With the conftest state guard both
pass; with --noconftest the second fails."""
import importlib.util, sys
from pathlib import Path

_BODIES = sorted(Path(__file__).resolve().parents[1].glob("GRIBS_v*.py"))
if len(_BODIES) != 1:
    raise RuntimeError(f"expected exactly one GRIBS_v*.py, found {[b.name for b in _BODIES]}")
spec = importlib.util.spec_from_file_location("gribs_state_guard_cases", _BODIES[0])
gribs = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = gribs
spec.loader.exec_module(gribs)
BEFORE = gribs._thermo_state()


def test_1_installs_a_backend_and_leaves_it():
    gribs._set_thermo(object())
    assert gribs._thermo_state() is not BEFORE


def test_2_sees_the_state_from_before_test_1():
    assert gribs._thermo_state() is BEFORE
