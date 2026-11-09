"""Test helper (not collected; run by tests/test_a6_thermo_state_guard.py in a child pytest): nothing is imported at
collection; the first test loads the entry (and with it gribs.thermo) and installs a backend without undoing it; the
second expects the state of a freshly imported gribs.thermo (None). With the conftest state guard both pass; with
--noconftest the second fails."""
import importlib.util, sys
from pathlib import Path


def test_1_loads_the_entry_and_installs_a_backend():
    assert "gribs.thermo" not in sys.modules
    bodies = sorted(Path(__file__).resolve().parents[1].glob("GRIBS_v*.py"))
    if len(bodies) != 1:
        raise RuntimeError(f"expected exactly one GRIBS_v*.py, found {[b.name for b in bodies]}")
    spec = importlib.util.spec_from_file_location("gribs_state_guard_late", bodies[0])
    gribs = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = gribs
    spec.loader.exec_module(gribs)
    gribs._set_thermo(object())
    assert gribs._thermo_state() is not None


def test_2_sees_the_state_of_a_fresh_import():
    assert sys.modules["gribs.thermo"]._thermo_state() is None
