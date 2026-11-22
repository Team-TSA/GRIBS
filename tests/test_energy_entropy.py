"""A-6 stage 2, U-H25-5 stage A (reports 272, 273; PRE_REGISTRATION_R545): the specific entropy in the TP state of the provider.
Structure and finiteness only: no tolerance is asserted for a CEA value (decision 110, T3: the basis of the entropy is recorded by
the run's tool, not decided here).  Needs the official NASA CEA package (the whole module is skipped without it (A-1 is in test_energy_isentrope.py, which needs no CEA), recorded as "not
executable (reason: cea not installed)").  No monkeypatch, no assignment.
"""
import json
import math

import pytest

cea = pytest.importorskip("cea")

from gribs.energy.cea_tp import TP_CACHE_PREFIX, TP_CACHE_SCHEMA, CEATPProvider, TPCache  # noqa: E402
from test_energy_cea_tp import B_E, ELEMENTS, NAMES, WT  # noqa: E402


def provider(cache=None):
    return CEATPProvider(NAMES, WT, ELEMENTS, cache=cache)


@pytest.mark.fast
def test_a02_the_cea_entropy_is_finite_and_survives_the_cache_round_trip_bit_for_bit(tmp_path):
    a = provider(cache=TPCache(tmp_path)).tp(3300.0, 2.0e6, B_E)
    second = TPCache(tmp_path)
    b = provider(cache=second).tp(3300.0, 2.0e6, B_E)
    assert second.hits == 1 and math.isfinite(a.entropy_J_per_kgK) and a.entropy_J_per_kgK == b.entropy_J_per_kgK
    assert provider().tp_detail(3300.0, 2.0e6, B_E).entropy_J_per_kgK == a.entropy_J_per_kgK


@pytest.mark.fast
def test_a04_one_tp_call_is_one_solve():
    prov = provider()
    real = prov._cea
    made = []

    class CountingModule:
        def __getattr__(self, name):
            return getattr(real, name)

        def EqSolution(self, solver):                              # one new EqSolution per solve (C18)
            made.append(1)
            return real.EqSolution(solver)
    prov._cea = CountingModule()
    prov.tp(3300.0, 2.0e6, B_E)
    assert len(made) == 1


@pytest.mark.fast
def test_a05_a_file_of_the_schema_3_is_not_read_and_the_state_is_solved_again(tmp_path):
    assert TP_CACHE_SCHEMA == "gribs-energy-tp-cache-4"
    cache = TPCache(tmp_path)
    prov = provider(cache=cache)
    prov.tp(3300.0, 2.0e6, B_E)
    (path,) = list(tmp_path.glob(TP_CACHE_PREFIX + "*.json"))
    doc = json.loads(path.read_text(encoding="utf-8"))
    path.unlink()
    doc["schema"] = "gribs-energy-tp-cache-3"
    del doc["state"]["entropy_J_per_kgK"]
    old = tmp_path / (TP_CACHE_PREFIX + "0123456789abcdef01234567.json")
    old.write_text(json.dumps(doc), encoding="utf-8")
    before = old.read_bytes()
    fresh = TPCache(tmp_path)
    provider(cache=fresh).tp(3300.0, 2.0e6, B_E)
    assert fresh.hits == 0 and fresh.stores == 1 and old.read_bytes() == before
