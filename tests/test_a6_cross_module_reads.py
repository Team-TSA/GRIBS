"""A-6 unit 4 (ADR decision 1): a package module reads a patched name of another package module through the owner
module attribute, so that one patch on the owner reaches every reader. equilibrium_pressure (gribs.chamber.equilibrium)
reads burn_rate (gribs.propellant.burning) and nozzle_state (gribs.nozzle)."""
import dataclasses
import importlib
import math

import pytest


def _config():
    gribs_config = importlib.import_module("gribs.config")
    return dataclasses.replace(gribs_config.Config(), rho_p=1000.0)  # the default runtime density is 0 (set from the case)


@pytest.mark.fast
def test_patch_on_the_owner_reaches_equilibrium_pressure(monkeypatch):
    eq = importlib.import_module("gribs.chamber.equilibrium")
    c = _config()
    seen = {"burn_rate": 0, "nozzle_state": 0}

    def fake_burn_rate(p, config, Ab, Ri):
        seen["burn_rate"] += 1
        return 1.0 / (config.rho_p * Ab)          # generation = 1 kg/s at every pressure

    def fake_nozzle_state(p, At, Ae, config):
        seen["nozzle_state"] += 1
        return {"mdot_total": p / 4.0e6}          # outflow = 1 kg/s at 4 MPa

    monkeypatch.setattr("gribs.propellant.burning.burn_rate", fake_burn_rate)
    monkeypatch.setattr("gribs.nozzle.nozzle_state", fake_nozzle_state)
    p = eq.equilibrium_pressure(1.0e-2, 1.0e-2, 1.0e-4, 4.0e-4, c)
    assert seen["burn_rate"] > 0 and seen["nozzle_state"] > 0
    assert p is not None and math.isclose(p, 4.0e6, rel_tol=1e-6)
