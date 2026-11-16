"""P5-4 (F-02): the erosive-burning burn rate is solved to convergence, or the call raises.

With erosive burning, r = r0 + K*r**0.8 (K does not depend on r, because the port mass flux is proportional to r).
For an erosive fraction phi = (r - r0)/r the exact root is r = r0/(1 - phi), reached with K = phi*r**0.2; the tests
choose alpha so that K takes that value. Before P5-4 a damped fixed-point iteration stopped after 60 steps and
returned its last value, which was off by up to about 1 % near phi = 1.  The configuration is golden g1
(reference_constant, no CEA) loaded through the JSON loader.
"""
import copy
import json
import math
from pathlib import Path

import pytest

import gribs.propellant.burning as burning
from gribs.config.loader import configuration_from_document
from gribs.errors import ErosiveRootNotConvergedError, GribsError

ROOT = Path(__file__).resolve().parents[1]
G1 = json.loads((ROOT / "cases" / "golden" / "g1_reference_constant.json").read_text(encoding="utf-8"))
RI, LP, RP = 4.0e-3, 0.10, 0.015
AB = 2.0 * math.pi * RI * LP + math.pi * (RP ** 2 - RI ** 2)


def _config(erosive):
    doc = copy.deepcopy(G1)
    doc["propellant"]["burn_law"]["erosive_burning"] = erosive
    return configuration_from_document(doc).config


def _alpha_for_phi(c, p0, phi, beta):
    """alpha that gives the erosive fraction phi at (p0, AB, RI); the exact root is r0/(1 - phi)."""
    r0 = burning.base_burn_rate(p0, c)
    r = r0 / (1.0 - phi)
    K = phi * r ** 0.2
    A_port = math.pi * RI ** 2
    unit = (c.rho_p * AB / A_port) ** 0.8 / (2.0 * RI) ** 0.2 * math.exp(-beta * A_port / AB)
    return K / unit, r


@pytest.mark.fast
@pytest.mark.parametrize("phi", [1e-6, 0.1, 0.35, 0.48, 0.75, 0.94, 0.999])
@pytest.mark.parametrize("beta", [53.0, 0.0])
def test_erosive_root_matches_closed_form(phi, beta):
    p0 = 5.0e6
    alpha, r_exact = _alpha_for_phi(_config({"enabled": False, "alpha": 0.0, "beta": beta}), p0, phi, beta)
    c = _config({"enabled": True, "alpha": alpha, "beta": beta})
    r = burning.burn_rate(p0, c, AB, RI)
    # alpha is rounded to a double, so the exact root of the rounded problem differs from r_exact by a few ulps
    # times the condition number phi/(1 - 0.8*phi) <= 5; 1e-12 is the B2 criterion (LEDGER section 3).
    assert abs(r - r_exact) <= 1.0e-12 * r_exact


@pytest.mark.fast
def test_erosive_root_satisfies_the_unchanged_equation():
    """The returned r solves r = r0 + alpha*G**0.8/(2Ri)**0.2*exp(-beta*rho_p*r/G), G = rho_p*r*Ab/(pi*Ri**2)."""
    c = _config({"enabled": True, "alpha": 3.2e-5, "beta": 53.0})
    for p0 in (1.0e6, 1.5e6, 7.0e6):
        r0 = burning.base_burn_rate(p0, c)
        r = burning.burn_rate(p0, c, AB, RI)
        G = c.rho_p * r * AB / (math.pi * RI ** 2)
        rhs = r0 + c.ero_alpha * G ** 0.8 / (2.0 * RI) ** 0.2 * math.exp(-c.ero_beta * c.rho_p * r / G)
        assert r > r0
        assert abs(r - rhs) <= 1.0e-13 * r


@pytest.mark.fast
def test_disabled_erosion_returns_the_base_rate_bit_for_bit():
    on = _config({"enabled": True, "alpha": 1.0e-4, "beta": 53.0})
    off = _config({"enabled": False, "alpha": 1.0e-4, "beta": 53.0})
    zero = _config({"enabled": True, "alpha": 0.0, "beta": 53.0})
    for p0 in (1.0e5, 1.0e6, 2.0e7):
        r0 = burning.base_burn_rate(p0, on)
        for c, ab, ri in ((off, AB, RI), (zero, AB, RI), (on, 0.0, RI), (on, AB, 0.0)):
            assert burning.burn_rate(p0, c, ab, ri).hex() == r0.hex()


@pytest.mark.fast
def test_non_convergence_raises_and_returns_no_value():
    c = _config({"enabled": True, "alpha": 3.2e-5, "beta": 53.0})
    r0 = burning.base_burn_rate(1.5e6, c)
    A_port = math.pi * RI ** 2
    with pytest.raises(ErosiveRootNotConvergedError) as info:
        burning._erosive_root(r0, c, AB, A_port, 2.0 * RI, 0)
    err = info.value
    assert isinstance(err, GribsError)                      # the command line reports GribsError from the solver as exit 4
    assert err.code == "E_EROSIVE_ROOT_NOT_CONVERGED" and str(err).startswith(err.code)
    assert err.context["max_iter"] == 0 and err.context["r0_m_s"] == r0
    assert burning._erosive_root(r0, c, AB, A_port, 2.0 * RI, burning.EROSIVE_ROOT_MAX_ITER) > r0
