"""A-6 stage 2, C6: the initial fill and its checks at t = 0 (contract 164 sections 2b-FIX-2, 2b-4.7, 2b-4.9;
PRE_REGISTRATION_C6).  A synthetic basis and a synthetic provider (assumed inputs, not data; held items H-04, H-08); the
tolerance is an argument.  No assignment, no monkeypatch."""
import dataclasses

import pytest

from gribs.constants import R_UNIVERSAL
from gribs.energy import initial
from gribs.energy.config_model import CondensedData, CondensedPhaseSpec
from gribs.energy.exceptions import InitialFillError
from gribs.energy.initial import ElementBasis, check_initial_fill, element_weights, initial_state_from_fill
from gribs.energy.provider import TPState
from gribs.energy.state import ChamberState

TOL = 1.0e-6
R1 = {"C": 0.5, "H": 0.1, "O": 0.4}
R2 = {"C": 0.2, "H": 0.05, "O": 0.15, "Al": 0.6}
R3 = {"C": 0.1, "H": 0.1, "O": 0.8}
BASIS = ElementBasis.from_mapping({"R1": R1, "R2": R2})
# the main-propellant products: 0.7 R1 + 0.3 R2 (element mass fractions)
PRODUCTS = {e: 0.7 * R1.get(e, 0.0) + 0.3 * R2.get(e, 0.0) for e in ("C", "H", "O", "Al")}


def mix(rows):
    return {e: sum(w * r.get(e, 0.0) for w, r in rows) for e in ("C", "H", "O", "Al")}


@pytest.mark.fast
def test_i01_an_air_fill_is_rejected():
    air = {"N": 0.77, "O": 0.23}
    with pytest.raises(InitialFillError, match=r"\(i\) .*N"):
        check_initial_fill(air, PRODUCTS, BASIS, element_tolerance=TOL)


@pytest.mark.fast
def test_i02_an_initial_composition_that_differs_from_the_products_is_rejected():
    other = mix([(0.6, R1), (0.4, R2)])                                # covered and reproducible, but not the products
    assert element_weights(BASIS, other, element_tolerance=TOL)["R1"] == pytest.approx(0.6)
    with pytest.raises(InitialFillError, match=r"\(iii\)"):
        check_initial_fill(other, PRODUCTS, BASIS, element_tolerance=TOL)


@pytest.mark.fast
def test_i03_a_fill_of_the_products_is_accepted():
    weights = check_initial_fill({e: 2.5 * f for e, f in PRODUCTS.items()}, PRODUCTS, BASIS, element_tolerance=TOL)   # scale is free
    assert weights["R1"] == pytest.approx(0.7, abs=1e-9) and weights["R2"] == pytest.approx(0.3, abs=1e-9)


@pytest.mark.fast
def test_i04_a_basis_that_does_not_cover_a_positive_element_is_rejected():
    only_r1 = ElementBasis.from_mapping({"R1": R1})
    with pytest.raises(InitialFillError, match=r"\(i\) .*Al"):
        element_weights(only_r1, PRODUCTS, element_tolerance=TOL)


@pytest.mark.fast
def test_i05_weights_that_cannot_be_non_negative_are_rejected():
    basis = ElementBasis.from_mapping({"R1": {**R1, "Al": 0.0}, "R3": {**R3, "Al": 0.0}})
    target = mix([(1.5, R1), (-0.5, R3)])                              # needs w = (1.5, -0.5); every fraction is positive
    assert all(v >= 0.0 for v in target.values())
    with pytest.raises(InitialFillError, match="negative weight"):
        element_weights(basis, target, element_tolerance=TOL)


@pytest.mark.fast
def test_i06_a_basis_with_non_unique_weights_is_rejected():
    dependent = {"R1": R1, "R2": R2, "R4": {e: 0.5 * R1.get(e, 0.0) + 0.5 * R2.get(e, 0.0) for e in ("C", "H", "O", "Al")}}
    with pytest.raises(InitialFillError, match="unique"):
        element_weights(ElementBasis.from_mapping(dependent), PRODUCTS, element_tolerance=TOL)


@dataclasses.dataclass(frozen=True)
class Products:
    """A synthetic equilibrium: ideal gas plus one condensed species (assumed inputs)."""
    converged: bool = True

    def tp(self, T, p, B_e):
        return TPState(u_J_per_kg=-3.0e6 + 900.0 * (T - 300.0), n_gas_per_kg=0.6 / 28.0, Yg=0.6,
                       condensed_mass_fraction={"OX(s)": 0.4}, converged=self.converged)


CONDENSED = CondensedData("neglected", (CondensedPhaseSpec("OX(s)", "solid", 3990.0, 1300.0, -1.6e7, 298.15),), ())


@pytest.mark.fast
def test_i07_the_state_follows_from_the_fill_temperature_pressure_and_composition():
    V, T, p = 0.02, 300.0, 2.0e5
    state = initial_state_from_fill(T, p, PRODUCTS, V, Products(), CONDENSED)
    M = state.mass_kg
    assert state.B_e["Al"] / M == pytest.approx(PRODUCTS["Al"], rel=1e-12)
    v_gas = M * (0.6 / 28.0) * R_UNIVERSAL * T / p
    assert v_gas == pytest.approx(V, rel=1e-12)                        # V0: the gas fills V_free (the condensed volume is 0)
    assert M == pytest.approx(V * p / ((0.6 / 28.0) * R_UNIVERSAL * T), rel=1e-12)      # M = V p/(n Ru T)
    assert state.E_cv == pytest.approx(M * (-3.0e6), rel=1e-12)        # E_cv = M u with u at T_fill_0 = 300 K
    assert initial._RU == R_UNIVERSAL
    with pytest.raises(InitialFillError, match="did not converge"):
        initial_state_from_fill(T, p, PRODUCTS, V, Products(converged=False), CONDENSED)
    with pytest.raises(ValueError):
        ChamberState(float("nan"), {"C": 1.0})
    with pytest.raises(ValueError):
        ChamberState(0.0, {"C": -1.0})
    with pytest.raises(ValueError):
        ChamberState(0.0, {})


@pytest.mark.fast
def test_i08_the_element_tolerance_separates_accepted_from_rejected():
    shift = lambda d: {**PRODUCTS, "C": PRODUCTS["C"] + d, "O": PRODUCTS["O"] - d}      # the sum of the fractions stays 1
    check_initial_fill(shift(0.5 * TOL), PRODUCTS, BASIS, element_tolerance=TOL)         # inside the tolerance: accepted
    with pytest.raises(InitialFillError):
        check_initial_fill(shift(5.0 * TOL), PRODUCTS, BASIS, element_tolerance=TOL)     # outside: rejected
