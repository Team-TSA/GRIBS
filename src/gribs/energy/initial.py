"""The initial fill of energy mode and its checks at t = 0 (A-6 stage 2, C6; contract 164 sections 2b-FIX-2, 2b-4.9).

energy mode starts from a charge that is the equilibrium products of the main propellant (an air fill is out of scope), at
the temperature T_fill_0 and the fill pressure.  At t = 0 it checks (i) that the reactant basis covers every element with a
positive inventory, (ii) that non-negative weights reproduce the element inventory (a basis whose weights are not unique is
rejected, contract 2b-4.7), and (iii) that the element fractions equal those of the main-propellant products within the
configured element_tolerance.  The state then follows from T_fill_0, the fill pressure and the composition.

The tolerance is an argument (the configured energy_model.initial_fill.element_tolerance, held item H-04): nothing here
defaults it.  numpy is used (energy path only).
"""
import math
from dataclasses import dataclass
from typing import Dict, Mapping, Tuple

import numpy as np

from gribs.energy.config_model import CondensedData
from gribs.energy.exceptions import InitialFillError
from gribs.energy.provider import EquilibriumProvider
from gribs.energy.state import ChamberState

#: Universal gas constant [J/(kmol K)] (the value of gribs.constants.R_UNIVERSAL; a test compares them).
_RU = 8314.51


@dataclass(frozen=True)
class ElementBasis:
    """The reactants of the basis and the mass fraction of each element in each reactant (rows: reactants)."""
    reactants: Tuple[str, ...]
    elements: Tuple[str, ...]
    a_ej: Tuple[Tuple[float, ...], ...]

    @classmethod
    def from_mapping(cls, composition: Mapping[str, Mapping[str, float]]) -> "ElementBasis":
        reactants = tuple(composition)
        elements = tuple(sorted({e for row in composition.values() for e in row}))
        return cls(reactants, elements, tuple(tuple(float(composition[r].get(e, 0.0)) for e in elements) for r in reactants))


def _fractions(inventory: Mapping[str, float], what: str) -> Dict[str, float]:
    for element, value in inventory.items():
        if not (math.isfinite(value) and value >= 0.0):
            raise InitialFillError(f"{what}[{element!r}] must be finite and not negative, got {value!r}")
    total = math.fsum(inventory.values())
    if not total > 0.0:
        raise InitialFillError(f"{what} has no mass")
    return {e: v / total for e, v in inventory.items()}


def element_weights(basis: ElementBasis, element_fractions: Mapping[str, float], *, element_tolerance: float) -> Dict[str, float]:
    """The weights of the reactants that reproduce the element fractions: (i) every element with a positive fraction must be
    in the basis; the basis must give unique weights (contract 2b-4.7); (ii) the weights must be non-negative (none is
    clipped) and reproduce the fractions within element_tolerance."""
    uncovered = sorted(e for e, f in element_fractions.items() if f > 0.0 and e not in basis.elements)
    if uncovered:
        raise InitialFillError(f"(i) the reactant basis does not cover the element(s) {uncovered} that the initial fill contains")
    A = np.array(basis.a_ej, dtype=float).T                           # elements x reactants
    b = np.array([element_fractions.get(e, 0.0) for e in basis.elements], dtype=float)
    if np.linalg.matrix_rank(A) < A.shape[1]:
        raise InitialFillError("(ii) the reactant basis does not give unique weights (its reactants are linearly dependent); "
                               "a non-unique weight is rejected (contract 2b-4.7)")
    w, *_ = np.linalg.lstsq(A, b, rcond=None)
    if float(np.max(np.abs(A @ w - b))) > element_tolerance:
        raise InitialFillError("(ii) no weights reproduce the element fractions of the initial fill within the tolerance")
    if float(np.min(w)) < -element_tolerance:
        raise InitialFillError(f"(ii) reproducing the element fractions needs a negative weight ({float(np.min(w)):.6g}); "
                               f"non-negative weights do not reproduce the initial fill")
    return dict(zip(basis.reactants, (float(x) for x in w)))


def check_initial_fill(B_e0: Mapping[str, float], product_element_fractions: Mapping[str, float], basis: ElementBasis,
                       *, element_tolerance: float) -> Dict[str, float]:
    """The t = 0 checks (i), (ii) and (iii); returns the reactant weights.  product_element_fractions are the element mass
    fractions of the main-propellant products."""
    fill = _fractions(B_e0, "B_e(0)")
    weights = element_weights(basis, fill, element_tolerance=element_tolerance)
    products = _fractions(product_element_fractions, "the product element fractions")
    differ = {e: (fill.get(e, 0.0), products.get(e, 0.0)) for e in set(fill) | set(products)
              if abs(fill.get(e, 0.0) - products.get(e, 0.0)) > element_tolerance}
    if differ:
        raise InitialFillError(
            f"(iii) the element fractions of the initial fill differ from those of the main-propellant products by more than "
            f"{element_tolerance:g}: {dict(sorted(differ.items()))}; energy mode accepts only a fill of the products "
            f"(an air fill is out of scope)")
    return weights


def initial_state_from_fill(T_fill_0: float, p_fill: float, element_fractions: Mapping[str, float], V_free: float,
                            provider: EquilibriumProvider, condensed: CondensedData) -> ChamberState:
    """The state at t = 0 from T_fill_0, the fill pressure and the composition (the supply form (a)): the TP equilibrium of
    the composition gives the gas volume per kg of total (V0: the condensed volume is taken as 0), so M = V_free / (n Ru T / p),
    B_e(0) = M times the element fractions and E_cv = M u.  A p0-only start is not possible (the temperature and the composition are required)."""
    for name, value in (("T_fill_0", T_fill_0), ("p_fill", p_fill), ("V_free", V_free)):
        if not (math.isfinite(value) and value > 0.0):
            raise InitialFillError(f"{name} must be positive and finite, got {value!r}")
    fractions = _fractions(element_fractions, "the element fractions")
    st = provider.tp(T_fill_0, p_fill, fractions)
    if not st.converged:
        raise InitialFillError(f"the equilibrium of the initial fill did not converge at T = {T_fill_0:g} K, p = {p_fill:g} Pa")
    rho = {p.species: p.density_kg_m3 for p in condensed.phases}
    for species in st.condensed_mass_fraction:
        if species not in rho:
            raise InitialFillError(f"no density for the condensed species {species!r} in energy_model.condensed")
    M = V_free / (st.n_gas_per_kg * _RU * T_fill_0 / p_fill)
    return ChamberState(E_cv=M * st.u_J_per_kg, B_e={e: M * f for e, f in fractions.items()})
