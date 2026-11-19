"""Element mass fractions from explicit reactant formulas (A-6 stage 2, C9; contract 164 sections 2b-FIX-2, 2b-4.9).

The element inventory of the initial fill and the element basis of the reactants (energy/initial.py) need the mass fraction of
each element in each reactant.  The CEA package does not return the element composition of a mixture, and a formula is never
inferred from a CEA species name: the caller states it (energy_model.propellant_formulas: reactant name -> {element symbol:
atoms}).  Element symbols are case-sensitive ("Cl" is chlorine, "CL" is not an element; "Co" is cobalt, "CO" is not an element;
carbon monoxide is {"C": 1, "O": 1}).

Atomic masses: IUPAC standard atomic weights as conventionally tabulated.  They are NOT audited against the atomic masses CEA
uses internally (held item H-17); the difference does not enter a synthetic fixture, but it must be checked before a real case.

Standard library and gribs.errors only.
"""
import math
from typing import Dict, Mapping

from gribs.errors import ConfigurationError

#: Atomic mass [kg/kmol] by element symbol (case-sensitive).  Unaudited against CEA (H-17).
ATOMIC_MASS_KG_PER_KMOL: Dict[str, float] = {
    "H": 1.00794, "Li": 6.941, "B": 10.811, "C": 12.0107, "N": 14.0067, "O": 15.9994, "F": 18.9984, "Na": 22.9898,
    "Mg": 24.305, "Al": 26.9815, "Si": 28.0855, "P": 30.9738, "S": 32.065, "Cl": 35.453, "K": 39.0983, "Ca": 40.078,
    "Ti": 47.867, "Fe": 55.845, "Co": 58.9332, "Cu": 63.546, "Zn": 65.38, "Sr": 87.62, "Zr": 91.224, "Ba": 137.327,
}

_CASE_NOTE = "element symbols are case-sensitive (for example Cl, Co) and formulas are never inferred from species names"


def check_formulas(formulas: Mapping[str, Mapping[str, float]], path: str = "propellant_formulas") -> None:
    """Raise ConfigurationError unless formulas is a non-empty table of non-empty formulas whose element symbols are in
    ATOMIC_MASS_KG_PER_KMOL and whose atom counts are finite positive numbers."""
    if not isinstance(formulas, Mapping) or not formulas:
        raise ConfigurationError(f"{path}: expected a non-empty table reactant name -> {{element: atoms}}")
    for name, formula in formulas.items():
        where = f"{path}[{name!r}]"
        if not isinstance(name, str) or not name:
            raise ConfigurationError(f"{path}: a reactant name must be a non-empty string, got {name!r}")
        if not isinstance(formula, Mapping) or not formula:
            raise ConfigurationError(f"{where}: expected a non-empty {{element: atoms}} formula")
        for symbol, atoms in formula.items():
            if symbol not in ATOMIC_MASS_KG_PER_KMOL:
                raise ConfigurationError(f"{where}: unknown element symbol {symbol!r}; {_CASE_NOTE}")
            if isinstance(atoms, bool) or not isinstance(atoms, (int, float)) or not math.isfinite(atoms) or not atoms > 0.0:
                raise ConfigurationError(f"{where}[{symbol!r}]: the atom count must be a finite number > 0, got {atoms!r}")


def reactant_element_rows(formulas: Mapping[str, Mapping[str, float]]) -> Dict[str, Dict[str, float]]:
    """Per reactant, the mass fraction of each element in that reactant (each row sums to 1): the rows of the element
    basis of energy/initial.py (ElementBasis.from_mapping)."""
    check_formulas(formulas)
    rows: Dict[str, Dict[str, float]] = {}
    for name, formula in formulas.items():
        masses = {e: ATOMIC_MASS_KG_PER_KMOL[e] * float(n) for e, n in formula.items()}
        total = math.fsum(masses.values())
        rows[name] = {e: m / total for e, m in masses.items()}
    return rows


def element_mass_fractions(reactant_wt_percent: Mapping[str, float],
                           formulas: Mapping[str, Mapping[str, float]]) -> Dict[str, float]:
    """The element mass fractions of the mixture (they sum to 1) from the weight percent of each reactant (normalised by
    their sum) and the formulas.  The reactant names of the two tables must be the same set."""
    check_formulas(formulas)
    if not isinstance(reactant_wt_percent, Mapping) or not reactant_wt_percent:
        raise ConfigurationError("the weight percent of the reactants must be a non-empty table")
    missing, extra = sorted(set(reactant_wt_percent) - set(formulas)), sorted(set(formulas) - set(reactant_wt_percent))
    if missing or extra:
        raise ConfigurationError(f"the formulas must be given for exactly the reactants: no formula for {missing}, "
                                 f"formula without a reactant {extra}")
    for name, wt in reactant_wt_percent.items():
        if isinstance(wt, bool) or not isinstance(wt, (int, float)) or not math.isfinite(wt) or wt < 0.0:
            raise ConfigurationError(f"the weight percent of {name!r} must be a finite number >= 0, got {wt!r}")
    total = math.fsum(float(w) for w in reactant_wt_percent.values())
    if not total > 0.0:
        raise ConfigurationError("the weight percent of the reactants has no mass")
    rows = reactant_element_rows(formulas)
    parts: Dict[str, list] = {}
    for name, wt in reactant_wt_percent.items():
        for element, fraction in rows[name].items():
            parts.setdefault(element, []).append(float(wt) / total * fraction)
    return {e: math.fsum(v) for e, v in sorted(parts.items())}
