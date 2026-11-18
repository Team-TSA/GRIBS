"""State of the chamber control volume (A-6 stage 2, C4; contract 164 section 2b-4.1, 2b-4.2).

The general form is (E_cv, B_e): the total energy of the contents [J] and the inventory of each element [kg-element]
(the mass of the contents is the sum of B_e).  No reduction to the internal energy U_cv is implemented.
"""
import math
from dataclasses import dataclass
from typing import Dict


@dataclass(frozen=True)
class ChamberState:
    E_cv: float
    B_e: Dict[str, float]

    def __post_init__(self):
        if not math.isfinite(self.E_cv):
            raise ValueError(f"E_cv must be finite, got {self.E_cv!r}")
        if not self.B_e:
            raise ValueError("B_e must hold at least one element")
        for element, kg in self.B_e.items():
            if not (math.isfinite(kg) and kg >= 0.0):
                raise ValueError(f"B_e[{element!r}] must be finite and not negative, got {kg!r}")

    @property
    def mass_kg(self) -> float:
        """Total mass of the contents [kg] (the sum of the element inventories)."""
        return math.fsum(self.B_e.values())
