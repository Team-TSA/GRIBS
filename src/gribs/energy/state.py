"""State of the chamber control volume (A-6 stage 2, C4; contract 164 section 2b-4.1, 2b-4.2).

The general form is (E_cv, B_e): the total energy of the contents [J] and the inventory of each element [kg-element]
(the mass of the contents is the sum of B_e).  No reduction to the internal energy U_cv is implemented.
RecoveredState (C5) is what the state recovery returns for such a state.
"""
import math
from dataclasses import dataclass, field
from typing import Dict, Optional


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


@dataclass(frozen=True)
class RecoveredState:
    """The state recovered from (E_cv, V_free, B_e): pressure and temperature, and on a phase-transition plateau the phase
    fraction phi (0 = all solid, 1 = all liquid; None off a plateau).  V_gas + V_condensed = V_free."""
    p_Pa: float
    T_K: float
    phi: Optional[float]
    on_plateau: bool
    m_total_kg: float
    V_gas_m3: float
    V_condensed_m3: float
    n_gas_kmol: float
    Yg: float
    condensed_mass_fraction: Dict[str, float] = field(default_factory=dict)
    n_starts_agreed: int = 1
