"""The equilibrium provider the state recovery calls (A-6 stage 2, C5; contract 164 section 3-A).

A provider answers a TP equilibrium at (T, p) for the element inventory B_e.  Quantities are per kg of TOTAL products
(gas plus condensed): u_J_per_kg is the internal energy of the mixture, n_gas_per_kg = 1 / M_CEA [kmol-gas / kg-total]
(so that p V_gas = m_total R_eff T with R_eff = Ru n_gas_per_kg), condensed_mass_fraction maps each condensed species to
its mass fraction, gamma_s is the equilibrium isentropic exponent (> 1 when given).  The condensed-phase volume is taken as 0, as in CEA (RP-1311 p.4): V, n and the densities are
those of the gas, the mass is that of the whole mixture; the constant densities of energy_model.condensed give a diagnostic only.

The analytic provider (reference_constant) and the synthetic providers of the tests are for tests: a result obtained with
them is not a result for a propellant.
"""
from dataclasses import dataclass, field
from typing import Dict, Mapping, Optional, Protocol


@dataclass(frozen=True)
class TPState:
    u_J_per_kg: float
    n_gas_per_kg: float
    Yg: float
    condensed_mass_fraction: Dict[str, float] = field(default_factory=dict)
    converged: bool = True
    #: Equilibrium isentropic exponent of the mixture (CEA EqSolution.gamma_s, the definition of the legacy gamma_s); None when
    #: the provider does not give it (a run that needs it - the nozzle flow - then fails fast, C11).
    gamma_s: Optional[float] = None
    #: Specific entropy of the mixture [J/(kg K)] per kg of TOTAL products, the value of CEA (EqSolution.entropy x 1000; the pressure
    #: term is included).  None when the provider does not give it; the numerical isentrope (gribs.energy.isentrope, A-6 stage 2,
    #: U-H25-5) then fails fast.  The basis (per kg of total products, like energy and enthalpy) is the assumption of this field.
    entropy_J_per_kgK: Optional[float] = None


class EquilibriumProvider(Protocol):
    def tp(self, T: float, p: float, B_e: Mapping[str, float]) -> TPState:
        """TP equilibrium at temperature T [K] and pressure p [Pa] for the element inventory B_e [kg-element]."""
        ...
