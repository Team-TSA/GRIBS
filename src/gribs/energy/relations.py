"""Pure relations between the mass-basis quantities of the chamber model (contract 164, section 1; F-21).

  R_eff = Ru / M_CEA = Yg * R_gas        [J/(kg-total K)]  (the value gas_props returns)
  R_gas = Ru / M_gas = R_eff / Yg        [J/(kg-gas K)]
  M_gas = Yg * M_CEA                     [kg-gas/kmol-gas]  (M_CEA = m_total / n_gas, kg-total/kmol-gas)
  Psi   = R_eff * T_used                 [J/kg]

Standard library only; no import of gribs modules.
"""


def r_gas(R_eff: float, Yg: float) -> float:
    """Gas-phase specific gas constant R_gas = R_eff / Yg [J/(kg-gas K)]."""
    if not Yg > 0.0:
        raise ValueError(f"Yg must be positive, got {Yg!r}")
    return R_eff / Yg


def m_gas(M_cea: float, Yg: float) -> float:
    """Mean molecular weight of the gas phase M_gas = Yg * M_CEA [kg-gas/kmol-gas]."""
    return Yg * M_cea


def psi(R_eff: float, T_used: float) -> float:
    """Psi = R_eff * T_used [J/kg] (= Pa m^3/kg)."""
    return R_eff * T_used
