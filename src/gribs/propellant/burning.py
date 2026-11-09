"""Burn rate (A-6 unit 4): Saint-Robert law with temperature sensitivity and the optional erosive contribution.
Moved verbatim from section 4 of the entry. Tests patch burn_rate here: "gribs.propellant.burning.burn_rate";
the module is not named burn_rate so that the module and the function never share a name (guard T7).
"""
from __future__ import annotations

import math

from gribs.config import Config


def base_burn_rate(p0: float, c: Config) -> float:
    return (c.a_burn * (max(p0, 0.0) / c.p_ref) ** c.n_burn
            * math.exp(c.sigma_p * (c.T_grain - c.T_ref)))


def burn_rate(p0: float, c: Config, Ab: float, Ri: float) -> float:
    """Total regression rate including the optional erosive contribution."""
    r0 = base_burn_rate(p0, c)
    if c.ero_alpha <= 0.0 or Ri <= 0.0 or Ab <= 0.0:
        return r0
    A_port = math.pi * Ri ** 2
    D_h = 2.0 * Ri
    r = r0
    for _ in range(60):                     # damped fixed-point iteration
        G = c.rho_p * r * Ab / A_port       # port mass flux at the aft end
        if G <= 1.0e-12:
            return r0
        r_new = r0 + c.ero_alpha * G ** 0.8 / D_h ** 0.2 \
            * math.exp(-c.ero_beta * c.rho_p * r / G)
        if abs(r_new - r) <= 1.0e-13 * max(r, 1.0e-12):
            return r_new
        r = 0.5 * (r + r_new)
    return r
