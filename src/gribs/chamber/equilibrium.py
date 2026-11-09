"""Quasi-steady equilibrium pressure (A-6 unit 4, diagnostic only). Moved verbatim from section 6 of the entry
except that burn_rate and nozzle_state, which tests patch and which live in other modules, are read as
gribs.propellant.burning.burn_rate and gribs.nozzle.nozzle_state (ADR decision 1). Tests patch
equilibrium_pressure here: "gribs.chamber.equilibrium.equilibrium_pressure".
"""
from __future__ import annotations

from typing import Optional

import numpy as np
from scipy.optimize import brentq

import gribs.nozzle
import gribs.propellant.burning
from gribs.config import Config
from gribs.errors import GribsError


def equilibrium_pressure(Ab: float, Ri: float, At: float, Ae: float,
                         c: Config) -> Optional[float]:
    """Lowest pressure satisfying mdot_gen_total(p) = mdot_out_total(p).

    Both sides use the TOTAL product mass basis, consistent with the pressure
    ODE; None if no crossing is found.  Diagnostic only.
    """
    if not np.isfinite(Ab) or Ab <= 0.0:
        return None

    def bal(p):
        try:
            v = (c.rho_p * gribs.propellant.burning.burn_rate(p, c, Ab, Ri) * Ab
                 - gribs.nozzle.nozzle_state(p, At, Ae, c)["mdot_total"])
            return float(v) if np.isfinite(v) else np.nan
        except (ValueError, OverflowError, FloatingPointError, GribsError):
            return np.nan

    grid = np.geomspace(c.p_a * (1.0 + 1e-9), 1.0e11, 140)
    vals = np.array([bal(p) for p in grid])
    for p1, p2, f1, f2 in zip(grid[:-1], grid[1:], vals[:-1], vals[1:]):
        if not (np.isfinite(f1) and np.isfinite(f2)):
            continue
        if f1 == 0.0:
            return float(p1)
        if f1 * f2 < 0.0:
            return float(brentq(bal, p1, p2, xtol=1.0, rtol=1e-12))
    return None
