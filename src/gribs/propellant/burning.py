"""Burn rate (A-6 unit 4): Saint-Robert law with temperature sensitivity and the optional erosive contribution.
Moved verbatim from section 4 of the entry. Tests patch burn_rate here: "gribs.propellant.burning.burn_rate";
the module is not named burn_rate so that the module and the function never share a name (guard T7).
"""
from __future__ import annotations

import math

from gribs.config import Config
from gribs.errors import ErosiveRootNotConvergedError

#: Newton iteration limit of the erosive-burning root (P5-4, F-02). A fixed code value, not a configuration item:
#: from the upper bound of the root the iteration needs about 10 steps (at most a few tens near an erosive
#: fraction of 1). Read at each call.
EROSIVE_ROOT_MAX_ITER = 100
#: Convergence test of the root: |r_new - r| <= EROSIVE_ROOT_RTOL * r_new, followed by one more Newton step.
EROSIVE_ROOT_RTOL = 1.0e-13


def base_burn_rate(p0: float, c: Config) -> float:
    return (c.a_burn * (max(p0, 0.0) / c.p_ref) ** c.n_burn
            * math.exp(c.sigma_p * (c.T_grain - c.T_ref)))


def burn_rate(p0: float, c: Config, Ab: float, Ri: float) -> float:
    """Total regression rate including the optional erosive contribution.

    With erosive burning, r solves r = r0 + e(r), where
    e(r) = alpha * G**0.8 / D_h**0.2 * exp(-beta * rho_p * r / G) and G = rho_p * r * Ab / A_port is the port
    mass flux at the aft end. Because G is proportional to r, e(r) = K * r**0.8 with K independent of r, so
    g(r) = r - r0 - e(r) is convex with one positive root. Newton's method is started at the upper bound
    (r0**0.2 + K)**5 of the root and decreases monotonically to it. If it does not converge, the function raises
    ErosiveRootNotConvergedError and returns no value (P5-4, F-02).
    """
    r0 = base_burn_rate(p0, c)
    if c.ero_alpha <= 0.0 or Ri <= 0.0 or Ab <= 0.0:
        return r0
    A_port = math.pi * Ri ** 2
    D_h = 2.0 * Ri
    if c.rho_p * r0 * Ab / A_port <= 1.0e-12:   # port mass flux at r = r0
        return r0
    return _erosive_root(r0, c, Ab, A_port, D_h, EROSIVE_ROOT_MAX_ITER)


def _erosive_root(r0: float, c: Config, Ab: float, A_port: float, D_h: float, max_iter: int) -> float:
    def erosive(r):
        G = c.rho_p * r * Ab / A_port       # port mass flux at the aft end
        return c.ero_alpha * G ** 0.8 / D_h ** 0.2 * math.exp(-c.ero_beta * c.rho_p * r / G)

    r = math.nan
    try:
        K = erosive(r0) / r0 ** 0.8
        r = (r0 ** 0.2 + K) ** 5            # upper bound of the root
        done = False
        for _ in range(max_iter):
            e = erosive(r)
            g = r - r0 - e
            dg = 1.0 - 0.8 * e / r
            if not (math.isfinite(g) and dg > 0.0):
                break
            r_new = r - g / dg
            if not (r_new > 0.0 and math.isfinite(r_new)):
                break
            if done:
                return r_new
            done = abs(r_new - r) <= EROSIVE_ROOT_RTOL * r_new
            r = r_new
    except (OverflowError, ZeroDivisionError, ValueError):
        pass
    raise ErosiveRootNotConvergedError(
        f"{ErosiveRootNotConvergedError.code}: the erosive-burning burn rate did not converge "
        f"in {max_iter} Newton steps (r0 = {r0:.9g} m/s, last r = {r:.9g} m/s).",
        context={"r0_m_s": r0, "last_r_m_s": r, "max_iter": max_iter, "rtol": EROSIVE_ROOT_RTOL})
