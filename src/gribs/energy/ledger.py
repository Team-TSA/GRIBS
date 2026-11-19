"""Closure residuals of an energy-mode run (A-6 stage 2, C10; contract 164 sections 1b, 2b-4.1, 2b-4.2).

Four residuals, each reported twice, as an ABSOLUTE value and normalised by a SCALE of the quantity (the two are kept apart):

R_M  mass:     M - (M0 + m_gen - m_out)                    the mass of the contents against the generated and discharged mass
R_B  elements: sum_e B_e - M                               the element inventories add up to the mass
R_E  energy:   E_cv - (E0 + Q + W_burn - H_out)            the energy of the contents against the accumulated energy terms
               (Q the wall heat, W_burn the net burn-surface term, H_out the energy carried out by the discharge; each
               accumulated by its own state variable, not by integrating E_cv)
R_V  volume:   V_gas + V_condensed - V_free                the volumes of the recovered state fill the free volume

check_closure raises ClosureError for a residual above its tolerance.  Standard library only.
"""
import math
from dataclasses import dataclass
from typing import Mapping

from gribs.errors import GribsError


class ClosureError(GribsError):
    """A closure residual of an energy-mode run exceeds its tolerance (contract 2b-4.1: no silent loss of mass, elements,
    energy or volume)."""


@dataclass(frozen=True)
class ClosureResiduals:
    R_M: float
    R_B: float
    R_E: float
    R_V: float
    R_M_normalised: float
    R_B_normalised: float
    R_E_normalised: float
    R_V_normalised: float


def _scale(*values: float) -> float:
    return max(max(abs(v) for v in values), 1.0e-300)


def closure_residuals(*, M: float, M0: float, m_gen: float, m_out: float, B_e: Mapping[str, float],
                      E_cv: float, E0: float, Q_cum: float, W_burn_cum: float, H_out_cum: float,
                      V_gas: float, V_condensed: float, V_free: float) -> ClosureResiduals:
    """The four closure residuals at one time (the quantities are the state of the run and its accumulators)."""
    values = dict(M=M, M0=M0, m_gen=m_gen, m_out=m_out, E_cv=E_cv, E0=E0, Q_cum=Q_cum, W_burn_cum=W_burn_cum,
                  H_out_cum=H_out_cum, V_gas=V_gas, V_condensed=V_condensed, V_free=V_free)
    for name, value in values.items():
        if not math.isfinite(value):
            raise ClosureError(f"{name} must be finite, got {value!r}")
    if not B_e or any(not (math.isfinite(v) and v >= 0.0) for v in B_e.values()):
        raise ClosureError("B_e must be a non-empty table of finite, non-negative inventories")
    R_M = M - (M0 + m_gen - m_out)
    R_B = math.fsum(B_e.values()) - M
    R_E = E_cv - (E0 + Q_cum + W_burn_cum - H_out_cum)
    R_V = V_gas + V_condensed - V_free
    s_M = _scale(M, M0, m_gen, m_out)
    s_E = _scale(E_cv, E0, Q_cum, W_burn_cum, H_out_cum)
    s_V = _scale(V_free)
    return ClosureResiduals(R_M=R_M, R_B=R_B, R_E=R_E, R_V=R_V, R_M_normalised=R_M / s_M, R_B_normalised=R_B / s_M,
                            R_E_normalised=R_E / s_E, R_V_normalised=R_V / s_V)


def check_closure(r: ClosureResiduals, *, tolerance: float) -> None:
    """Raise ClosureError unless every NORMALISED residual is within tolerance (a positive, finite number: the caller
    states it, nothing here defaults it)."""
    if not (math.isfinite(tolerance) and tolerance > 0.0):
        raise ClosureError(f"the closure tolerance must be positive and finite, got {tolerance!r}")
    for name in ("R_M", "R_B", "R_E", "R_V"):
        value = getattr(r, name + "_normalised")
        if not abs(value) <= tolerance:
            raise ClosureError(f"closure residual {name} = {getattr(r, name)!r} (normalised {value!r}) exceeds {tolerance:g}")
