"""The isentrope of the chamber contents through a phase-transition plateau (A-6 stage 2, U-H25-5 group A-D; reports 272, 273).

Building blocks of a numerical isentropic expansion that does not assume a gamma_s constant along the path (H-25).  NOTHING in
this module is called by the energy run yet: no other module imports it (the connection is a later group, stage I).

Quantities are per kg of TOTAL products, as in gribs.energy.provider.  The condensed volume is 0 (V0), so v = n R_u T / p with
n = n_gas_per_kg and the enthalpy is h = u + n R_u T.  The provider must give the specific entropy (TPState.entropy_J_per_kgK);
without it every function here raises StateRecoveryError.

Plateau (a declared transition at T_m; the provider's two states at T_m -+ delta, delta = d_relative T_m): on the plateau T = T_m,
the phase fraction is phi = (s_0 - s_S)/(s_L - s_S), h = (1 - phi) h_S + phi h_L, n likewise, v = n R_u T_m / p, u = sqrt(2 (h_0 - h)).
The gas of the two sides is the same at T_m, so n(p) is a function of p alone on the plateau (a derivation, D2a/D2b of report 271;
NOT a result of RP-1311).  gamma_{s,T} = 1 / (1 - (d ln n / d ln p)_T) is formed from the number of kmol of gas n of the provider's
states at p exp(-+h); it is NOT read from TPState.gamma_s, so it does not pass the checks gamma_s > 1 of
energy_integrator.gamma_of and cea_tp.tp_detail (K-H25-38 (a), K-H25-40 (1)).  The states used are the provider's own states at
T_m -+ delta, which pass the check of cea_tp.tp_detail on THEIR gamma_s as every state of the provider does.

Tolerances: the roots are found by bracketing (brentq) with a relative tolerance ROOT_RTOL = 1e-12 (the temperature or the
pressure; the pressure is solved in ln p, whose absolute tolerance is the relative tolerance of p).
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Mapping, Optional, Tuple

from scipy.optimize import brentq

from gribs.energy.exceptions import StateRecoveryError
from gribs.energy.provider import EquilibriumProvider, TPState

#: Universal gas constant [J/(kmol K)] (the value of gribs.constants.R_UNIVERSAL; recovery.py has the same constant).
_RU = 8314.51
#: Relative tolerance of every root (temperature, and pressure through ln p).  Decision 110, item 5 (T1): at most 1e-12.
ROOT_RTOL = 1.0e-12
#: Rounding allowance of the check h <= h_0: a state with h above h_0 by at most _H_SLACK_REL (|h_0| + |h|) is the stagnation
#: state to rounding (u = 0, e.g. at p = p_0, where a root in T carries the root tolerance); a larger excess raises.
_H_SLACK_REL = 1.0e-9


def default_d_relative(root_tolerance: float) -> float:
    """The relative half-width of the plateau ends of recovery.plateau_bounds: max(1e-9, 10 root_tolerance) (recovery.py:163; a
    test compares the two)."""
    return max(1.0e-9, 10.0 * root_tolerance)


def _delta_K(T_m: float, d_relative: float) -> float:
    if not (math.isfinite(T_m) and T_m > 0.0):
        raise StateRecoveryError(f"T_m must be positive and finite, got {T_m!r}")
    if not (math.isfinite(d_relative) and d_relative > 0.0) or T_m * (1.0 - d_relative) <= 0.0:
        raise StateRecoveryError(f"d_relative must be positive, finite and below 1, got {d_relative!r}")
    return d_relative * T_m


def _tp(provider: EquilibriumProvider, T: float, p: float, B_e: Mapping[str, float]) -> TPState:
    if not (math.isfinite(p) and p > 0.0):
        raise StateRecoveryError(f"p must be positive and finite, got {p!r}")
    st = provider.tp(T, p, B_e)
    if not st.converged:
        raise StateRecoveryError(f"the equilibrium solve did not converge at T = {T:.10g} K, p = {p:.10g} Pa")
    return st


def _entropy(st: TPState, T: float, p: float) -> float:
    s = st.entropy_J_per_kgK
    if s is None or not math.isfinite(s):
        raise StateRecoveryError(f"the provider gave no finite entropy at T = {T:.10g} K, p = {p:.10g} Pa (got {s!r}); the "
                                 f"numerical isentrope needs TPState.entropy_J_per_kgK")
    return s


def _velocity(h_0: float, h: float, p: float) -> float:
    excess = h - h_0
    if excess > _H_SLACK_REL * (abs(h_0) + abs(h)):
        raise StateRecoveryError(f"h = {h!r} J/kg at p = {p:.10g} Pa is above the stagnation enthalpy h_0 = {h_0!r}")
    return math.sqrt(2.0 * max(h_0 - h, 0.0))


def _enthalpy(st: TPState, T: float) -> float:
    """h = u + n R_u T (V0: the condensed volume is 0, so p v = n R_u T)."""
    return st.u_J_per_kg + st.n_gas_per_kg * _RU * T


# ------------------------------------------------------------------------------------------------ stage B: gamma_{s,T}
@dataclass(frozen=True)
class PlateauGammaSide:
    """gamma_{s,T} from the provider's states at one side T_side = T_m -+ delta, at p exp(-h) and p exp(+h)."""
    T_K: float
    n_minus: float                  # n_gas_per_kg at p exp(-h)
    n_plus: float                   # n_gas_per_kg at p exp(+h)
    dlnn_dlnp: float                # central difference (ln n_plus - ln n_minus) / (2 h)
    gamma_sT: float                 # 1 / (1 - dlnn_dlnp)


@dataclass(frozen=True)
class PlateauGamma:
    """Both sides (T_m - delta, T_m + delta); no judgement of their agreement is made (decision 110, K-H25-47 (b))."""
    p_Pa: float
    dlnp: float
    lower: PlateauGammaSide
    upper: PlateauGammaSide


def _side_gamma(provider, B_e, T: float, p: float, dlnp: float) -> PlateauGammaSide:
    n_m = _tp(provider, T, p * math.exp(-dlnp), B_e).n_gas_per_kg
    n_p = _tp(provider, T, p * math.exp(dlnp), B_e).n_gas_per_kg
    if not (math.isfinite(n_m) and math.isfinite(n_p) and n_m > 0.0 and n_p > 0.0):
        raise StateRecoveryError(f"n of the gas must be positive and finite at T = {T:.10g} K, p = {p:.10g} Pa: {n_m!r}, {n_p!r}")
    d = (math.log(n_p) - math.log(n_m)) / (2.0 * dlnp)
    denom = 1.0 - d
    if not (math.isfinite(denom) and denom > 0.0):
        raise StateRecoveryError(f"1 - d ln n / d ln p = {denom!r} is not positive at T = {T:.10g} K, p = {p:.10g} Pa")
    return PlateauGammaSide(T_K=T, n_minus=n_m, n_plus=n_p, dlnn_dlnp=d, gamma_sT=1.0 / denom)


def plateau_gamma_sT(provider: EquilibriumProvider, B_e: Mapping[str, float], T_m: float, p: float, *, d_relative: float,
                     dlnp: float) -> PlateauGamma:
    """gamma_{s,T} on the plateau at pressure p from the provider's n at T_m -+ delta and p exp(-+dlnp) (4 TP calls).  dlnp (the
    step h of the central difference in ln p) has NO default (decision 110).  Reads TPState.n_gas_per_kg only."""
    delta = _delta_K(T_m, d_relative)
    if not (math.isfinite(dlnp) and dlnp > 0.0):
        raise StateRecoveryError(f"dlnp must be positive and finite, got {dlnp!r}")
    if not (math.isfinite(p) and p > 0.0):
        raise StateRecoveryError(f"p must be positive and finite, got {p!r}")
    return PlateauGamma(p_Pa=p, dlnp=dlnp, lower=_side_gamma(provider, B_e, T_m - delta, p, dlnp),
                        upper=_side_gamma(provider, B_e, T_m + delta, p, dlnp))


# ---------------------------------------------------------------------------------------------- stage C: plateau algebra
@dataclass(frozen=True)
class PlateauPoint:
    """The plateau state at pressure p for the stagnation entropy s_0 and enthalpy h_0.  region is "plateau" (0 <= phi <= 1),
    "low_side" (phi < 0: s_0 below the entropy of the low-T end, the path is on the low-T side) or "high_side" (phi > 1).
    h, v, u, n, G are None unless region == "plateau"."""
    p_Pa: float
    T_K: float
    region: str
    phi: float
    s_low: float
    s_high: float
    h_low: float
    h_high: float
    n_gas_per_kg: Optional[float] = None
    h_J_per_kg: Optional[float] = None
    v_m3_per_kg: Optional[float] = None
    u_m_per_s: Optional[float] = None
    G: Optional[float] = None       # u / v [kg/(s m^2)]


def plateau_point(p: float, s_0: float, h_0: float, B_e: Mapping[str, float], T_m: float, provider: EquilibriumProvider, *,
                  d_relative: float) -> PlateauPoint:
    """The state on the plateau at p (2 TP calls, at T_m -+ delta)."""
    delta = _delta_K(T_m, d_relative)
    for name, value in (("s_0", s_0), ("h_0", h_0)):
        if not math.isfinite(value):
            raise StateRecoveryError(f"{name} must be finite, got {value!r}")
    lo, hi = _tp(provider, T_m - delta, p, B_e), _tp(provider, T_m + delta, p, B_e)
    s_lo, s_hi = _entropy(lo, T_m - delta, p), _entropy(hi, T_m + delta, p)
    width = s_hi - s_lo
    if not width > 0.0:
        raise StateRecoveryError(f"the entropy does not increase across T_m at p = {p:.10g} Pa: s(T_m - delta) = {s_lo!r}, "
                                 f"s(T_m + delta) = {s_hi!r}")
    h_lo, h_hi = _enthalpy(lo, T_m - delta), _enthalpy(hi, T_m + delta)
    phi = (s_0 - s_lo) / width
    base = dict(p_Pa=p, T_K=T_m, phi=phi, s_low=s_lo, s_high=s_hi, h_low=h_lo, h_high=h_hi)
    if phi < 0.0:
        return PlateauPoint(region="low_side", **base)
    if phi > 1.0:
        return PlateauPoint(region="high_side", **base)
    n = (1.0 - phi) * lo.n_gas_per_kg + phi * hi.n_gas_per_kg
    h = (1.0 - phi) * h_lo + phi * h_hi
    v = n * _RU * T_m / p
    u = _velocity(h_0, h, p)
    return PlateauPoint(region="plateau", n_gas_per_kg=n, h_J_per_kg=h, v_m3_per_kg=v, u_m_per_s=u, G=u / v, **base)


# ------------------------------------------------------------------------------------------ stage D: single-phase roots
@dataclass(frozen=True)
class SinglePhasePoint:
    p_Pa: float
    T_K: float
    s_J_per_kgK: float
    n_gas_per_kg: float
    h_J_per_kg: float
    v_m3_per_kg: float
    u_m_per_s: float
    G: float


@dataclass(frozen=True)
class SinglePhaseResult:
    """status "root" (point given), "below_bracket" (s(T_lo, p) > s_0: the root is below T_lo) or "above_bracket"
    (s(T_hi, p) < s_0).  point is None unless status == "root"."""
    status: str
    point: Optional[SinglePhasePoint] = None


def single_phase_point(p: float, s_0: float, h_0: float, B_e: Mapping[str, float], provider: EquilibriumProvider,
                       T_bracket: Tuple[float, float]) -> SinglePhaseResult:
    """The state at p with s(T, p) = s_0 (a root in T inside T_bracket; one TP call per evaluation).  The provider's entropy
    increases with T at fixed p."""
    T_lo, T_hi = T_bracket
    if not (math.isfinite(T_lo) and math.isfinite(T_hi) and 0.0 < T_lo < T_hi):
        raise StateRecoveryError(f"T_bracket must satisfy 0 < T_lo < T_hi, got {T_bracket!r}")
    for name, value in (("s_0", s_0), ("h_0", h_0)):
        if not math.isfinite(value):
            raise StateRecoveryError(f"{name} must be finite, got {value!r}")
    cache = {}

    def state(T):
        if T not in cache:
            cache[T] = _tp(provider, T, p, B_e)
        return cache[T]

    def f(T):
        return _entropy(state(T), T, p) - s_0

    f_lo, f_hi = f(T_lo), f(T_hi)
    if f_lo > 0.0:
        return SinglePhaseResult("below_bracket")
    if f_hi < 0.0:
        return SinglePhaseResult("above_bracket")
    T = T_lo if f_lo == 0.0 else T_hi if f_hi == 0.0 else brentq(f, T_lo, T_hi, xtol=ROOT_RTOL * T_lo, rtol=ROOT_RTOL)
    st = state(T) if T in cache else _tp(provider, T, p, B_e)
    h = _enthalpy(st, T)
    v = st.n_gas_per_kg * _RU * T / p
    u = _velocity(h_0, h, p)
    return SinglePhaseResult("root", SinglePhasePoint(p_Pa=p, T_K=T, s_J_per_kgK=_entropy(st, T, p), n_gas_per_kg=st.n_gas_per_kg,
                                                       h_J_per_kg=h, v_m3_per_kg=v, u_m_per_s=u, G=u / v))


@dataclass(frozen=True)
class PressureRoot:
    """status "root" (p_Pa given), "above_bracket" (the root is above p_hi: s(p_hi) > s_0) or "below_bracket" (below p_lo)."""
    status: str
    p_Pa: Optional[float] = None


@dataclass(frozen=True)
class PlateauPressures:
    """p_high: the pressure at which the path reaches T_m from the high-T side (s(T_m + delta, p) = s_0); p_low: where it leaves the
    plateau to the low-T side (s(T_m - delta, p) = s_0).  p_low < p_high when the entropy of the high-T end is larger."""
    p_high: PressureRoot
    p_low: PressureRoot


def _pressure_root(provider, B_e, T: float, s_0: float, p_bracket: Tuple[float, float]) -> PressureRoot:
    p_lo, p_hi = p_bracket

    def g(x):                                               # x = ln p; s decreases with p at fixed T
        p = math.exp(x)
        return _entropy(_tp(provider, T, p, B_e), T, p) - s_0

    x_lo, x_hi = math.log(p_lo), math.log(p_hi)
    g_lo, g_hi = g(x_lo), g(x_hi)
    if g_hi > 0.0:
        return PressureRoot("above_bracket")
    if g_lo < 0.0:
        return PressureRoot("below_bracket")
    if g_lo == 0.0:
        return PressureRoot("root", p_lo)
    if g_hi == 0.0:
        return PressureRoot("root", p_hi)
    return PressureRoot("root", math.exp(brentq(g, x_lo, x_hi, xtol=ROOT_RTOL, rtol=ROOT_RTOL)))


def plateau_pressures(s_0: float, B_e: Mapping[str, float], T_m: float, provider: EquilibriumProvider,
                      p_bracket: Tuple[float, float], *, d_relative: float) -> PlateauPressures:
    """The two pressures at which the isentrope of entropy s_0 meets the plateau ends (roots in ln p inside p_bracket)."""
    delta = _delta_K(T_m, d_relative)
    p_lo, p_hi = p_bracket
    if not (math.isfinite(p_lo) and math.isfinite(p_hi) and 0.0 < p_lo < p_hi):
        raise StateRecoveryError(f"p_bracket must satisfy 0 < p_lo < p_hi, got {p_bracket!r}")
    if not math.isfinite(s_0):
        raise StateRecoveryError(f"s_0 must be finite, got {s_0!r}")
    return PlateauPressures(p_high=_pressure_root(provider, B_e, T_m + delta, s_0, p_bracket),
                            p_low=_pressure_root(provider, B_e, T_m - delta, s_0, p_bracket))
