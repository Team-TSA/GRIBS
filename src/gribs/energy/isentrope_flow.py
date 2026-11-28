"""The flow along the isentrope of the chamber contents: throat, exit state, thrust (A-6 stage 2, U-H25-5 group E-G; reports 277, 278).

Built on gribs.energy.isentrope (stages A-D).  The energy run calls it through gribs.energy.numerical_path (group H-I, which decides when
the numerical flow replaces the closed form).  Quantities are per kg of TOTAL products; the condensed volume is 0 (V0).

The mass flux G(p) = u / v = u p / (n R_u T) of the isentrope (u = sqrt(2 (h_0 - h)), the stagnation enthalpy h_0) is a function of
the pressure alone once the stagnation state is fixed.  By continuity G A = mdot, so the area ratio at p is A / A_t = G* / G(p) with
G* = G(p_t) at the throat.  The throat is the FIRST maximum of G met from p_0 downward (K-H25-32 (c)); the speed of sound is not used
to find it.  The first maximum is searched on a geometric grid and refined; it may be a smooth maximum (u = a) or a kink at the entry
of a phase-transition plateau (an inequality a_{s,T} <= u <= a_I at p_high, derived only; D2a/D2b, report 271, NOT RP-1311).

No default is given to any numerical setting (grid points, search tolerances, the probe step of the kink test): every one is a
required argument (decision 113, K-H25-52).  A state that cannot be recovered is returned as a classification, never filled in.
The momentum equation F = eta mdot u_e + (p_e - p_a) A_e does not depend on gamma.  The discharge coefficient is not applied (mdot
here is the ideal one); the separation of the flow is NOT computed (decision 113, K-H25-41 (c): the thrust is the one at A_e).
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Mapping, Optional, Tuple

from scipy.optimize import brentq, minimize_scalar

from gribs.energy.exceptions import StateRecoveryError
from gribs.energy.isentrope import (ROOT_RTOL, plateau_gamma_sT, plateau_point, plateau_pressures,
                                    single_phase_point)
from gribs.energy.provider import EquilibriumProvider

#: Step of the central difference in ln p used to polish a maximum (internal; the polished value is accepted only if G does not
#: become smaller than the value of the first search).  A property of the polish, not of the physics.
_POLISH_DLNP = 1.0e-5
#: Half-width (in ln p) of the interval of the polish around the first maximum.
_POLISH_HALF_WIDTH = 1.0e-4
#: Relative allowance by which G of the polished maximum may be below the value of the first search and still be accepted.
_POLISH_ACCEPT_REL = 1.0e-12
#: An area ratio at or below 1 + this is taken as 1 (the same threshold as gribs.energy.nozzle_from_state).
_EPS_ONE = 1.0e-12
#: Relative slack above p_0 for a requested pressure (rounding).
_P0_SLACK = 1.0e-12


class _Failed(Exception):
    """A state of the path could not be recovered (internal: carries the classification text)."""


@dataclass(frozen=True)
class FlowState:
    """The state of the isentrope at pressure p.  status is "ok" (the fields are given), "above_stagnation" (p > p_0),
    "bracket_below" / "bracket_above" (the root in T is outside T_bracket), "plateau_phase_fraction_outside" (p inside the plateau
    interval but the phase fraction is outside [0, 1]) or "state_error" (the provider or a check raised StateRecoveryError).
    region is "single_phase" or "plateau"; phi is the phase fraction on the plateau."""
    status: str
    p_Pa: float
    region: Optional[str] = None
    T_K: Optional[float] = None
    n_gas_per_kg: Optional[float] = None
    h_J_per_kg: Optional[float] = None
    v_m3_per_kg: Optional[float] = None
    u_m_per_s: Optional[float] = None
    G: Optional[float] = None               # u / v [kg/(s m^2)]
    phi: Optional[float] = None
    message: str = ""


class Expansion:
    """The isentropic expansion from the stagnation state (s_0, h_0, p_0) of the chamber contents: state_at(p) for any p in
    [p_floor, p_0].  T_bracket must hold the root in T of s(T, p) = s_0 for every single-phase p of interest.  With T_m (the
    declared transition; d_relative required) the plateau interval (p_low, p_high) is found once from the roots of the entropies of
    the two ends (plateau_pressures, p_bracket = (p_floor, p_0)); p inside the open interval is on the plateau (plateau_point), every
    other p is single phase (single_phase_point; p equal to an end is single phase).  The states are cached by p."""

    def __init__(self, provider: EquilibriumProvider, B_e: Mapping[str, float], s_0: float, h_0: float, p_0: float, p_floor: float,
                 T_bracket: Tuple[float, float], *, T_m: Optional[float] = None, d_relative: Optional[float] = None):
        for name, value in (("s_0", s_0), ("h_0", h_0)):
            if not math.isfinite(value):
                raise StateRecoveryError(f"{name} must be finite, got {value!r}")
        if not (math.isfinite(p_0) and math.isfinite(p_floor) and 0.0 < p_floor < p_0):
            raise StateRecoveryError(f"0 < p_floor < p_0 is required, got p_floor = {p_floor!r}, p_0 = {p_0!r}")
        T_lo, T_hi = T_bracket
        if not (math.isfinite(T_lo) and math.isfinite(T_hi) and 0.0 < T_lo < T_hi):
            raise StateRecoveryError(f"T_bracket must satisfy 0 < T_lo < T_hi, got {T_bracket!r}")
        if T_m is not None and d_relative is None:
            raise StateRecoveryError("d_relative is required when T_m is given")
        self.provider, self.B_e, self.s_0, self.h_0, self.p_0, self.p_floor = provider, B_e, s_0, h_0, p_0, p_floor
        self.T_bracket, self.T_m, self.d_relative = (T_lo, T_hi), T_m, d_relative
        self._cache: dict = {}
        #: "none" (no transition declared), "absent_above" (the path stays above the plateau in [p_floor, p_0]),
        #: "absent_below" (below it), or "present".
        self.plateau = "none"
        self.p_high = math.inf              # the plateau interval (p_low, p_high); meaningful when plateau == "present"
        self.p_low = 0.0
        if T_m is not None:
            pp = plateau_pressures(s_0, B_e, T_m, provider, (p_floor, p_0), d_relative=d_relative)
            hi, lo = pp.p_high, pp.p_low
            if hi.status == "below_bracket":
                self.plateau = "absent_above"
            elif lo.status == "above_bracket":
                self.plateau = "absent_below"
            else:
                self.plateau = "present"
                self.p_high = math.inf if hi.status == "above_bracket" else hi.p_Pa
                self.p_low = 0.0 if lo.status == "below_bracket" else lo.p_Pa
                if not self.p_low < self.p_high:
                    raise StateRecoveryError(f"the plateau pressures are not ordered: p_low = {self.p_low!r}, p_high = {self.p_high!r}")

    @property
    def has_plateau(self) -> bool:
        return self.plateau == "present"

    def _on_plateau(self, p: float) -> bool:
        return self.has_plateau and self.p_low < p < self.p_high

    def state_at(self, p: float) -> FlowState:
        """The state at p (cached).  A state that cannot be recovered is returned with a status, not raised."""
        if not (math.isfinite(p) and p > 0.0):
            raise StateRecoveryError(f"p must be positive and finite, got {p!r}")
        if p not in self._cache:
            self._cache[p] = self._compute(p)
        return self._cache[p]

    def _compute(self, p: float) -> FlowState:
        if p > self.p_0 * (1.0 + _P0_SLACK):
            return FlowState("above_stagnation", p, message=f"p = {p!r} is above p_0 = {self.p_0!r}")
        try:
            if self._on_plateau(p):
                pt = plateau_point(p, self.s_0, self.h_0, self.B_e, self.T_m, self.provider, d_relative=self.d_relative)
                if pt.region != "plateau":
                    return FlowState("plateau_phase_fraction_outside", p, phi=pt.phi, message=f"phase fraction {pt.phi!r} at p = {p!r}")
                return FlowState("ok", p, "plateau", pt.T_K, pt.n_gas_per_kg, pt.h_J_per_kg, pt.v_m3_per_kg, pt.u_m_per_s, pt.G, pt.phi)
            res = single_phase_point(p, self.s_0, self.h_0, self.B_e, self.provider, self.T_bracket)
        except StateRecoveryError as exc:
            return FlowState("state_error", p, message=str(exc))
        if res.status == "below_bracket":
            return FlowState("bracket_below", p, message=f"the root in T at p = {p!r} is below T_bracket")
        if res.status == "above_bracket":
            return FlowState("bracket_above", p, message=f"the root in T at p = {p!r} is above T_bracket")
        q = res.point
        return FlowState("ok", p, "single_phase", q.T_K, q.n_gas_per_kg, q.h_J_per_kg, q.v_m3_per_kg, q.u_m_per_s, q.G)

    def _G(self, p: float) -> float:
        if p >= self.p_0:
            return 0.0                       # the stagnation state: u = 0
        st = self.state_at(p)
        if st.status != "ok":
            raise _Failed(f"{st.status}: {st.message}")
        return st.G


def _check_positive(name: str, value, integer: bool = False, minimum: float = 0.0) -> None:
    ok = math.isfinite(value) and value > minimum and (not integer or float(value) == int(value))
    if not ok:
        raise StateRecoveryError(f"{name} must be {'an integer ' if integer else ''}finite and above {minimum}, got {value!r}")


def _grid_x(x_start: float, x_end: float, n_grid: int):
    return [x_start + (x_end - x_start) * i / (n_grid - 1) for i in range(n_grid)]


# ------------------------------------------------------------------------------------------------ stage E: the throat
@dataclass(frozen=True)
class Throat:
    """classification: "single_phase_high" (p_t above the plateau, with the plateau below), "kink_at_plateau_entry" (G has its
    maximum at p_high), "inside_plateau", "single_phase_low" (p_t below p_low), "single_phase" (no transition declared),
    "single_phase_above_plateau" / "single_phase_below_plateau" (a transition is declared but the path does not meet the plateau in
    [p_floor, p_0]), "not_found" (no maximum on the grid) or "state_error" (a state needed could not be recovered).
    G_star, p_t_Pa and state are None for the last two."""
    classification: str
    G_star: Optional[float] = None
    p_t_Pa: Optional[float] = None
    state: Optional[FlowState] = None
    n_grid: int = 0
    message: str = ""


def _polish_maximum(exp: Expansion, x_m: float, G_m: float, lo: float, hi: float) -> Tuple[float, float]:
    """Polish a maximum of G in x = ln p: a root of the central difference of G in x near x_m.  Falls back to (x_m, G_m)."""
    a, b = max(lo, x_m - _POLISH_HALF_WIDTH), min(hi, x_m + _POLISH_HALF_WIDTH)
    h = _POLISH_DLNP

    def d(x):
        return exp._G(math.exp(x + h)) - exp._G(math.exp(x - h))

    try:
        if not (a - h > math.log(exp.p_floor) - 1.0e-9 and b + h <= math.log(exp.p_0)):
            return x_m, G_m
        da, db = d(a), d(b)
        if not (da * db < 0.0):
            return x_m, G_m
        x = brentq(d, a, b, xtol=1.0e-13, rtol=1.0e-14, maxiter=200)
        g = exp._G(math.exp(x))
    except (_Failed, ValueError, RuntimeError):
        return x_m, G_m
    return (x, g) if g >= G_m * (1.0 - _POLISH_ACCEPT_REL) else (x_m, G_m)


def find_throat(exp: Expansion, *, n_grid: int, kink_dlnp: float, ln_p_xtol: float) -> Throat:
    """The first maximum of G from p_0 downward.  n_grid points, equally spaced in ln p from p_0 to p_floor (G(p_0) = 0 is not
    evaluated); the first grid point i with G_i > G_{i-1} and G_i >= G_{i+1} brackets it in [p_{i+1}, p_{i-1}].  If the plateau
    entry p_high is inside that interval the maximum is a kink at p_high when G(p_high) exceeds G at p_high exp(+-kink_dlnp);
    otherwise the maximum is refined by a bounded search in ln p (absolute tolerance ln_p_xtol) and polished.  kink_dlnp must keep
    both probe points inside the bracketing interval.  All three arguments are required (no default)."""
    _check_positive("n_grid", n_grid, integer=True, minimum=3.0)
    _check_positive("kink_dlnp", kink_dlnp)
    _check_positive("ln_p_xtol", ln_p_xtol)
    x0, x1 = math.log(exp.p_0), math.log(exp.p_floor)
    xs = _grid_x(x0, x1, n_grid)
    ps = [exp.p_0] + [math.exp(x) for x in xs[1:]]
    Gs = {0: 0.0}
    try:
        def g(i):
            if i not in Gs:
                Gs[i] = exp._G(ps[i])
            return Gs[i]
        found = None
        for i in range(1, n_grid - 1):
            if g(i) > g(i - 1) and g(i) >= g(i + 1):
                found = i
                break
        if found is None:
            return Throat("not_found", n_grid=n_grid, message="no maximum of G on the grid between p_0 and p_floor")
        i = found
        lo, hi = xs[i + 1], xs[i - 1]                                   # x decreases with i
        if exp.has_plateau and math.isfinite(exp.p_high) and lo < math.log(exp.p_high) < hi:
            x_h = math.log(exp.p_high)
            if not (lo < x_h - kink_dlnp and x_h + kink_dlnp < hi):
                raise StateRecoveryError(f"kink_dlnp = {kink_dlnp!r} puts a probe point of the kink test outside the interval of the "
                                         f"first maximum [{math.exp(lo)!r}, {math.exp(hi)!r}] Pa")
            G_h = exp._G(exp.p_high)
            if G_h > exp._G(math.exp(x_h + kink_dlnp)) and G_h > exp._G(math.exp(x_h - kink_dlnp)):
                st = exp.state_at(exp.p_high)
                return Throat("kink_at_plateau_entry", G_h, exp.p_high, st, n_grid)
        m = minimize_scalar(lambda x: -exp._G(math.exp(x)), bounds=(lo, hi), method="bounded", options={"xatol": ln_p_xtol})
        x_m, G_m = (m.x, -m.fun) if -m.fun > Gs[i] else (xs[i], Gs[i])
        x_m, G_m = _polish_maximum(exp, x_m, G_m, lo, hi)
        p_t = exp.p_0 if x_m >= x0 else math.exp(x_m)
        st = exp.state_at(p_t) if p_t < exp.p_0 else None
    except _Failed as exc:
        return Throat("state_error", n_grid=n_grid, message=str(exc))
    if st is None or st.status != "ok":
        return Throat("state_error", n_grid=n_grid, message="the state at the throat could not be recovered")
    return Throat(_throat_class(exp, p_t), G_m, p_t, st, n_grid)


def _throat_class(exp: Expansion, p_t: float) -> str:
    if exp.plateau == "none":
        return "single_phase"
    if exp.plateau == "absent_above":
        return "single_phase_above_plateau"
    if exp.plateau == "absent_below":
        return "single_phase_below_plateau"
    if exp.p_low < p_t < exp.p_high:
        return "inside_plateau"
    return "single_phase_high" if p_t >= exp.p_high else "single_phase_low"


# ------------------------------------------------------------------------------------ stage F: exit state, choke, thrust
@dataclass(frozen=True)
class ExitState:
    """classification: "unchoked_subsonic" (p_e = p_a), "choked_supersonic", "choked_sonic" (A_e = A_t), "no_flow" (p_a >= p_0),
    "out_of_range" (no continuous isentropic solution: G stops decreasing, a fold, before G* / G reaches the area ratio; stage G,
    K-H25-36), "not_found" (the grid does not reach the area ratio), "no_throat" (the throat was not found; message says why) or
    "state_error".  The fields are None when no state is given.  p_sub_Pa: the exit pressure of the subsonic solution of the
    area ratio (the numerical choke threshold: choked when p_a < p_sub_Pa).  mdot is the ideal flow (no discharge coefficient)."""
    classification: str
    choked: Optional[bool] = None
    p_e_Pa: Optional[float] = None
    u_e_m_per_s: Optional[float] = None
    T_e_K: Optional[float] = None
    mdot_kg_per_s: Optional[float] = None
    F_N: Optional[float] = None
    p_sub_Pa: Optional[float] = None
    state: Optional[FlowState] = None
    fold_p_Pa: Optional[float] = None
    message: str = ""


def _root_ratio(exp: Expansion, G_star: float, eps: float, x_a: float, x_b: float) -> float:
    """The pressure in [exp(x_b), exp(x_a)] with G_star / G(p) = eps (root in ln p, the ratio is monotone in the bracket)."""
    def f(x):
        G = exp._G(math.exp(x))
        return (G_star / G if G > 0.0 else math.inf) - eps
    return math.exp(brentq(f, min(x_a, x_b), max(x_a, x_b), xtol=ROOT_RTOL, rtol=ROOT_RTOL, maxiter=300))


def _subsonic_pressure(exp: Expansion, throat: Throat, eps: float, n_grid: int) -> float:
    """The root p_sub in (p_t, p_0) of G* / G = eps, the first crossing from p_0 downward on n_grid points."""
    x0, xt = math.log(exp.p_0), math.log(throat.p_t_Pa)
    xs = _grid_x(x0, xt, n_grid)
    prev_x = xs[0]
    for x in xs[1:]:
        G = exp._G(math.exp(x)) if x > xt else throat.G_star
        if G > 0.0 and throat.G_star / G <= eps:
            return _root_ratio(exp, throat.G_star, eps, prev_x, x)
        prev_x = x
    # practically unreachable: G* / G >= 1 at the throat and the scan from p_0 passes the throat pressure (report 390 s.4)
    raise _Failed("not_found: the subsonic branch does not reach the area ratio on the grid")


def _supersonic_pressure(exp: Expansion, throat: Throat, eps: float, n_grid: int) -> Tuple[Optional[float], Optional[float]]:
    """(p_e, None) the first root below p_t of G* / G = eps; (None, p_fold) if G stops decreasing before it (a fold);
    raises _Failed("not_found...") if the grid ends first."""
    xt, x1 = math.log(throat.p_t_Pa), math.log(exp.p_floor)
    xs = _grid_x(xt, x1, n_grid)
    prev_x, prev_G = xs[0], throat.G_star
    for x in xs[1:]:
        G = exp._G(math.exp(x))
        if throat.G_star / G >= eps:
            return _root_ratio(exp, throat.G_star, eps, prev_x, x), None
        if G > prev_G:
            return None, math.exp(prev_x)
        prev_x, prev_G = x, G
    note = f"; the grid ends at p_floor = {exp.p_floor!r} Pa, the requested area ratio is {eps!r}"
    try:                                                    # informative only; the not_found itself is raised below
        st = exp.state_at(exp.p_floor)
        if st.status == "ok":
            note += f"; last grid point: p = {exp.p_floor!r} Pa, T = {st.T_K!r} K, area ratio = {throat.G_star / exp._G(exp.p_floor)!r}"
    except Exception:
        pass
    raise _Failed("not_found: the supersonic branch does not reach the area ratio on the grid" + note)


def exit_state(exp: Expansion, throat: Throat, *, eps: float, p_a: float, eta: float, A_t: float, n_grid: int) -> ExitState:
    """The exit state of a nozzle of area ratio eps = A_e / A_t into the ambient pressure p_a (throat area A_t, thrust efficiency
    eta).  The mass flow is set at the throat, mdot = G* A_t, when the nozzle is choked (p_a < p_sub); otherwise the exit pressure is
    p_a and mdot = G(p_a) A_e.  F = eta mdot u_e + (p_e - p_a) A_e (p_e = p_a when unchoked).  Flow separation is NOT computed.
    n_grid (required): points of the branch scans."""
    _check_positive("n_grid", n_grid, integer=True, minimum=2.0)
    for name, value in (("eps", eps), ("p_a", p_a), ("eta", eta), ("A_t", A_t)):
        _check_positive(name, value)
    if eps < 1.0 - _EPS_ONE:
        raise StateRecoveryError(f"the area ratio must be at least 1, got {eps!r}")
    if throat.classification in ("not_found", "state_error"):
        return ExitState("no_throat", message=f"{throat.classification}: {throat.message}")
    A_e = eps * A_t
    if p_a >= exp.p_0:
        return ExitState("no_flow", choked=False, p_e_Pa=p_a, u_e_m_per_s=0.0, mdot_kg_per_s=0.0, F_N=0.0)
    try:
        sonic = eps <= 1.0 + _EPS_ONE
        p_sub = throat.p_t_Pa if sonic else _subsonic_pressure(exp, throat, eps, n_grid)
        if p_a >= p_sub:                                                # not choked: p_e = p_a, mdot = G(p_a) A_e
            st = exp.state_at(p_a)
            if st.status != "ok":
                raise _Failed(f"{st.status}: {st.message}")
            mdot = st.G * A_e
            return ExitState("unchoked_subsonic", False, p_a, st.u_m_per_s, st.T_K, mdot, eta * mdot * st.u_m_per_s, p_sub, st)
        if sonic:
            p_e, fold = throat.p_t_Pa, None
        else:
            p_e, fold = _supersonic_pressure(exp, throat, eps, n_grid)
        if p_e is None:
            return ExitState("out_of_range", True, p_sub_Pa=p_sub, fold_p_Pa=fold,
                             message="G stops decreasing before G* / G reaches the area ratio: no continuous isentropic solution")
        st = exp.state_at(p_e)
        if st.status != "ok":
            raise _Failed(f"{st.status}: {st.message}")
        mdot = throat.G_star * A_t
        F = eta * mdot * st.u_m_per_s + (p_e - p_a) * A_e
        return ExitState("choked_sonic" if sonic else "choked_supersonic", True, p_e, st.u_m_per_s, st.T_K, mdot, F, p_sub, st)
    except _Failed as exc:
        text = str(exc)
        return ExitState("not_found" if text.startswith("not_found") else "state_error", message=text)


# ------------------------------------------------------------------------------ stage G: the limit of the supersonic branch
@dataclass(frozen=True)
class BranchLimit:
    """The fold of the supersonic branch: the first minimum of G below p_t.  classification: "kink_at_plateau_exit" (the minimum is at
    p_low, G at p_low exp(+-kink_dlnp) both above), "smooth_minimum", "none_in_range" (G decreases to p_floor: every area ratio up to
    G* / G(p_floor) has a solution), "no_throat" or "state_error".  eps_max = G* / G_min is the largest area ratio with a solution
    on the branch (None for the last two and, for "none_in_range", the value at p_floor is given in eps_max as a lower bound)."""
    classification: str
    p_Pa: Optional[float] = None
    G_min: Optional[float] = None
    eps_max: Optional[float] = None
    message: str = ""


def supersonic_limit(exp: Expansion, throat: Throat, *, n_grid: int, kink_dlnp: float, ln_p_xtol: float) -> BranchLimit:
    """Detect the fold (stage G; K-H25-36).  Detection and classification only: no exception class, no warning is defined here."""
    _check_positive("n_grid", n_grid, integer=True, minimum=3.0)
    _check_positive("kink_dlnp", kink_dlnp)
    _check_positive("ln_p_xtol", ln_p_xtol)
    if throat.classification in ("not_found", "state_error"):
        return BranchLimit("no_throat", message=f"{throat.classification}: {throat.message}")
    xt, x1 = math.log(throat.p_t_Pa), math.log(exp.p_floor)
    xs = _grid_x(xt, x1, n_grid)
    try:
        Gs = [throat.G_star]
        fold_i = None
        for i in range(1, n_grid):
            Gs.append(exp._G(math.exp(xs[i])))
            if Gs[i] > Gs[i - 1]:
                fold_i = i
                break
        if fold_i is None:
            return BranchLimit("none_in_range", exp.p_floor, Gs[-1], throat.G_star / Gs[-1])
        i = fold_i
        lo, hi = xs[i], xs[i - 2] if i >= 2 else xs[0]                   # G_{i-1} is the smallest so far: bracket [x_i, x_{i-2}]
        x_c = xs[i - 1]
        edge = math.log(exp.p_low) if exp.has_plateau and exp.p_low > 0.0 else None
        if edge is not None and lo < edge < hi:
            if not (lo < edge - kink_dlnp and edge + kink_dlnp < hi):
                raise StateRecoveryError(f"kink_dlnp = {kink_dlnp!r} puts a probe point of the kink test outside the interval of the "
                                         f"fold [{math.exp(lo)!r}, {math.exp(hi)!r}] Pa")
            G_e = exp._G(exp.p_low)
            if G_e < exp._G(math.exp(edge + kink_dlnp)) and G_e < exp._G(math.exp(edge - kink_dlnp)):
                return BranchLimit("kink_at_plateau_exit", exp.p_low, G_e, throat.G_star / G_e)
        m = minimize_scalar(lambda x: exp._G(math.exp(x)), bounds=(lo, hi), method="bounded", options={"xatol": ln_p_xtol})
        G_m = min(m.fun, Gs[i - 1])
        p_m = math.exp(m.x) if m.fun <= Gs[i - 1] else math.exp(x_c)
        return BranchLimit("smooth_minimum", p_m, G_m, throat.G_star / G_m)
    except _Failed as exc:
        return BranchLimit("state_error", message=str(exc))


# ------------------------------------------------------------------------- stage H: continuation from the previous evaluation
def find_throat_near(exp: Expansion, *, p_hint: float, width: float, ln_p_xtol: float, kink_dlnp: float) -> Optional[Throat]:
    """The throat searched only in [p_hint exp(-width), p_hint exp(+width)] (ln p; clipped to [p_floor, p_0]), a continuation of the
    previous evaluation (decision 117, K-H25-57 (2)).  Returns None, and the caller must use find_throat (the cold search), when the
    interval does not bracket a maximum (G at either end is not below the maximum found), when the kink test cannot be made inside it,
    or when a state of the interval cannot be recovered.  The classification is NOT compared here (the caller compares it with the
    previous one).  The maximum is the one of the bracket, taken to be the first from p_0 because the bracket is around the previous
    first maximum; this is an assumption of the continuation, and find_throat is the check (a test compares the two)."""
    _check_positive("width", width)
    _check_positive("kink_dlnp", kink_dlnp)
    _check_positive("ln_p_xtol", ln_p_xtol)
    if not (math.isfinite(p_hint) and exp.p_floor < p_hint < exp.p_0):
        return None
    x0, xf, xh0 = math.log(exp.p_0), math.log(exp.p_floor), math.log(p_hint)
    lo, hi = max(xh0 - width, xf), min(xh0 + width, x0)
    try:
        G_lo, G_hi = exp._G(math.exp(lo)), exp._G(math.exp(hi))
        if exp.has_plateau and math.isfinite(exp.p_high) and lo < math.log(exp.p_high) < hi:
            x_h = math.log(exp.p_high)
            if not (lo < x_h - kink_dlnp and x_h + kink_dlnp < hi):
                return None
            G_h = exp._G(exp.p_high)
            if G_h > exp._G(math.exp(x_h + kink_dlnp)) and G_h > exp._G(math.exp(x_h - kink_dlnp)) and G_h > G_lo and G_h > G_hi:
                st = exp.state_at(exp.p_high)
                return Throat("kink_at_plateau_entry", G_h, exp.p_high, st, 0)
        m = minimize_scalar(lambda x: -exp._G(math.exp(x)), bounds=(lo, hi), method="bounded", options={"xatol": ln_p_xtol})
        x_m, G_m = m.x, -m.fun
        if not (G_m > G_lo and G_m > G_hi):
            return None
        x_m, G_m = _polish_maximum(exp, x_m, G_m, lo, hi)
        p_t = math.exp(x_m)
        st = exp.state_at(p_t)
    except _Failed:
        return None
    if st.status != "ok":
        return None
    return Throat(_throat_class(exp, p_t), G_m, p_t, st, 0)


def exit_state_near(exp: Expansion, throat: Throat, *, eps: float, p_a: float, eta: float, A_t: float, p_exit_hint: Optional[float],
                    width: float) -> Optional[ExitState]:
    """The exit state by a continuation (decision 117, K-H25-57 (2)): the choke test from G(p_a) alone (choked when p_a <= p_t or
    G* / G(p_a) < eps; this is the test of exit_state, p_a < p_sub, written with the monotone G on the subsonic branch; p_sub_Pa is
    None) and the exit pressure as the root of G* / G = eps in [p_hint exp(-width), p_hint exp(+width)] below p_t.  Returns None, and the
    caller must use exit_state (the cold scan), when the throat is not usable, the area ratio is 1 or p_a >= p_0 (those are cheap in
    exit_state), the interval does not bracket the root, or a state could not be recovered.  The branch between p_t and the root is
    taken to be monotone (no fold): a fold makes the root disappear from the interval or is found by exit_state."""
    _check_positive("width", width)
    for name, value in (("eps", eps), ("p_a", p_a), ("eta", eta), ("A_t", A_t)):
        _check_positive(name, value)
    if eps < 1.0 - _EPS_ONE:
        raise StateRecoveryError(f"the area ratio must be at least 1, got {eps!r}")
    if throat.classification in ("not_found", "state_error") or p_a >= exp.p_0 or eps <= 1.0 + _EPS_ONE:
        return None
    A_e = eps * A_t
    try:
        if p_a > throat.p_t_Pa:
            st = exp.state_at(p_a)
            if st.status != "ok":
                return None
            if throat.G_star / st.G >= eps:                             # not choked: p_e = p_a, mdot = G(p_a) A_e
                mdot = st.G * A_e
                return ExitState("unchoked_subsonic", False, p_a, st.u_m_per_s, st.T_K, mdot, eta * mdot * st.u_m_per_s, None, st)
        if p_exit_hint is None or not (math.isfinite(p_exit_hint) and exp.p_floor < p_exit_hint < throat.p_t_Pa):
            return None
        xt, xh0 = math.log(throat.p_t_Pa), math.log(p_exit_hint)
        xa, xb = max(xh0 - width, math.log(exp.p_floor)), min(xh0 + width, xt)

        def f(x):
            G = throat.G_star if x >= xt else exp._G(math.exp(x))
            return (throat.G_star / G if G > 0.0 else math.inf) - eps
        if not (f(xa) >= 0.0 >= f(xb)):
            return None
        p_e = _root_ratio(exp, throat.G_star, eps, xa, xb)
        st = exp.state_at(p_e)
        if st.status != "ok":
            return None
        mdot = throat.G_star * A_t
        return ExitState("choked_supersonic", True, p_e, st.u_m_per_s, st.T_K, mdot, eta * mdot * st.u_m_per_s + (p_e - p_a) * A_e, None, st)
    except (_Failed, ValueError):
        return None


# --------------------------------------------------------------------------------------------------- the speed of sound
def sound_speed_single_phase(exp: Expansion, state: FlowState) -> float:
    """a = sqrt(gamma_s p v) of a single-phase state (gamma_s = -(d ln p / d ln v)_s, the quantity of RP-1311 eq. 2.71), from the
    gamma_s of the provider's state at (T, p) (one more TP call).  The provider must give gamma_s above 1 (the gate)."""
    if state.status != "ok" or state.region != "single_phase":
        raise StateRecoveryError("sound_speed_single_phase needs a single-phase state with status ok")
    gs = exp.provider.tp(state.T_K, state.p_Pa, exp.B_e).gamma_s
    if gs is None or not (math.isfinite(gs) and gs > 1.0):
        raise StateRecoveryError(f"the provider gave gamma_s = {gs!r} at T = {state.T_K!r} K, p = {state.p_Pa!r} Pa (a value above 1 is required)")
    return math.sqrt(gs * state.p_Pa * state.v_m3_per_kg)


def sound_speed_plateau_sides(exp: Expansion, state: FlowState, *, dlnp: float) -> Tuple[float, float]:
    """(a_lower, a_upper) = sqrt(gamma_{s,T} p v) of a plateau state, with gamma_{s,T} of the low- and high-T end of the plateau
    (plateau_gamma_sT; both are returned, their agreement is not judged).  dlnp is the step of the difference (required)."""
    if state.status != "ok" or state.region != "plateau":
        raise StateRecoveryError("sound_speed_plateau_sides needs a plateau state with status ok")
    g = plateau_gamma_sT(exp.provider, exp.B_e, exp.T_m, state.p_Pa, d_relative=exp.d_relative, dlnp=dlnp)
    return (math.sqrt(g.lower.gamma_sT * state.p_Pa * state.v_m3_per_kg),
            math.sqrt(g.upper.gamma_sT * state.p_Pa * state.v_m3_per_kg))
