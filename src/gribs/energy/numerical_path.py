"""The numerical nozzle flow of an energy run (A-6 stage 2, U-H25-5 group H-I; reports 283, 284; decisions 109, 113, 115, 116, 117).

nozzle_state_from_chamber (the closed form) holds a gamma_s that is constant along the isentrope.  A run that crosses a phase-transition
plateau, or starts on one, does not satisfy that (H-25).  For such evaluations this module replaces the closed form by the numerical
isentrope of gribs.energy.isentrope / isentrope_flow: the throat as the first maximum of the mass flux G, the exit state and the thrust.
Every other evaluation, and every run without a declared transition, uses the closed form unchanged (same values, same TP calls).

The rule that selects the numerical flow (the gate; K-H25-55 (B), K-H25-63 (g2)) is evaluated at every call and is deterministic:
  - no transition declared, or no flow (p_stag <= p_a): the closed form;
  - the chamber is on the plateau (RecoveredState.on_plateau): numerical;
  - T_stag <= T_m: the closed form (an isentropic expansion cools, so the path stays below T_m);
  - T_stag > T_m: x = ln(p_e,closed / p_m,closed) = gamma_s / (gamma_s - 1) ln(T_e,closed / T_m) from the closed form (p_m,closed is the
    pressure at which the closed-form isentrope reaches T_m, the TP-free proxy of p_high).  x > band: the closed form; x < -band:
    numerical; |x| <= band: the path is built (Expansion, numerical exit) and the numerical flow is used only if the numerical exit
    pressure is at or below p_high (the path enters the plateau; a throat that is the kink at p_high does).
At the switch the output has a jump (the difference of the two flows at that point; decision 117 K-H25-66 (J1)); no blending is done.

The numerical flow uses a continuation from the previous evaluation (K-H25-57 (2)): the throat is searched in a bracket of +-bracket_width
(ln p) around the previous throat (scaled with p_0), the exit pressure likewise, and the cold search (find_throat, exit_state, n_grid
points) is used when there is no previous evaluation, no bracket is found, the maximum is not interior, or the classification of the
throat changes.  The cache is the last chamber state's Expansion and the last result (a call with the same chamber state and area).
There is no cache of states between evaluations.

h_0 and s_0 are those of the provider's TP state at the chamber (T, p) (h_0 = u + n R_u T, V0), or, for a chamber on the plateau, the
phase-fraction mixture of the two ends at p with the phase fraction of the recovery; they differ from E_cv / M + p V_free / M by the
root tolerance of the recovery only.  The provider must give the entropy (TPState.entropy_J_per_kgK) whenever the gate needs it; the
gate is evaluated without TP calls, except in the band.
The mass flow is the ideal G* A_t (or G(p_a) A_e when not choked) times Cd; the thrust is eta_thrust mdot u_e + (p_e - p_a) A_e; the
phase split is that of the chamber (nozzle_state_from_chamber).  The Mach number is not computed (None): the speed of sound on a plateau
is not single-valued (D2a/D2b, report 271).  Flow separation is NOT computed: the thrust is the one at A_e, regime stays "choked", and
W_ENERGY_SEPARATION_NOT_MODELLED is recorded when use_separation is true and the numerical exit pressure is below sep_ratio p_a.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

import gribs.nozzle
from gribs.energy.exceptions import ExitOutOfRangeError, StateRecoveryError
from gribs.energy.isentrope import default_d_relative
from gribs.energy.isentrope_flow import Expansion, ExitState, Throat, exit_state, exit_state_near, find_throat, find_throat_near
from gribs.energy.nozzle_from_state import nozzle_state_from_chamber

#: Universal gas constant [J/(kmol K)] (the value of gribs.constants.R_UNIVERSAL; recovery.py and isentrope.py have the same constant).
_RU = 8314.51
SEPARATION_WARNING_CODE = "W_ENERGY_SEPARATION_NOT_MODELLED"
KINK_WARNING_CODE = "W_ENERGY_THROAT_AT_PLATEAU_ENTRY"


@dataclass(frozen=True)
class NumericalPathSettings:
    """The numerical settings (decision 117 item 2; the values of r547, r549 and r550).  plateau_dlnp is the step of the difference
    of gamma_{s,T} (sound_speed_plateau_sides); the run does not use it (no Mach number is computed), the diagnostics do."""
    n_grid: int = 400                   # points of the throat search and of the exit scans (cold)
    p_floor_ratio: float = 1.0e-3       # p_floor = p_0 * p_floor_ratio
    ln_p_xtol: float = 1.0e-10
    kink_dlnp: float = 1.0e-4
    plateau_dlnp: float = 0.025
    gate_band: float = 0.1              # |ln(p_e,closed / p_m,closed)| of the exact judgement
    bracket_width: float = 0.05         # ln p, each side, of the continuation
    T_ladder_K: Tuple[float, ...] = (500.0, 800.0, 1000.0, 1200.0, 1500.0, 1800.0)
    T_upper_K: float = 4500.0           # the upper end of T_bracket (raised to 1.01 T_stag when that is above it)
    continuation: bool = True


DEFAULT_SETTINGS = NumericalPathSettings()


@dataclass
class _Record:
    """The decision and the numbers of one call (also the cache entry)."""
    numeric: bool
    nz: Optional[Dict[str, Any]] = None
    margin: Optional[float] = None
    throat_class: Optional[str] = None
    exit_class: Optional[str] = None
    p_t: Optional[float] = None
    p_e: Optional[float] = None
    separation_ratio: Optional[float] = None       # p_e / p_a when the separation warning condition holds
    gate_x: Optional[float] = None
    gate_reason: str = ""


class NumericalPath:
    """The dispatcher of one run (state: the counters, the last evaluation).  nozzle() returns (nz, margin): nz has the keys of
    nozzle_state_from_chamber; margin is None for the closed form (the caller computes p - p_thr) and, for the numerical flow, a
    number that is positive when the nozzle is choked: eps - G* / G(p_a) for p_a > p_t, (eps - 1) + ln(p_t / p_a) for p_a <= p_t
    (continuous at p_a = p_t; the two have the same sign as p_sub - p_a)."""

    def __init__(self, c, inputs, settings: NumericalPathSettings = DEFAULT_SETTINGS):
        trs = c.energy.condensed.transitions
        self.provider = inputs.provider
        self.declared = bool(trs)
        self.T_m = trs[0].melting_temperature_K if trs else None
        self.d_relative = default_d_relative(c.energy.recovery.root_tolerance)
        self.s = settings
        self.n_numeric = self.n_closed = self.n_switches = self.n_kink = 0
        self._last_numeric: Optional[bool] = None
        self.separation: Dict[str, Any] = {"count": 0, "first_time_s": None, "min_pe_over_pa": None}
        self.kink: Dict[str, Any] = {"count": 0, "first_time_s": None}
        #: counters of the search (not provenance): how often the continuation was used and how often the cold search
        self.stats = {"throat_warm": 0, "throat_cold": 0, "exit_warm": 0, "exit_cold": 0, "exp_built": 0, "gate_exact": 0}
        self.last: Optional[_Record] = None                 # the record of the last call (tools and tests read it)
        self._exp_key = None
        self._exp: Optional[Expansion] = None
        self._res_key = None
        self._res: Optional[_Record] = None
        self._hint: Optional[Dict[str, Any]] = None        # throat class and the ratios p_t / p_0, p_e / p_0 of the last numerical call

    # ------------------------------------------------------------------------------------------------------ the call
    def nozzle(self, *, rec, p_stag: float, T_stag: float, R_eff: float, gamma_s: float, Yg: float, At: float, Ae: float, c, B_e,
               key, t: Optional[float] = None) -> Tuple[Dict[str, Any], Optional[float]]:
        closed = nozzle_state_from_chamber(p_stag, T_stag, R_eff, gamma_s, Yg, At, Ae, c)
        if not self.declared or closed["regime"] == "no-flow":
            self._count(False, None, None)
            return closed, None
        rkey = (key, At)
        if rkey != self._res_key:
            self._res, self._res_key = self._decide(rec, p_stag, T_stag, gamma_s, Yg, At, Ae, c, B_e, key, closed), rkey
        r = self._res
        self.last = r
        self._count(r.numeric, r, t)
        return (r.nz, r.margin) if r.numeric else (closed, None)

    def _count(self, numeric: bool, r: Optional[_Record], t: Optional[float]) -> None:
        if numeric:
            self.n_numeric += 1
            if r.throat_class == "kink_at_plateau_entry":
                if self.kink["count"] == 0:
                    self.kink["first_time_s"] = t
                self.kink["count"] += 1
                self.n_kink += 1
            if r.separation_ratio is not None:
                if self.separation["count"] == 0:
                    self.separation["first_time_s"] = t
                self.separation["count"] += 1
                m = self.separation["min_pe_over_pa"]
                self.separation["min_pe_over_pa"] = r.separation_ratio if m is None else min(m, r.separation_ratio)
        else:
            self.n_closed += 1
        if self._last_numeric is not None and self._last_numeric != numeric:
            self.n_switches += 1
        self._last_numeric = numeric

    # ------------------------------------------------------------------------------------------------------ the gate
    def _decide(self, rec, p_stag, T_stag, gamma_s, Yg, At, Ae, c, B_e, key, closed) -> _Record:
        Tm, band = self.T_m, self.s.gate_band
        if rec.on_plateau:
            return self._numerical(rec, p_stag, T_stag, Yg, At, Ae, c, B_e, key, None, "on_plateau")
        if not T_stag > Tm:
            return _Record(False, gate_reason="T_stag <= T_m")
        pe = self._closed_form_exit_pressure(closed, p_stag, gamma_s, Ae / At, c)
        if not (0.0 < pe < p_stag):
            return _Record(False, gate_reason="no exit pressure below p_stag")
        Te = T_stag * (pe / p_stag) ** ((gamma_s - 1.0) / gamma_s)
        x = gamma_s / (gamma_s - 1.0) * math.log(Te / Tm)
        if x > band:
            return _Record(False, gate_x=x, gate_reason="the closed-form exit is above the band")
        if x < -band:
            return self._numerical(rec, p_stag, T_stag, Yg, At, Ae, c, B_e, key, x, "the closed-form exit is below the band")
        self.stats["gate_exact"] += 1
        return self._numerical(rec, p_stag, T_stag, Yg, At, Ae, c, B_e, key, x, "band", exact=True)

    @staticmethod
    def _closed_form_exit_pressure(closed, p_stag: float, gamma_s: float, eps: float, c) -> float:
        """The exit pressure of the isentropic closed form at the area ratio (p_a when the nozzle is not choked).  closed["pe"] is not used
        for a choked nozzle: with use_separation the closed form returns the pressure at the separation point instead."""
        if closed["regime"] == "subsonic":
            return c.p_a
        Me = 1.0 if eps <= 1.0 + 1.0e-12 else gribs.nozzle.mach_from_area(eps, gamma_s, True)
        return p_stag * (1.0 + 0.5 * (gamma_s - 1.0) * Me * Me) ** (-gamma_s / (gamma_s - 1.0))

    # --------------------------------------------------------------------------------------------- the numerical flow
    def _expansion(self, rec, p_stag, T_stag, B_e, key) -> Expansion:
        if self._exp_key == key and self._exp is not None:
            return self._exp
        Tm, d = self.T_m, self.d_relative
        if rec.on_plateau:
            delta = d * Tm
            lo, hi = self._tp(Tm - delta, p_stag, B_e), self._tp(Tm + delta, p_stag, B_e)
            phi = rec.phi
            s_0 = (1.0 - phi) * self._entropy(lo, Tm - delta, p_stag) + phi * self._entropy(hi, Tm + delta, p_stag)
            h_0 = (1.0 - phi) * (lo.u_J_per_kg + lo.n_gas_per_kg * _RU * (Tm - delta)) \
                + phi * (hi.u_J_per_kg + hi.n_gas_per_kg * _RU * (Tm + delta))
        else:
            st = self._tp(T_stag, p_stag, B_e)
            s_0 = self._entropy(st, T_stag, p_stag)
            h_0 = st.u_J_per_kg + st.n_gas_per_kg * _RU * T_stag
        p_floor = p_stag * self.s.p_floor_ratio
        T_lo = self._T_lo(s_0, p_stag, p_floor, B_e)
        exp = Expansion(self.provider, B_e, s_0, h_0, p_stag, p_floor, (T_lo, max(self.s.T_upper_K, 1.01 * T_stag)), T_m=Tm, d_relative=d)
        self.stats["exp_built"] += 1
        self._exp_key, self._exp = key, exp
        return exp

    def _tp(self, T, p, B_e):
        st = self.provider.tp(T, p, B_e)
        if not st.converged:
            raise StateRecoveryError(f"the equilibrium solve did not converge at T = {T:.10g} K, p = {p:.10g} Pa")
        return st

    @staticmethod
    def _entropy(st, T, p) -> float:
        s = st.entropy_J_per_kgK
        if s is None or not math.isfinite(s):
            raise StateRecoveryError(f"the numerical nozzle flow needs the entropy of the provider's state; none at T = {T:.10g} K, "
                                     f"p = {p:.10g} Pa (got {s!r})")
        return s

    def _T_lo(self, s_0: float, p_0: float, p_floor: float, B_e) -> float:
        """The smallest T of the ladder at which the provider converges at (T, p_0) and (T, p_floor) and s(T, p_floor) < s_0."""
        for cand in self.s.T_ladder_K:
            try:
                a, b = self.provider.tp(cand, p_0, B_e), self.provider.tp(cand, p_floor, B_e)
            except StateRecoveryError:
                continue
            if a.converged and b.converged and b.entropy_J_per_kgK is not None and b.entropy_J_per_kgK < s_0:
                return cand
        raise StateRecoveryError(f"no lower end of T_bracket in {self.s.T_ladder_K} K satisfies s(T, p_floor) < s_0 with a converged "
                                 f"equilibrium at both pressures (s_0 = {s_0!r} J/(kg K))")

    def _path(self, exp: Expansion, p_stag: float, At: float, Ae: float, c) -> Tuple[Throat, ExitState]:
        s = self.s
        eps = Ae / At
        h = self._hint
        th = None
        if s.continuation and h is not None:
            th = find_throat_near(exp, p_hint=p_stag * h["r_t"], width=s.bracket_width, ln_p_xtol=s.ln_p_xtol, kink_dlnp=s.kink_dlnp)
            if th is not None and th.classification != h["throat_class"]:
                th = None
        if th is None:
            th = find_throat(exp, n_grid=s.n_grid, kink_dlnp=s.kink_dlnp, ln_p_xtol=s.ln_p_xtol)
            self.stats["throat_cold"] += 1
        else:
            self.stats["throat_warm"] += 1
        ex = None
        if s.continuation and h is not None and h.get("r_e") is not None and th.p_t_Pa is not None:
            ex = exit_state_near(exp, th, eps=eps, p_a=c.p_a, eta=c.eta_thrust, A_t=At, p_exit_hint=p_stag * h["r_e"], width=s.bracket_width)
            if ex is not None:
                self.stats["exit_warm"] += 1
        if ex is None:
            ex = exit_state(exp, th, eps=eps, p_a=c.p_a, eta=c.eta_thrust, A_t=At, n_grid=s.n_grid)
            self.stats["exit_cold"] += 1
        return th, ex

    def _numerical(self, rec, p_stag, T_stag, Yg, At, Ae, c, B_e, key, x, reason, exact: bool = False) -> _Record:
        exp = self._expansion(rec, p_stag, T_stag, B_e, key)
        th, ex = self._path(exp, p_stag, At, Ae, c)
        if ex.classification == "out_of_range":
            raise ExitOutOfRangeError(
                f"{ExitOutOfRangeError.code}: the isentropic flow of the chamber state (p = {p_stag:.9g} Pa, T = {T_stag:.9g} K) has no "
                f"continuous solution at the area ratio {Ae / At:.9g}: G stops decreasing at p = {ex.fold_p_Pa!r} Pa, before G* / G reaches "
                f"it (the fold of the supersonic branch at the exit of the plateau)")
        if ex.classification not in ("choked_supersonic", "choked_sonic", "unchoked_subsonic", "no_flow") or th.classification in (
                "not_found", "state_error"):
            raise StateRecoveryError(f"the numerical nozzle flow failed (throat: {th.classification} {th.message}; exit: "
                                     f"{ex.classification} {ex.message})")
        if exact and not (exp.plateau == "present" and ex.p_e_Pa <= exp.p_high):       # <=: a throat that is the kink at p_high enters it
            return _Record(False, gate_x=x, gate_reason="the numerical path does not enter the plateau")
        # the hint of the next call, as ratios to p_0
        self._hint = {"throat_class": th.classification, "r_t": th.p_t_Pa / p_stag,
                      "r_e": (ex.p_e_Pa / p_stag) if ex.classification == "choked_supersonic" else None}
        eps = Ae / At
        pe, ue = ex.p_e_Pa, ex.u_e_m_per_s
        split = gribs.nozzle.NozzleTwoPhaseModel(c.two_phase_mode)
        mdot_ref = c.Cd * ex.mdot_kg_per_s
        mdot_total, mdot_gas, mdot_cond = split.split(mdot_ref, Yg)
        F = c.eta_thrust * mdot_total * ue + (pe - c.p_a) * Ae
        nz = dict(mdot=mdot_total, mdot_total=mdot_total, mdot_gas=mdot_gas, mdot_condensed=mdot_cond, mdot_total_reference=mdot_ref,
                  gas_mass_fraction=Yg, condensed_mass_fraction=1.0 - Yg, F=F, Me=None, pe=pe, ve=ue,
                  regime="choked" if ex.choked else "subsonic", A_eff=Ae)
        if not ex.choked:
            margin = eps - th.G_star / ex.state.G if ex.state is not None else (p_stag - c.p_a) / c.p_a
        elif c.p_a <= th.p_t_Pa:
            margin = (eps - 1.0) + math.log(th.p_t_Pa / c.p_a)
        else:
            margin = eps - th.G_star / exp.state_at(c.p_a).G
        sep = pe / c.p_a if (ex.choked and c.use_separation and eps > 1.0 and pe < c.sep_ratio * c.p_a) else None
        return _Record(True, nz, margin, th.classification, ex.classification, th.p_t_Pa, pe, sep, x, reason)

    # ------------------------------------------------------------------------------------------------- the run's record
    def warnings(self) -> List[Dict[str, Any]]:
        """The warnings of the run (one entry per code; empty when none)."""
        out = []
        if self.kink["count"]:
            out.append({
                "code": KINK_WARNING_CODE, "severity": "warning",
                "message": ("the throat of the numerical isentropic flow is at the entry of the phase-transition plateau (the maximum of "
                            "G is a kink at p_high). The state of the flow there rests on an inequality that is derived only "
                            "(a_{s,T} <= u <= a_I at p_high; D2a/D2b of report 271, not RP-1311) and has no independent confirmation; "
                            "the difference from the CEA rocket problem is recorded in report 282; the size of the effect on the run "
                            "has not been evaluated"),
                "context": {"evaluations": self.kink["count"], "first_time_s": self.kink["first_time_s"]}})
        if self.separation["count"]:
            out.append({
                "code": SEPARATION_WARNING_CODE, "severity": "warning",
                "message": ("flow separation is not computed by the numerical nozzle flow: the thrust is the one at the exit area A_e "
                            "and the regime stays 'choked', although use_separation is true and the numerical exit pressure was below "
                            "sep_ratio * p_a"),
                "context": {"evaluations": self.separation["count"], "first_time_s": self.separation["first_time_s"],
                            "min_pe_over_pa": self.separation["min_pe_over_pa"]}})
        return out

    def provenance(self) -> Dict[str, Any]:
        return {"numerical_path_evaluations": self.n_numeric, "closed_form_evaluations": self.n_closed,
                "numerical_path_switches": self.n_switches, "throat_kink_evaluations": self.n_kink}
