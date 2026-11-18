"""State recovery across a phase-transition plateau (A-6 stage 2, C5; contract 164 section 2b-FIX-1).

Given the total energy E_cv of the contents, the free volume V_free and the element inventory B_e, find (p, T) and, on a
plateau, the phase fraction.  Method (the contract leaves it to the implementation; this is "accept by inclusion of E"):
for each T the pressure that fills the volume is found (the condensed volume uses constant densities), then the residual
g(T) = M u(T, p(T)) - E_cv is solved in T from each start point.  A sign change of g with a small residual is an ordinary
root.  A sign change with a large residual is a jump of the energy: if a declared transition has its T_m there, the state
is on the plateau, T = T_m and the phase fraction is the LINEAR APPROXIMATION (E_cv - E-)/(E+ - E-) (E-, E+ = the energy just
below and just above T_m at the same V_free and B_e); it ignores the change of p along the plateau that the melting causes
(of order 1e-4 .. 1e-3 in ln p), so a tolerance on the phase fraction must be larger than that.  A jump where no transition
is declared, or no sign change in the searched range, is NoRootOutsidePlateauError.  The plateau is never rejected as
"no root".  The contents are taken to be at rest (E_cv is the internal energy; the kinetic energy at the nozzle inlet is in
H0, contract 2b-4.4).  There is no direct solve of V and U (UV): the report 181 reasons are in contract 2b-FIX-1.

Start points whose results disagree in (p, T) or in the phase fraction are an error; the module never picks one root.
"""
import dataclasses
import math
from typing import Dict, Mapping, NamedTuple, Optional, Sequence, Tuple

from scipy.optimize import brentq

from gribs.energy.config_model import CondensedData, CondensedPhaseSpec, PhaseTransitionSpec, RecoverySettings
from gribs.energy.exceptions import NoRootOutsidePlateauError, StateRecoveryError
from gribs.energy.provider import EquilibriumProvider, TPState
from gribs.energy.state import RecoveredState

#: Universal gas constant [J/(kmol K)] (the value of gribs.constants.R_UNIVERSAL; a test compares them).
_RU = 8314.51
_MAX_DOUBLINGS = 200
#: Start points agree when T and ln p differ by less than this factor times root_tolerance.
_AGREE_FACTOR = 1.0e4
_VOLUME_RESIDUAL_TOL = 1.0e-8
_LATENT_HEAT_RTOL = 1.0e-6


class PlateauBounds(NamedTuple):
    E_minus: float          # M u just below T_m [J]
    E_plus: float           # M u just above T_m [J]
    p_minus: float
    p_plus: float
    state_minus: TPState
    state_plus: TPState


def latent_heat(transition: PhaseTransitionSpec, phases: Sequence[CondensedPhaseSpec]) -> float:
    """L [J/kg] = h_liquid(T_m) - h_solid(T_m) from the constant-cp data of the two phases (the one source of L)."""
    by_species = {p.species: p for p in phases}
    missing = [s for s in (transition.solid_species, transition.liquid_species) if s not in by_species]
    if missing:
        raise StateRecoveryError(
            f"phase transition {transition.solid_species!r} -> {transition.liquid_species!r}: the data of {missing} is "
            f"missing (T_m and the enthalpy of both phases are required; a transition without them is out of scope)")
    Tm = transition.melting_temperature_K

    def h(phase: CondensedPhaseSpec) -> float:
        return phase.enthalpy_reference_J_kg + phase.cp_J_kgK * (Tm - phase.reference_temperature_K)
    return h(by_species[transition.liquid_species]) - h(by_species[transition.solid_species])


def plateau_phase_fraction(E_cv: float, E_minus: float, E_plus: float) -> float:
    """Phase fraction on a plateau, the linear APPROXIMATION (E_cv - E-)/(E+ - E-) (0 = all solid, 1 = all liquid)."""
    width = E_plus - E_minus
    if not width > 0.0:
        raise StateRecoveryError(f"the plateau has no width: E- = {E_minus!r}, E+ = {E_plus!r}")
    return (E_cv - E_minus) / width


class _Solver:
    def __init__(self, E_cv, V_free, B_e, provider, condensed, recovery):
        self.E, self.V, self.B = E_cv, V_free, B_e
        self.provider, self.cond, self.rec = provider, condensed, recovery
        self.M = math.fsum(B_e.values())
        self.rho = {p.species: p.density_kg_m3 for p in condensed.phases}
        if not (math.isfinite(self.M) and self.M > 0.0):
            raise StateRecoveryError(f"the element inventory must have a positive mass, got {self.M!r} kg")
        if not (math.isfinite(V_free) and V_free > 0.0):
            raise StateRecoveryError(f"V_free must be positive and finite, got {V_free!r}")
        self.scale = max(abs(E_cv), 1.0e3 * self.M)

    def _tp(self, T: float, p: float) -> TPState:
        st = self.provider.tp(T, p, self.B)
        if not st.converged:
            raise StateRecoveryError(f"the equilibrium solve did not converge at T = {T:.10g} K, p = {p:.10g} Pa")
        return st

    def _v_condensed(self, st: TPState) -> float:
        total = 0.0
        for species, y in st.condensed_mass_fraction.items():
            rho = self.rho.get(species)
            if rho is None:
                raise StateRecoveryError(f"no density for the condensed species {species!r} (no common constant rho_c is used)")
            total += y / rho
        return self.M * total

    def _v_gas(self, st: TPState, T: float, p: float) -> float:
        return self.M * st.n_gas_per_kg * _RU * T / p

    def _volume_residual(self, T: float, p: float) -> float:
        st = self._tp(T, p)
        return self._v_gas(st, T, p) + self._v_condensed(st) - self.V

    def pressure(self, T: float) -> Tuple[float, TPState]:
        """The pressure at which the gas and the condensed phases fill V_free at temperature T, and the state there."""
        st0 = self._tp(T, 1.0e5)
        # the root lies at or above the pressure that the gas alone would have (the condensed volume only raises it); the
        # bracket starts at half of it, so that the root is inside the bracket and never on an end (no condensed phase:
        # the root is exactly that pressure)
        lo = 0.5 * self.M * st0.n_gas_per_kg * _RU * T / self.V
        for _ in range(_MAX_DOUBLINGS):
            if self._volume_residual(T, lo) > 0.0:
                break
            lo *= 0.5
        else:
            raise StateRecoveryError(f"no pressure bracket at T = {T:.10g} K")
        hi = 4.0 * lo
        for _ in range(_MAX_DOUBLINGS):
            if self._volume_residual(T, hi) <= 0.0:
                break
            lo, hi = hi, 2.0 * hi
        else:
            raise StateRecoveryError(f"V_free is not larger than the condensed volume at T = {T:.10g} K")
        lp = brentq(lambda x: self._volume_residual(T, math.exp(x)), math.log(lo), math.log(hi),
                    xtol=1.0e-14, rtol=max(self.rec.root_tolerance, 1.0e-14))
        p = math.exp(lp)
        st = self._tp(T, p)
        if abs(self._v_gas(st, T, p) + self._v_condensed(st) - self.V) > _VOLUME_RESIDUAL_TOL * self.V:
            raise StateRecoveryError(f"the volume residual is above {_VOLUME_RESIDUAL_TOL:g} at T = {T:.10g} K")
        return p, st

    def residual(self, T: float) -> float:
        _, st = self.pressure(T)
        return self.M * st.u_J_per_kg - self.E

    def bracket(self, T0: float, T_lo: float, T_hi: float) -> Tuple[float, float]:
        g0 = self.residual(T0)
        if g0 == 0.0:
            return T0, T0
        delta, sign = 0.01 * T0, (1.0 if g0 < 0.0 else -1.0)
        a = T0
        while True:
            b = a + sign * delta
            b = min(b, T_hi) if sign > 0 else max(b, T_lo)
            gb = self.residual(b)
            if (gb >= 0.0) if sign > 0 else (gb <= 0.0):
                return (a, b) if sign > 0 else (b, a)
            if b in (T_lo, T_hi):
                raise NoRootOutsidePlateauError(
                    f"no sign change of the energy residual between {T_lo:g} K and {T_hi:g} K (start {T0:g} K): "
                    f"E_cv = {self.E:.10g} J is out of the searched range")
            a, delta = b, delta * 1.5

    def plateau_bounds(self, transition: PhaseTransitionSpec) -> PlateauBounds:
        Tm = transition.melting_temperature_K
        d = max(1.0e-9, 10.0 * self.rec.root_tolerance) * Tm
        p_minus, s_minus = self.pressure(Tm - d)
        p_plus, s_plus = self.pressure(Tm + d)
        return PlateauBounds(self.M * s_minus.u_J_per_kg, self.M * s_plus.u_J_per_kg, p_minus, p_plus, s_minus, s_plus)

    def _state(self, p: float, T: float, st: TPState, phi: Optional[float], on_plateau: bool) -> RecoveredState:
        v_c = self._v_condensed(st)
        return RecoveredState(p_Pa=p, T_K=T, phi=phi, on_plateau=on_plateau, m_total_kg=self.M,
                              V_gas_m3=self._v_gas(st, T, p), V_condensed_m3=v_c, n_gas_kmol=self.M * st.n_gas_per_kg,
                              Yg=st.Yg, condensed_mass_fraction=dict(st.condensed_mass_fraction))

    def on_plateau(self, transition: PhaseTransitionSpec) -> RecoveredState:
        b = self.plateau_bounds(transition)
        frac = plateau_phase_fraction(self.E, b.E_minus, b.E_plus)
        tol = self.rec.plateau_edge_tolerance
        if frac < -tol or frac > 1.0 + tol:
            raise StateRecoveryError(
                f"E_cv = {self.E:.10g} J is outside the plateau [{b.E_minus:.10g}, {b.E_plus:.10g}] J although the "
                f"energy residual jumps at T_m = {transition.melting_temperature_K:g} K")
        phi = min(max(frac, 0.0), 1.0)
        sm, sp = b.state_minus, b.state_plus
        n = (1.0 - phi) * sm.n_gas_per_kg + phi * sp.n_gas_per_kg
        yg = (1.0 - phi) * sm.Yg + phi * sp.Yg
        cond: Dict[str, float] = {k: (1.0 - phi) * sm.condensed_mass_fraction.get(k, 0.0)
                                  + phi * sp.condensed_mass_fraction.get(k, 0.0)
                                  for k in set(sm.condensed_mass_fraction) | set(sp.condensed_mass_fraction)}
        interpolated = TPState(u_J_per_kg=self.E / self.M, n_gas_per_kg=n, Yg=yg, condensed_mass_fraction=cond)
        v_c = self._v_condensed(interpolated)
        if not self.V - v_c > 0.0:
            raise StateRecoveryError("V_free is not larger than the condensed volume on the plateau")
        Tm = transition.melting_temperature_K
        p = self.M * n * _RU * Tm / (self.V - v_c)
        return self._state(p, Tm, interpolated, phi, True)

    def solve_from(self, T0: float, T_lo: float, T_hi: float) -> RecoveredState:
        a, b = self.bracket(T0, T_lo, T_hi)
        T = a if a == b else brentq(self.residual, a, b, xtol=1.0e-12, rtol=max(self.rec.root_tolerance, 1.0e-14))
        jump = abs(self.residual(T)) > max(1.0e3 * self.rec.root_tolerance, 1.0e-9) * self.scale
        if not jump:
            p, st = self.pressure(T)
            return self._state(p, T, st, None, False)
        for tr in self.cond.transitions:
            if abs(tr.melting_temperature_K - T) <= 1.0e-6 * tr.melting_temperature_K:
                return self.on_plateau(tr)
        raise NoRootOutsidePlateauError(
            f"the energy residual jumps at T = {T:.10g} K, where no phase transition is declared (a transition needs "
            f"T_m and the data of both phases in energy_model.condensed)")


def plateau_bounds(V_free: float, B_e: Mapping[str, float], transition: PhaseTransitionSpec,
                   provider: EquilibriumProvider, condensed: CondensedData, recovery: RecoverySettings) -> PlateauBounds:
    """E-, E+ and the pressures just below and just above T_m, at the same V_free and B_e."""
    latent_heat(transition, condensed.phases)
    return _Solver(0.0, V_free, B_e, provider, condensed, recovery).plateau_bounds(transition)


def recover_state(E_cv: float, V_free: float, B_e: Mapping[str, float], provider: EquilibriumProvider,
                  condensed: CondensedData, recovery: RecoverySettings, *, T_search_K: Tuple[float, float],
                  starts_K: Sequence[float], latent_heat_J_per_kg: Optional[float] = None) -> RecoveredState:
    """Recover (p, T, phase fraction) from (E_cv, V_free, B_e); see the module docstring.

    T_search_K = (T_lo, T_hi) is the searched temperature range, starts_K the temperatures the search starts from (every
    start must give the same state).  latent_heat_J_per_kg, if given, is only checked against the L the two phases' data
    give (L is never defined twice)."""
    for tr in condensed.transitions:
        L = latent_heat(tr, condensed.phases)                       # raises when a phase's data is missing
        if latent_heat_J_per_kg is not None and \
                abs(L - latent_heat_J_per_kg) > _LATENT_HEAT_RTOL * max(abs(L), abs(latent_heat_J_per_kg), 1.0):
            raise StateRecoveryError(
                f"latent heat {latent_heat_J_per_kg!r} J/kg does not match the {L!r} J/kg that the data of the two phases "
                f"of {tr.solid_species!r} -> {tr.liquid_species!r} give")
    T_lo, T_hi = T_search_K
    if not (0.0 < T_lo < T_hi):
        raise StateRecoveryError(f"T_search_K must satisfy 0 < T_lo < T_hi, got {T_search_K!r}")
    if not starts_K or any(not (T_lo < t < T_hi) for t in starts_K):
        raise StateRecoveryError(f"every start must lie inside the searched range ({T_lo:g}, {T_hi:g}) K, got {list(starts_K)!r}")
    solver = _Solver(E_cv, V_free, B_e, provider, condensed, recovery)
    results = [solver.solve_from(t, T_lo, T_hi) for t in starts_K]
    ref, agree = results[0], _AGREE_FACTOR * max(recovery.root_tolerance, 1.0e-15)
    for other in results[1:]:
        dT = abs(other.T_K - ref.T_K) / ref.T_K
        dlnp = abs(math.log(other.p_Pa / ref.p_Pa))
        dphi = abs(other.phi - ref.phi) if (ref.phi is not None and other.phi is not None) else 0.0
        if dT > agree or dlnp > agree or dphi > recovery.phase_fraction_tolerance:
            raise StateRecoveryError(
                f"the start points disagree: (T, p, phi) = ({ref.T_K:.10g} K, {ref.p_Pa:.10g} Pa, {ref.phi!r}) and "
                f"({other.T_K:.10g} K, {other.p_Pa:.10g} Pa, {other.phi!r}); no root is chosen")
    return dataclasses.replace(ref, n_starts_agreed=len(results))
