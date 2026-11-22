"""A-6 stage 2, U-H25-5 group A-D (stages B, C, D; reports 272, 273; PRE_REGISTRATION_R545): the isentrope building blocks of
gribs.energy.isentrope with SYNTHETIC providers (ideal gas, constant cp, an optional latent jump at T_m; every number is an
ASSUMED INPUT, not data).  The providers are classes defined here (decision 110, K-H25-46 (b)): ReferenceConstantProvider is not
used and gets no entropy.  No CEA, no monkeypatch, no assignment.

Tolerances (decision 110, item 5): T1 a relative difference of at most 1e-8 against a closed form (algebra and roots), at most 1e-6
for a difference quotient (B-1, C-8); the roots are found to ROOT_RTOL <= 1e-12.  T2 the observed order of convergence of a
refinement sequence of at least 3 refinements is at least the theoretical order - 0.3 in the part above the rounding floor.  A
check that does not pass is not loosened.  The theoretical orders (registered before the measurement): central difference in
ln p 2, trapezoid 2, Simpson 4, the limit kappa -> 0 of the ideal gas 1.
"""
import dataclasses
import math

import pytest

from gribs.energy import recovery
from gribs.energy.config_model import CondensedData, CondensedPhaseSpec, PhaseTransitionSpec, RecoverySettings
from gribs.energy.exceptions import StateRecoveryError
from gribs.energy.isentrope import (ROOT_RTOL, default_d_relative, plateau_gamma_sT, plateau_point, plateau_pressures,
                                    single_phase_point)
from gribs.energy.provider import TPState

RU = 8314.51
PREF, TREF, HREF, SREF = 1.0e6, 300.0, -2.0e6, 5000.0
T_M, LATENT = 2000.0, 1.0e6
B_E = {"X": 1.0}
D_REL = 1.0e-9
T1_TOL, T1_DIFF_TOL, T2_MARGIN, FLOOR = 1.0e-8, 1.0e-6, 0.3, 1.0e-10


def rel(a, b):
    return abs(a - b) / abs(b)


@dataclasses.dataclass(frozen=True)
class Gas:
    """An ideal gas of molar mass M and ratio gamma, constant cp; with T_m not None the phase above T_m has the extra enthalpy
    LATENT (and the extra entropy LATENT / T_m).  n(p) = n0 exp(alpha x + c3 x^3), x = ln(p / PREF) (alpha = c3 = 0: n constant).
    side_factor multiplies n above T_m (a gas that differs between the two ends).  report_gamma_s is the gamma_s the state carries
    (None, or a value that the code under test must not read)."""
    M: float = 24.0
    gamma: float = 1.4
    T_m: float = None
    alpha: float = 0.0
    c3: float = 0.0
    side_factor: float = 1.0
    report_gamma_s: float = 1.4
    with_entropy: bool = True

    def n(self, p, T=None):
        x = math.log(p / PREF)
        f = self.side_factor if (self.T_m is not None and T is not None and T > self.T_m) else 1.0
        return f * (1.0 / self.M) * math.exp(self.alpha * x + self.c3 * x ** 3)

    def cp(self):
        return self.gamma / (self.gamma - 1.0) * RU / self.M

    def tp(self, T, p, B_e):
        n = self.n(p, T)
        R = RU * n
        jump = self.T_m is not None and T > self.T_m
        h = HREF + self.cp() * (T - TREF) + (LATENT if jump else 0.0)
        s = SREF + self.cp() * math.log(T / TREF) - R * math.log(p / PREF) + (LATENT / self.T_m if jump else 0.0)
        return TPState(u_J_per_kg=h - n * RU * T, n_gas_per_kg=n, Yg=1.0, converged=True, gamma_s=self.report_gamma_s,
                       entropy_J_per_kgK=s if self.with_entropy else None)


class Counting:
    def __init__(self, inner):
        self.inner, self.calls = inner, []

    def tp(self, T, p, B_e):
        self.calls.append((T, p))
        return self.inner.tp(T, p, B_e)


def orders(errors, theory, floor=FLOOR):
    """Observed orders log2(e_k / e_{k+1}) of a sequence whose step halves, for the pairs whose finer error is above the floor."""
    pairs = [(a, b) for a, b in zip(errors, errors[1:]) if b > floor]
    assert len(pairs) >= 3, errors
    for a, b in pairs:
        assert math.log2(a / b) >= theory - T2_MARGIN, (errors, theory)


@pytest.mark.fast
def test_a01_the_entropy_of_a_tp_state_is_optional_and_defaults_to_none():
    st = TPState(u_J_per_kg=1.0, n_gas_per_kg=0.04, Yg=1.0)
    assert st.entropy_J_per_kgK is None and st.gamma_s is None
    assert TPState(u_J_per_kg=1.0, n_gas_per_kg=0.04, Yg=1.0, entropy_J_per_kgK=7.5).entropy_J_per_kgK == 7.5


# --------------------------------------------------------------------------------------------------- stage B
@pytest.mark.fast
def test_b01_gamma_sT_of_a_power_law_n_is_one_over_one_minus_alpha_on_both_sides():
    g = plateau_gamma_sT(Gas(T_m=T_M, alpha=0.3), B_E, T_M, 2.0e6, d_relative=D_REL, dlnp=0.1)
    for side in (g.lower, g.upper):
        assert rel(side.gamma_sT, 1.0 / (1.0 - 0.3)) <= T1_DIFF_TOL
        assert rel(side.dlnn_dlnp, 0.3) <= T1_DIFF_TOL
    assert g.lower.T_K == T_M - D_REL * T_M and g.upper.T_K == T_M + D_REL * T_M


@pytest.mark.fast
def test_b02_gamma_sT_of_a_constant_n_is_one():
    g = plateau_gamma_sT(Gas(T_m=T_M), B_E, T_M, 2.0e6, d_relative=D_REL, dlnp=0.1)
    assert rel(g.lower.gamma_sT, 1.0) <= T1_TOL and rel(g.upper.gamma_sT, 1.0) <= T1_TOL


@pytest.mark.fast
def test_b03_gamma_sT_does_not_read_the_gamma_s_of_the_state():
    results = [plateau_gamma_sT(Gas(T_m=T_M, alpha=0.3, report_gamma_s=v), B_E, T_M, 2.0e6, d_relative=D_REL, dlnp=0.1)
               for v in (None, 0.5, 2.0, 1.4)]
    assert all(r == results[0] for r in results)


@pytest.mark.fast
def test_b04_gamma_sT_uses_four_tp_calls_at_the_two_ends_and_two_pressures():
    counting = Counting(Gas(T_m=T_M, alpha=0.3))
    plateau_gamma_sT(counting, B_E, T_M, 2.0e6, d_relative=D_REL, dlnp=0.1)
    assert len(counting.calls) == 4
    assert {T for T, _ in counting.calls} == {T_M - D_REL * T_M, T_M + D_REL * T_M}
    assert {p for _, p in counting.calls} == {2.0e6 * math.exp(-0.1), 2.0e6 * math.exp(0.1)}


@pytest.mark.fast
def test_b05_the_two_sides_are_both_returned_and_their_difference_is_not_judged():
    g = plateau_gamma_sT(Gas(T_m=T_M, alpha=0.3, side_factor=1.1, c3=0.02), B_E, T_M, 2.0e6, d_relative=D_REL, dlnp=0.1)
    assert g.lower.n_plus != g.upper.n_plus and math.isfinite(g.lower.gamma_sT) and math.isfinite(g.upper.gamma_sT)


@pytest.mark.fast
def test_b06_invalid_inputs_raise_and_the_step_has_no_default():
    gas = Gas(T_m=T_M)
    with pytest.raises(TypeError):
        plateau_gamma_sT(gas, B_E, T_M, 2.0e6, d_relative=D_REL)                      # dlnp is required
    with pytest.raises(TypeError):
        plateau_gamma_sT(gas, B_E, T_M, 2.0e6, dlnp=0.1)                             # d_relative is required
    for kw in (dict(p=-1.0), dict(p=0.0), dict(d_relative=1.0), dict(d_relative=0.0), dict(dlnp=0.0), dict(dlnp=-0.1),
               dict(dlnp=math.nan), dict(T_m=-1.0)):
        args = dict(p=2.0e6, d_relative=D_REL, dlnp=0.1, T_m=T_M)
        args.update(kw)
        with pytest.raises(StateRecoveryError):
            plateau_gamma_sT(gas, B_E, args["T_m"], args["p"], d_relative=args["d_relative"], dlnp=args["dlnp"])
    with pytest.raises(StateRecoveryError):                                            # 1 - d ln n / d ln p not positive
        plateau_gamma_sT(Gas(T_m=T_M, alpha=1.5), B_E, T_M, 2.0e6, d_relative=D_REL, dlnp=0.1)


@pytest.mark.fast
def test_b07_the_central_difference_converges_with_order_two_in_the_step():
    alpha, c3, p = 0.3, 0.01, 2.0e6
    x = math.log(p / PREF)
    exact = alpha + 3.0 * c3 * x * x
    steps = [0.4 / 2 ** k for k in range(5)]
    errors = [abs(plateau_gamma_sT(Gas(T_m=T_M, alpha=alpha, c3=c3), B_E, T_M, p, d_relative=D_REL, dlnp=h).lower.dlnn_dlnp
                  - exact) / abs(exact) for h in steps]
    orders(errors, 2.0)


# --------------------------------------------------------------------------------------------------- stage C
GAS2 = Gas(T_m=T_M)
DELTA = D_REL * T_M


def _expected_ends():
    cp, R = GAS2.cp(), RU / GAS2.M
    t_lo, t_hi = T_M - DELTA, T_M + DELTA
    h_lo, h_hi = HREF + cp * (t_lo - TREF), HREF + cp * (t_hi - TREF) + LATENT
    s_lo = lambda p: SREF + cp * math.log(t_lo / TREF) - R * math.log(p / PREF)
    s_hi = lambda p: SREF + cp * math.log(t_hi / TREF) - R * math.log(p / PREF) + LATENT / T_M
    return h_lo, h_hi, s_lo, s_hi, R


@pytest.mark.fast
def test_c01_the_phase_fraction_is_the_algebra_of_the_entropies():
    h_lo, h_hi, s_lo, s_hi, R = _expected_ends()
    s_0 = 0.5 * (s_lo(PREF) + s_hi(PREF))                                              # phi = 0.5 at PREF
    for p in (PREF * math.exp(x) for x in (-0.5, -0.2, 0.0, 0.3, 0.5)):
        pt = plateau_point(p, s_0, h_hi + 5.0e5, B_e=B_E, T_m=T_M, provider=GAS2, d_relative=D_REL)
        assert pt.region == "plateau"
        assert rel(pt.phi, (s_0 - s_lo(p)) / (s_hi(p) - s_lo(p))) <= T1_TOL


@pytest.mark.fast
def test_c02_the_enthalpy_the_volume_and_the_velocity_follow_from_the_phase_fraction():
    h_lo, h_hi, s_lo, s_hi, R = _expected_ends()
    s_0, h_0, p = 0.5 * (s_lo(PREF) + s_hi(PREF)), h_hi + 5.0e5, PREF * math.exp(0.3)
    pt = plateau_point(p, s_0, h_0, B_e=B_E, T_m=T_M, provider=GAS2, d_relative=D_REL)
    phi = (s_0 - s_lo(p)) / (s_hi(p) - s_lo(p))
    h = h_lo + phi * (h_hi - h_lo)
    v = RU / GAS2.M * T_M / p
    assert rel(pt.h_J_per_kg, h) <= T1_TOL and rel(pt.v_m3_per_kg, v) <= T1_TOL
    assert rel(pt.u_m_per_s, math.sqrt(2.0 * (h_0 - h))) <= T1_TOL and rel(pt.G, pt.u_m_per_s / v) <= T1_TOL


@pytest.mark.fast
def test_c04_i7_the_enthalpy_jump_is_t_m_times_the_entropy_jump_for_a_provider_that_satisfies_it():
    h_lo, h_hi, s_lo, s_hi, R = _expected_ends()
    pt = plateau_point(PREF, 0.5 * (s_lo(PREF) + s_hi(PREF)), h_hi + 5.0e5, B_e=B_E, T_m=T_M, provider=GAS2, d_relative=D_REL)
    assert rel(pt.h_high - pt.h_low, T_M * (pt.s_high - pt.s_low)) <= T1_TOL


@pytest.mark.fast
def test_c05_outside_the_plateau_a_region_is_returned_and_nothing_is_raised():
    h_lo, h_hi, s_lo, s_hi, R = _expected_ends()
    below = plateau_point(PREF, s_lo(PREF) - 10.0, h_hi + 5.0e5, B_e=B_E, T_m=T_M, provider=GAS2, d_relative=D_REL)
    above = plateau_point(PREF, s_hi(PREF) + 10.0, h_hi + 5.0e5, B_e=B_E, T_m=T_M, provider=GAS2, d_relative=D_REL)
    assert (below.region, above.region) == ("low_side", "high_side") and below.phi < 0.0 < 1.0 < above.phi
    assert below.u_m_per_s is None and above.G is None


@pytest.mark.fast
def test_c06_a_provider_without_entropy_and_inputs_that_do_not_hold_raise():
    with pytest.raises(StateRecoveryError):
        plateau_point(PREF, 0.0, 0.0, B_e=B_E, T_m=T_M, provider=Gas(T_m=T_M, with_entropy=False), d_relative=D_REL)
    h_lo, h_hi, s_lo, s_hi, R = _expected_ends()
    with pytest.raises(StateRecoveryError):                                            # h above the stagnation enthalpy
        plateau_point(PREF, 0.5 * (s_lo(PREF) + s_hi(PREF)), h_lo, B_e=B_E, T_m=T_M, provider=GAS2, d_relative=D_REL)
    with pytest.raises(StateRecoveryError):
        plateau_point(PREF, math.nan, 0.0, B_e=B_E, T_m=T_M, provider=GAS2, d_relative=D_REL)


@pytest.mark.fast
def test_c07_a_plateau_point_uses_two_tp_calls():
    counting = Counting(GAS2)
    h_lo, h_hi, s_lo, s_hi, R = _expected_ends()
    plateau_point(PREF, 0.5 * (s_lo(PREF) + s_hi(PREF)), h_hi + 5.0e5, B_e=B_E, T_m=T_M, provider=counting, d_relative=D_REL)
    assert len(counting.calls) == 2 and {T for T, _ in counting.calls} == {T_M - DELTA, T_M + DELTA}


@pytest.mark.fast
def test_c08_i6_gamma_from_the_path_equals_gamma_sT_of_the_provider_for_a_power_law_n():
    """Gamma = -dln p / dln v along the plateau path (n(p) of plateau_point at a geometric grid) against plateau_gamma_sT.  Both
    use the n of the same provider: this checks the consistency of the difference quotient, not the provider."""
    alpha = 0.3
    gas = Gas(T_m=T_M, alpha=alpha)
    h_lo, h_hi, s_lo, s_hi, R = _expected_ends()
    s_0 = 0.5 * (s_lo(PREF) + s_hi(PREF))
    ps = [PREF * math.exp(-0.2 + 0.05 * i) for i in range(9)]
    pts = [plateau_point(p, s_0, h_hi + 5.0e5, B_e=B_E, T_m=T_M, provider=gas, d_relative=D_REL) for p in ps]
    assert all(q.region == "plateau" for q in pts)
    gammas = [-(math.log(b.p_Pa) - math.log(a.p_Pa)) / (math.log(b.v_m3_per_kg) - math.log(a.v_m3_per_kg)) for a, b in zip(pts, pts[1:])]
    for g in gammas:
        assert rel(g, 1.0 / (1.0 - alpha)) <= T1_DIFF_TOL
    local = plateau_gamma_sT(gas, B_E, T_M, ps[4], d_relative=D_REL, dlnp=0.05)
    assert rel(gammas[0], local.lower.gamma_sT) <= T1_DIFF_TOL


@pytest.mark.fast
def test_c09_the_plateau_ends_of_default_d_relative_are_those_of_the_recovery():
    @dataclasses.dataclass
    class Recorder:
        T: list

        def tp(self, T, p, B_e):
            self.T.append(T)
            return TPState(u_J_per_kg=1.0e6 + 1000.0 * T, n_gas_per_kg=1.0 / 28.0, Yg=1.0)

    rec = Recorder([])
    tm = 2327.0
    cond = CondensedData("neglected", (CondensedPhaseSpec("A(s)", "solid", 3990.0, 1300.0, -1.6e7, 298.15),
                                       CondensedPhaseSpec("A(L)", "liquid", 3000.0, 1300.0, -1.5e7, 298.15)),
                         (PhaseTransitionSpec("A(s)", "A(L)", tm),))
    settings = RecoverySettings(phase_fraction_tolerance=1.0e-6, root_tolerance=1.0e-10, plateau_edge_tolerance=1.0e-6)
    recovery.plateau_bounds(1.0, {"X": 1.0}, cond.transitions[0], rec, cond, settings)
    d = default_d_relative(settings.root_tolerance)
    assert tm - d * tm in rec.T and tm + d * tm in rec.T
    assert default_d_relative(1.0e-12) == 1.0e-9 and default_d_relative(1.0e-6) == pytest.approx(1.0e-5, rel=1.0e-12)


# --------------------------------------------------------------------------------------------------- stage D
def ideal(gamma=1.4):
    return Gas(gamma=gamma)


def _stagnation(gas, T0=3000.0, p0=5.0e6):
    st = gas.tp(T0, p0, B_E)
    return st.entropy_J_per_kgK, st.u_J_per_kg + st.n_gas_per_kg * RU * T0


@pytest.mark.fast
def test_d01_d02_the_single_phase_root_matches_the_closed_form_of_an_ideal_gas():
    gas, p0, T0 = ideal(), 5.0e6, 3000.0
    s_0, h_0 = _stagnation(gas, T0, p0)
    kappa = (gas.gamma - 1.0) / gas.gamma
    for p in (3.0e6, 1.0e6, 2.0e5, 5.0e4):
        res = single_phase_point(p, s_0, h_0, B_E, gas, (300.0, 4000.0))
        assert res.status == "root"
        T = T0 * (p / p0) ** kappa
        pt = res.point
        assert rel(pt.T_K, T) <= T1_TOL
        assert rel(pt.h_J_per_kg, HREF + gas.cp() * (T - TREF)) <= T1_TOL
        assert rel(pt.v_m3_per_kg, RU / gas.M * T / p) <= T1_TOL
        assert rel(pt.u_m_per_s, math.sqrt(2.0 * gas.cp() * (T0 - T))) <= T1_TOL


@pytest.mark.fast
def test_d03_i2_dh_equals_v_dp_with_the_orders_of_the_trapezoid_and_of_simpson():
    gas, p0, T0 = ideal(), 5.0e6, 3000.0
    s_0, h_0 = _stagnation(gas, T0, p0)
    ratio = 1000.0
    e_trap, e_simp = [], []
    for k in range(5):
        n = 4 * 2 ** k
        xs = [-math.log(ratio) * i / n for i in range(n + 1)]                          # ln(p / p0)
        pts = [single_phase_point(p0 * math.exp(x), s_0, h_0, B_E, gas, (100.0, 4000.0)).point for x in xs]
        f = [q.p_Pa * q.v_m3_per_kg for q in pts]                                      # v dp = v p d ln p
        dx = xs[1] - xs[0]
        dh = pts[-1].h_J_per_kg - pts[0].h_J_per_kg
        trap = dx * (0.5 * f[0] + sum(f[1:-1]) + 0.5 * f[-1])
        simp = dx / 3.0 * (f[0] + f[-1] + 4.0 * sum(f[1:-1:2]) + 2.0 * sum(f[2:-1:2]))
        e_trap.append(abs(trap - dh) / abs(dh))
        e_simp.append(abs(simp - dh) / abs(dh))
        # the integral is negative (p decreases): dx < 0, so trap and simp carry the sign of dh
    orders(e_trap, 2.0)
    orders(e_simp, 4.0)


@pytest.mark.fast
def test_d04_the_plateau_pressures_are_the_roots_of_the_end_entropies():
    gas = Gas(T_m=T_M)
    cp, R = gas.cp(), RU / gas.M
    s_0 = SREF + cp * math.log(2200.0 / TREF) - R * math.log(4.0e6 / PREF) + LATENT / T_M
    res = plateau_pressures(s_0, B_E, T_M, gas, (1.0e4, 1.0e8), d_relative=D_REL)
    t_lo, t_hi = T_M - D_REL * T_M, T_M + D_REL * T_M
    p_hi = PREF * math.exp((SREF + cp * math.log(t_hi / TREF) + LATENT / T_M - s_0) / R)
    p_lo = PREF * math.exp((SREF + cp * math.log(t_lo / TREF) - s_0) / R)
    assert res.p_high.status == "root" and res.p_low.status == "root"
    assert rel(res.p_high.p_Pa, p_hi) <= T1_TOL and rel(res.p_low.p_Pa, p_lo) <= T1_TOL and res.p_low.p_Pa < res.p_high.p_Pa


@pytest.mark.fast
def test_d05_a_bracket_that_does_not_hold_the_root_is_reported_not_raised():
    gas = ideal()
    s_0, h_0 = _stagnation(gas)
    assert single_phase_point(1.0e6, s_0, h_0, B_E, gas, (2900.0, 4000.0)).status == "below_bracket"
    assert single_phase_point(1.0e6, s_0, h_0, B_E, gas, (300.0, 1000.0)).status == "above_bracket"
    g2 = Gas(T_m=T_M)
    s_0 = SREF + g2.cp() * math.log(2200.0 / TREF) - RU / g2.M * math.log(4.0e6 / PREF) + LATENT / T_M
    assert plateau_pressures(s_0, B_E, T_M, g2, (1.0e4, 1.0e5), d_relative=D_REL).p_high.status == "above_bracket"
    assert plateau_pressures(s_0, B_E, T_M, g2, (1.0e9, 1.0e10), d_relative=D_REL).p_low.status == "below_bracket"


@pytest.mark.fast
def test_d06_a_provider_without_entropy_and_bad_inputs_raise():
    gas = Gas(with_entropy=False)
    with pytest.raises(StateRecoveryError):
        single_phase_point(1.0e6, 1.0, 1.0, B_E, gas, (300.0, 4000.0))
    with pytest.raises(StateRecoveryError):
        plateau_pressures(1.0, B_E, T_M, Gas(T_m=T_M, with_entropy=False), (1.0e5, 1.0e6), d_relative=D_REL)
    with pytest.raises(StateRecoveryError):
        single_phase_point(1.0e6, 1.0, 1.0, B_E, ideal(), (500.0, 400.0))
    with pytest.raises(StateRecoveryError):
        plateau_pressures(1.0, B_E, T_M, Gas(T_m=T_M), (2.0e6, 1.0e6), d_relative=D_REL)
    assert ROOT_RTOL <= 1.0e-12


@pytest.mark.fast
def test_d07_l1_the_isothermal_limit_of_the_ideal_gas_is_approached_with_order_one_in_kappa():
    T0, p0, p = 3000.0, 5.0e6, 1.0e6
    errors = []
    for k in range(5):
        kappa = 0.2 / 2 ** k
        gas = ideal(gamma=1.0 / (1.0 - kappa))
        s_0, h_0 = _stagnation(gas, T0, p0)
        pt = single_phase_point(p, s_0, h_0, B_E, gas, (300.0, 4000.0)).point
        limit = RU / gas.M * T0 * math.log(p0 / p)                                    # u^2 / 2 of an isothermal ideal gas
        errors.append(abs(pt.u_m_per_s ** 2 / 2.0 - limit) / limit)
    orders(errors, 1.0)


@pytest.mark.fast
def test_d08_h_above_h_0_within_the_rounding_allowance_gives_u_zero_and_beyond_it_raises():
    """Decision 111 item 3: h - h_0 <= 1e-9 (|h_0| + |h|) -> u = 0; larger -> StateRecoveryError (single phase and plateau)."""
    gas = ideal()
    s_0, h_0 = _stagnation(gas)
    h_at = gas.tp(3000.0, 5.0e6, B_E)
    h = h_at.u_J_per_kg + h_at.n_gas_per_kg * RU * 3000.0
    slack = 1.0e-9 * (abs(h) + abs(h))
    inside = single_phase_point(5.0e6, s_0, h - 0.5 * slack, B_E, gas, (300.0, 4000.0))
    assert inside.status == "root" and inside.point.u_m_per_s == 0.0 and inside.point.G == 0.0
    with pytest.raises(StateRecoveryError):
        single_phase_point(5.0e6, s_0, h - 2.0 * slack, B_E, gas, (300.0, 4000.0))
    h_lo, h_hi, s_lo, s_hi, R = _expected_ends()
    s_mid = 0.5 * (s_lo(PREF) + s_hi(PREF))
    h_plateau = 0.5 * (h_lo + h_hi)                                                    # phi = 0.5 at PREF: h = (h_lo + h_hi) / 2
    slack_p = 1.0e-9 * 2.0 * abs(h_plateau)
    inside_p = plateau_point(PREF, s_mid, h_plateau - 0.5 * slack_p, B_e=B_E, T_m=T_M, provider=GAS2, d_relative=D_REL)
    assert inside_p.region == "plateau" and inside_p.u_m_per_s == 0.0
    with pytest.raises(StateRecoveryError):
        plateau_point(PREF, s_mid, h_plateau - 2.0 * slack_p, B_e=B_E, T_m=T_M, provider=GAS2, d_relative=D_REL)
