"""A-6 stage 2, U-H25-5 group E-G (stages E, F, G; reports 277, 278; PRE_REGISTRATION_R547): the throat, the exit state, the thrust
and the fold of the supersonic branch of gribs.energy.isentrope_flow with SYNTHETIC providers (the ideal gas with an optional latent
jump at T_m of tests/test_energy_isentrope.py; every number is an ASSUMED INPUT, not data).  No CEA, no monkeypatch, no assignment.

Tolerances (decision 113 item 6 = decision 110 item 5): T1 a relative difference of at most 1e-8 against a closed form (algebra and
roots; the throat is a polished maximum, whose position is found to well below 1e-8); T2 the observed order of a refinement sequence
(at least 3 refinements, above the rounding floor) is at least the theoretical order - 0.3: the limit gamma -> 1 of the ideal gas 1
(the throat pressure and the flux G*), the trapezoid 2 and Simpson 4 for dh = v dp along the plateau.  The inequality a_{s,T} <= u <= a_I
(I4) is a structural check with no tolerance.  A check that does not pass is not loosened.
"""
import math

import pytest
from scipy.optimize import brentq

from gribs.energy.exceptions import StateRecoveryError
from gribs.energy.isentrope import plateau_gamma_sT
from gribs.energy.isentrope_flow import (Expansion, exit_state, find_throat, sound_speed_plateau_sides, sound_speed_single_phase,
                                         supersonic_limit)
from test_energy_isentrope import B_E, D_REL, RU, T1_TOL, T_M, Counting, Gas, _stagnation, orders, rel

N_GRID, KINK, XTOL = 200, 1.0e-3, 1.0e-10
DELTA = D_REL * T_M
M_MOL = 24.0
R_GAS = RU / M_MOL


def ideal_expansion(gamma, T0=3000.0, p0=5.0e6, p_floor=None, **kw):
    gas = Gas(gamma=gamma, report_gamma_s=gamma)
    s_0, h_0 = _stagnation(gas, T0, p0)
    return gas, Expansion(gas, B_E, s_0, h_0, p0, p_floor or 1.0e-4 * p0, (100.0, 5000.0), **kw)


def closed_throat(gamma, T0, p0):
    p_t = p0 * (2.0 / (gamma + 1.0)) ** (gamma / (gamma - 1.0))
    G_star = p0 * math.sqrt(gamma / (R_GAS * T0)) * (2.0 / (gamma + 1.0)) ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))
    return p_t, G_star


def area_ratio(M, g):
    return (1.0 / M) * ((2.0 / (g + 1.0)) * (1.0 + 0.5 * (g - 1.0) * M * M)) ** ((g + 1.0) / (2.0 * (g - 1.0)))


def closed_exit(gamma, T0, p0, eps, supersonic=True):
    g = gamma
    M = brentq(lambda m: area_ratio(m, g) - eps, 1.0 + 1.0e-12 if supersonic else 1.0e-9, 60.0 if supersonic else 1.0 - 1.0e-12,
               xtol=1.0e-15, rtol=1.0e-15)
    f = 1.0 + 0.5 * (g - 1.0) * M * M
    return p0 * f ** (-g / (g - 1.0)), M * math.sqrt(g * R_GAS * T0 / f)


# ------------------------------------------------------------------------------------------------- the plateau chambers
def plateau_chamber(gas, p0, phi0):
    lo, hi = gas.tp(T_M - DELTA, p0, B_E), gas.tp(T_M + DELTA, p0, B_E)
    h_lo, h_hi = lo.u_J_per_kg + lo.n_gas_per_kg * RU * (T_M - DELTA), hi.u_J_per_kg + hi.n_gas_per_kg * RU * (T_M + DELTA)
    return (1.0 - phi0) * lo.entropy_J_per_kgK + phi0 * hi.entropy_J_per_kgK, (1.0 - phi0) * h_lo + phi0 * h_hi


def isothermal_expansion(phi0, p0=5.0e6, p_floor=None, **gas_kw):
    gas = Gas(T_m=T_M, **gas_kw)
    s_0, h_0 = plateau_chamber(gas, p0, phi0)
    return gas, Expansion(gas, B_E, s_0, h_0, p0, p_floor or p0 * math.exp(-1.6), (100.0, 5000.0), T_m=T_M, d_relative=D_REL)


# ---------------------------------------------------------------------------------------------------------- stage E
@pytest.mark.fast
@pytest.mark.parametrize("gamma", [1.2, 1.4, 1.67])
def test_e01_l3_the_throat_of_a_single_phase_ideal_gas_is_the_closed_form(gamma):
    T0, p0 = 3000.0, 5.0e6
    _, exp = ideal_expansion(gamma, T0, p0)
    th = find_throat(exp, n_grid=N_GRID, kink_dlnp=KINK, ln_p_xtol=XTOL)
    p_t, G_star = closed_throat(gamma, T0, p0)
    assert th.classification == "single_phase" and th.state.status == "ok"
    assert rel(th.p_t_Pa, p_t) <= T1_TOL and rel(th.G_star, G_star) <= T1_TOL


@pytest.mark.fast
def test_e02_i3_the_flow_speed_equals_the_sound_speed_at_a_smooth_throat():
    _, exp = ideal_expansion(1.4)
    th = find_throat(exp, n_grid=N_GRID, kink_dlnp=KINK, ln_p_xtol=XTOL)
    a = sound_speed_single_phase(exp, th.state)
    assert abs(th.state.u_m_per_s ** 2 - a ** 2) / th.state.u_m_per_s ** 2 <= T1_TOL


@pytest.mark.fast
def test_e03_i5_the_throat_does_not_depend_on_the_grid():
    _, exp = ideal_expansion(1.4)
    found = [find_throat(exp, n_grid=n, kink_dlnp=KINK, ln_p_xtol=XTOL) for n in (50, 100, 200, 400)]
    for th in found[1:]:
        assert rel(th.p_t_Pa, found[0].p_t_Pa) <= T1_TOL and rel(th.G_star, found[0].G_star) <= T1_TOL


@pytest.mark.fast
def test_e04_l1_the_isothermal_limit_of_the_throat_is_approached_with_order_one_in_kappa():
    T0, p0 = 3000.0, 5.0e6
    e_p, e_G = [], []
    for k in range(5):
        kappa = 0.2 / 2 ** k
        gamma = 1.0 / (1.0 - kappa)
        _, exp = ideal_expansion(gamma, T0, p0)
        th = find_throat(exp, n_grid=N_GRID, kink_dlnp=KINK, ln_p_xtol=XTOL)
        e_p.append(abs(th.p_t_Pa - p0 * math.exp(-0.5)) / (p0 * math.exp(-0.5)))
        e_G.append(abs(th.G_star - p0 * math.exp(-0.5) / math.sqrt(R_GAS * T0)) / (p0 * math.exp(-0.5) / math.sqrt(R_GAS * T0)))
    orders(e_p, 1.0)
    orders(e_G, 1.0)


@pytest.mark.fast
def test_e05_the_search_does_not_read_the_gamma_s_of_the_state():
    results = []
    for g in (None, 0.5, 2.0, 1.4):
        gas = Gas(report_gamma_s=g)
        s_0, h_0 = _stagnation(gas)
        exp = Expansion(gas, B_E, s_0, h_0, 5.0e6, 5.0e2, (100.0, 5000.0))
        results.append(find_throat(exp, n_grid=N_GRID, kink_dlnp=KINK, ln_p_xtol=XTOL))
    assert all(r == results[0] for r in results) and results[0].classification == "single_phase"


@pytest.mark.fast
def test_e06_l2_the_throat_of_the_isothermal_plateau_flow_is_p_0_exp_minus_one_half():
    p0 = 5.0e6
    _, exp = isothermal_expansion(0.4, p0)
    th = find_throat(exp, n_grid=N_GRID, kink_dlnp=KINK, ln_p_xtol=XTOL)
    assert th.classification == "inside_plateau" and th.state.region == "plateau" and 0.0 < th.state.phi < 1.0
    assert rel(th.p_t_Pa, p0 * math.exp(-0.5)) <= T1_TOL
    assert rel(th.G_star, p0 * math.exp(-0.5) / math.sqrt(R_GAS * T_M)) <= T1_TOL
    assert exp.has_plateau and exp.p_low < th.p_t_Pa < exp.p_high


@pytest.mark.fast
def test_e07_l4_a_throat_below_the_plateau_is_single_phase_low():
    _, exp = isothermal_expansion(0.2)
    th = find_throat(exp, n_grid=N_GRID, kink_dlnp=KINK, ln_p_xtol=XTOL)
    assert th.classification == "single_phase_low" and th.p_t_Pa < exp.p_low and th.state.region == "single_phase"
    a = sound_speed_single_phase(exp, th.state)
    assert abs(th.state.u_m_per_s ** 2 - a ** 2) / th.state.u_m_per_s ** 2 <= T1_TOL


@pytest.mark.fast
def test_e08_l4_a_throat_above_the_plateau_is_single_phase_high_and_the_closed_form():
    T0, p0 = 3000.0, 5.0e6
    gas, exp = ideal_expansion(1.4, T0, p0, T_m=T_M, d_relative=D_REL)
    th = find_throat(exp, n_grid=N_GRID, kink_dlnp=KINK, ln_p_xtol=XTOL)
    p_t, G_star = closed_throat(1.4, T0, p0)
    assert exp.has_plateau and exp.p_high < th.p_t_Pa
    assert th.classification == "single_phase_high" and rel(th.p_t_Pa, p_t) <= T1_TOL and rel(th.G_star, G_star) <= T1_TOL


def _kink_expansion(T0=2340.0, p0=5.0e6):
    gas = Gas(T_m=T_M, report_gamma_s=1.4)
    s_0, h_0 = _stagnation(gas, T0, p0)
    return gas, Expansion(gas, B_E, s_0, h_0, p0, 1.0e-3 * p0, (100.0, 5000.0), T_m=T_M, d_relative=D_REL)


@pytest.mark.fast
def test_e09_l4_i4_the_maximum_of_g_is_a_kink_at_the_plateau_entry_between_the_two_sound_speeds():
    T0, p0 = 2340.0, 5.0e6
    gas, exp = _kink_expansion(T0, p0)
    th = find_throat(exp, n_grid=N_GRID, kink_dlnp=KINK, ln_p_xtol=XTOL)
    cp_over_R = gas.gamma / (gas.gamma - 1.0)
    p_high = p0 * ((T_M + DELTA) / T0) ** cp_over_R
    assert th.classification == "kink_at_plateau_entry" and rel(th.p_t_Pa, p_high) <= T1_TOL and rel(exp.p_high, p_high) <= T1_TOL
    u = th.state.u_m_per_s
    a_I = sound_speed_single_phase(exp, th.state)
    a_sT = math.sqrt(plateau_gamma_sT(gas, B_E, T_M, th.p_t_Pa, d_relative=D_REL, dlnp=0.05).lower.gamma_sT
                     * th.p_t_Pa * th.state.v_m3_per_kg)
    assert a_sT <= u <= a_I                                                            # I4: no tolerance
    below = exp.state_at(th.p_t_Pa * math.exp(-KINK))
    above = exp.state_at(th.p_t_Pa * math.exp(KINK))
    assert below.region == "plateau" and above.region == "single_phase" and th.G_star > below.G and th.G_star > above.G


@pytest.mark.fast
def test_e10_a_state_is_cached_and_a_plateau_state_uses_two_tp_calls():
    gas = Gas(T_m=T_M)
    s_0, h_0 = plateau_chamber(gas, 5.0e6, 0.4)
    counting = Counting(gas)
    exp = Expansion(counting, B_E, s_0, h_0, 5.0e6, 5.0e6 * math.exp(-1.6), (100.0, 5000.0), T_m=T_M, d_relative=D_REL)
    n0 = len(counting.calls)
    exp.state_at(4.0e6)
    assert len(counting.calls) == n0 + 2 and {T for T, _ in counting.calls[n0:]} == {T_M - DELTA, T_M + DELTA}
    exp.state_at(4.0e6)
    assert len(counting.calls) == n0 + 2


@pytest.mark.fast
def test_e11_i2_dh_equals_v_dp_on_each_segment_of_a_path_through_the_plateau():
    """The path of a hot chamber: single phase above the plateau, the plateau (n constant, T = T_m), single phase below it.  The
    rules are applied on each smooth segment separately (f = p v has kinks at the plateau ends).  Above the plateau: the orders of
    the trapezoid 2 and of Simpson 4 (T2).  On the plateau f = n R_u T_m is constant: dh = n R_u T_m ln(p_2 / p_1) (T1)."""
    gas = Gas(T_m=T_M, report_gamma_s=1.4)
    p0 = 5.0e6
    s_0, h_0 = _stagnation(gas, 3000.0, p0)
    exp = Expansion(gas, B_E, s_0, h_0, p0, 1.0e3, (100.0, 5000.0), T_m=T_M, d_relative=D_REL)
    assert exp.has_plateau and 1.0e3 < exp.p_low < exp.p_high < p0

    def rules(p_start, p_end, count):
        e_trap, e_simp = [], []
        for k in range(count):
            n = 4 * 2 ** k
            xs = [math.log(p_start) + (math.log(p_end) - math.log(p_start)) * i / n for i in range(n + 1)]
            pts = [exp.state_at(math.exp(x)) for x in xs]
            assert all(q.status == "ok" for q in pts)
            f = [q.p_Pa * q.v_m3_per_kg for q in pts]
            dx = xs[1] - xs[0]
            dh = pts[-1].h_J_per_kg - pts[0].h_J_per_kg
            e_trap.append(abs(dx * (0.5 * f[0] + sum(f[1:-1]) + 0.5 * f[-1]) - dh) / abs(dh))
            e_simp.append(abs(dx / 3.0 * (f[0] + f[-1] + 4.0 * sum(f[1:-1:2]) + 2.0 * sum(f[2:-1:2])) - dh) / abs(dh))
        return e_trap, e_simp

    e_trap, e_simp = rules(p0, exp.p_high, 5)
    orders(e_trap, 2.0)
    orders(e_simp, 4.0)
    p1, p2 = exp.p_high * math.exp(-1.0e-6), exp.p_low * math.exp(1.0e-6)
    s1, s2 = exp.state_at(p1), exp.state_at(p2)
    assert s1.region == "plateau" and s2.region == "plateau"
    assert rel(s2.h_J_per_kg - s1.h_J_per_kg, (1.0 / M_MOL) * RU * T_M * math.log(p2 / p1)) <= T1_TOL


# ---------------------------------------------------------------------------------------------------------- stage F
@pytest.mark.fast
@pytest.mark.parametrize("eps", [1.5, 4.0, 10.0])
def test_f01_l3_the_choked_supersonic_exit_and_the_thrust_are_the_closed_forms(eps):
    gamma, T0, p0, p_a, eta, A_t = 1.4, 3000.0, 5.0e6, 1.0e5, 0.97, 2.0e-3
    _, exp = ideal_expansion(gamma, T0, p0)
    th = find_throat(exp, n_grid=N_GRID, kink_dlnp=KINK, ln_p_xtol=XTOL)
    ex = exit_state(exp, th, eps=eps, p_a=p_a, eta=eta, A_t=A_t, n_grid=N_GRID)
    p_e, u_e = closed_exit(gamma, T0, p0, eps)
    _, G_star = closed_throat(gamma, T0, p0)
    mdot = G_star * A_t
    assert ex.classification == "choked_supersonic" and ex.choked is True
    assert rel(ex.p_e_Pa, p_e) <= T1_TOL and rel(ex.u_e_m_per_s, u_e) <= T1_TOL and rel(ex.mdot_kg_per_s, mdot) <= T1_TOL
    assert rel(ex.F_N, eta * mdot * u_e + (p_e - p_a) * eps * A_t) <= T1_TOL
    p_sub, _ = closed_exit(gamma, T0, p0, eps, supersonic=False)
    assert rel(ex.p_sub_Pa, p_sub) <= T1_TOL and p_a < ex.p_sub_Pa


@pytest.mark.fast
def test_f02_the_subsonic_exit_has_p_e_equal_to_p_a_and_mdot_from_g_at_p_a():
    gamma, T0, p0, eps, eta, A_t = 1.4, 3000.0, 5.0e6, 2.0, 0.97, 2.0e-3
    _, exp = ideal_expansion(gamma, T0, p0)
    th = find_throat(exp, n_grid=N_GRID, kink_dlnp=KINK, ln_p_xtol=XTOL)
    p_sub, _ = closed_exit(gamma, T0, p0, eps, supersonic=False)
    p_a = 0.5 * (p_sub + p0)
    ex = exit_state(exp, th, eps=eps, p_a=p_a, eta=eta, A_t=A_t, n_grid=N_GRID)
    M = math.sqrt(2.0 / (gamma - 1.0) * ((p0 / p_a) ** ((gamma - 1.0) / gamma) - 1.0))
    mdot = p0 * eps * A_t * math.sqrt(gamma / (R_GAS * T0)) * M * (1.0 + 0.5 * (gamma - 1.0) * M * M) ** (-(gamma + 1.0) / (2.0 * (gamma - 1.0)))
    u_e = M * math.sqrt(gamma * R_GAS * T0 / (1.0 + 0.5 * (gamma - 1.0) * M * M))
    assert ex.classification == "unchoked_subsonic" and ex.choked is False and ex.p_e_Pa == p_a
    assert rel(ex.mdot_kg_per_s, mdot) <= T1_TOL and rel(ex.u_e_m_per_s, u_e) <= T1_TOL and rel(ex.F_N, eta * mdot * u_e) <= T1_TOL


@pytest.mark.fast
def test_f03_the_choke_threshold_separates_the_two_regimes_and_the_flow_is_continuous():
    gamma, T0, p0, eps, A_t = 1.4, 3000.0, 5.0e6, 3.0, 2.0e-3
    _, exp = ideal_expansion(gamma, T0, p0)
    th = find_throat(exp, n_grid=N_GRID, kink_dlnp=KINK, ln_p_xtol=XTOL)
    p_sub = exit_state(exp, th, eps=eps, p_a=1.0e5, eta=1.0, A_t=A_t, n_grid=N_GRID).p_sub_Pa
    above = exit_state(exp, th, eps=eps, p_a=p_sub * (1.0 + 1.0e-11), eta=1.0, A_t=A_t, n_grid=N_GRID)
    below = exit_state(exp, th, eps=eps, p_a=p_sub * (1.0 - 1.0e-11), eta=1.0, A_t=A_t, n_grid=N_GRID)
    assert above.classification == "unchoked_subsonic" and below.classification == "choked_supersonic"
    assert rel(above.mdot_kg_per_s, below.mdot_kg_per_s) <= T1_TOL


@pytest.mark.fast
def test_f04_an_area_ratio_of_one_is_sonic_and_no_pressure_difference_is_no_flow():
    _, exp = ideal_expansion(1.4)
    th = find_throat(exp, n_grid=N_GRID, kink_dlnp=KINK, ln_p_xtol=XTOL)
    sonic = exit_state(exp, th, eps=1.0, p_a=1.0e5, eta=1.0, A_t=1.0e-3, n_grid=N_GRID)
    assert sonic.classification == "choked_sonic" and sonic.p_e_Pa == th.p_t_Pa and sonic.u_e_m_per_s == th.state.u_m_per_s
    assert rel(sonic.mdot_kg_per_s, th.G_star * 1.0e-3) <= T1_TOL
    assert exit_state(exp, th, eps=2.0, p_a=exp.p_0, eta=1.0, A_t=1.0e-3, n_grid=N_GRID).classification == "no_flow"
    assert exit_state(exp, th, eps=2.0, p_a=2.0 * exp.p_0, eta=1.0, A_t=1.0e-3, n_grid=N_GRID).mdot_kg_per_s == 0.0
    with pytest.raises(StateRecoveryError):
        exit_state(exp, th, eps=0.9, p_a=1.0e5, eta=1.0, A_t=1.0e-3, n_grid=N_GRID)


@pytest.mark.fast
def test_f05_the_sound_speed_helpers_are_the_closed_forms():
    gas, exp = ideal_expansion(1.4)
    th = find_throat(exp, n_grid=N_GRID, kink_dlnp=KINK, ln_p_xtol=XTOL)
    assert rel(sound_speed_single_phase(exp, th.state), math.sqrt(1.4 * R_GAS * th.state.T_K)) <= T1_TOL
    gas2, exp2 = isothermal_expansion(0.4, alpha=0.3)
    st = exp2.state_at(4.0e6)
    a_lo, a_hi = sound_speed_plateau_sides(exp2, st, dlnp=0.05)
    n = gas2.n(4.0e6)
    assert st.region == "plateau" and rel(a_lo, math.sqrt(RU * n * T_M / (1.0 - 0.3))) <= 1.0e-6 and rel(a_hi, a_lo) <= 1.0e-6
    with pytest.raises(StateRecoveryError):
        sound_speed_single_phase(exp2, st)
    with pytest.raises(StateRecoveryError):
        sound_speed_plateau_sides(exp, th.state, dlnp=0.05)


# ---------------------------------------------------------------------------------------------------------- stage G
def _isothermal_G(p, p0):
    return p * math.sqrt(2.0 * R_GAS * T_M * math.log(p0 / p)) / (R_GAS * T_M)


@pytest.mark.fast
def test_g01_the_fold_of_the_supersonic_branch_is_detected_at_the_plateau_exit():
    p0 = 5.0e6
    _, exp = isothermal_expansion(0.45, p0)
    th = find_throat(exp, n_grid=N_GRID, kink_dlnp=KINK, ln_p_xtol=XTOL)
    lim = supersonic_limit(exp, th, n_grid=N_GRID, kink_dlnp=KINK, ln_p_xtol=XTOL)
    assert lim.classification == "kink_at_plateau_exit" and lim.p_Pa == exp.p_low
    assert rel(lim.G_min, _isothermal_G(exp.p_low, p0)) <= T1_TOL and rel(lim.eps_max, th.G_star / lim.G_min) <= T1_TOL
    assert exp.state_at(exp.p_low * math.exp(-KINK)).G > lim.G_min and exp.state_at(exp.p_low * math.exp(KINK)).G > lim.G_min


@pytest.mark.fast
def test_g02_l2_an_area_ratio_inside_the_plateau_is_the_isothermal_closed_form_and_beyond_the_fold_is_out_of_range():
    p0, eta, A_t, p_a = 5.0e6, 1.0, 1.0e-3, 1.0e3
    _, exp = isothermal_expansion(0.45, p0)
    th = find_throat(exp, n_grid=N_GRID, kink_dlnp=KINK, ln_p_xtol=XTOL)
    lim = supersonic_limit(exp, th, n_grid=N_GRID, kink_dlnp=KINK, ln_p_xtol=XTOL)
    eps = 1.0 + 0.5 * (lim.eps_max - 1.0)
    ex = exit_state(exp, th, eps=eps, p_a=p_a, eta=eta, A_t=A_t, n_grid=N_GRID)
    p_e = math.exp(brentq(lambda x: th.G_star / _isothermal_G(math.exp(x), p0) - eps, math.log(exp.p_low), math.log(th.p_t_Pa),
                          xtol=1.0e-15, rtol=1.0e-15))
    assert ex.classification == "choked_supersonic" and ex.state.region == "plateau" and rel(ex.p_e_Pa, p_e) <= T1_TOL
    beyond = exit_state(exp, th, eps=1.0 + 1.5 * (lim.eps_max - 1.0), p_a=p_a, eta=eta, A_t=A_t, n_grid=N_GRID)
    assert beyond.classification == "out_of_range" and beyond.p_e_Pa is None and beyond.F_N is None and beyond.mdot_kg_per_s is None
    assert beyond.fold_p_Pa is not None and abs(math.log(beyond.fold_p_Pa / exp.p_low)) < 0.1


@pytest.mark.fast
def test_g03_a_branch_without_a_fold_in_the_range_is_reported_as_such():
    _, exp = ideal_expansion(1.4)
    th = find_throat(exp, n_grid=N_GRID, kink_dlnp=KINK, ln_p_xtol=XTOL)
    lim = supersonic_limit(exp, th, n_grid=N_GRID, kink_dlnp=KINK, ln_p_xtol=XTOL)
    assert lim.classification == "none_in_range" and lim.eps_max > 1.0 and lim.p_Pa == exp.p_floor


# ----------------------------------------------------------------------------------------------------- reverse scenarios
@pytest.mark.fast
def test_r01_a_floor_above_the_throat_gives_not_found_and_the_exit_stage_reports_no_throat():
    gas = Gas()
    s_0, h_0 = _stagnation(gas)
    exp = Expansion(gas, B_E, s_0, h_0, 5.0e6, 4.5e6, (100.0, 5000.0))
    th = find_throat(exp, n_grid=N_GRID, kink_dlnp=KINK, ln_p_xtol=XTOL)
    assert th.classification == "not_found" and th.p_t_Pa is None and th.G_star is None
    ex = exit_state(exp, th, eps=2.0, p_a=1.0e5, eta=1.0, A_t=1.0e-3, n_grid=N_GRID)
    assert ex.classification == "no_throat" and ex.F_N is None
    assert supersonic_limit(exp, th, n_grid=N_GRID, kink_dlnp=KINK, ln_p_xtol=XTOL).classification == "no_throat"


@pytest.mark.fast
def test_r02_a_provider_without_entropy_or_a_bracket_without_the_root_is_a_classification_not_a_value():
    for gas in (Gas(with_entropy=False),):
        exp = Expansion(gas, B_E, 5000.0, 1.0e6, 5.0e6, 5.0e2, (100.0, 5000.0))
        th = find_throat(exp, n_grid=N_GRID, kink_dlnp=KINK, ln_p_xtol=XTOL)
        assert th.classification == "state_error" and th.G_star is None
    gas = Gas()
    s_0, h_0 = _stagnation(gas)
    narrow = Expansion(gas, B_E, s_0, h_0, 5.0e6, 5.0e2, (2900.0, 3100.0))
    th = find_throat(narrow, n_grid=N_GRID, kink_dlnp=KINK, ln_p_xtol=XTOL)
    assert th.classification == "state_error" and "bracket" in th.message
    assert narrow.state_at(1.0e5).status == "bracket_below" and narrow.state_at(2.0 * 5.0e6).status == "above_stagnation"


@pytest.mark.fast
def test_r03_a_grid_that_ends_before_the_area_ratio_is_not_found():
    _, exp = ideal_expansion(1.4, p_floor=2.0e6)
    th = find_throat(exp, n_grid=N_GRID, kink_dlnp=KINK, ln_p_xtol=XTOL)
    assert th.classification == "single_phase"
    ex = exit_state(exp, th, eps=10.0, p_a=1.0e5, eta=1.0, A_t=1.0e-3, n_grid=N_GRID)
    assert ex.classification == "not_found" and ex.F_N is None


@pytest.mark.fast
def test_r03b_the_not_found_message_names_p_floor_the_requested_ratio_and_the_last_grid_point():
    """Synthetic ideal gas (ASSUMED INPUT, no CEA)."""
    _, exp = ideal_expansion(1.4, p_floor=2.0e6)
    th = find_throat(exp, n_grid=N_GRID, kink_dlnp=KINK, ln_p_xtol=XTOL)
    ex = exit_state(exp, th, eps=10.0, p_a=1.0e5, eta=1.0, A_t=1.0e-3, n_grid=N_GRID)
    msg = ex.message
    assert msg.startswith("not_found: the supersonic branch does not reach the area ratio on the grid")
    assert f"p_floor = {exp.p_floor!r} Pa" in msg and "the requested area ratio is 10.0" in msg
    last_ratio = th.G_star / exp._G(exp.p_floor)
    assert f"last grid point: p = {exp.p_floor!r} Pa, T = " in msg and f"area ratio = {last_ratio!r}" in msg
    assert last_ratio < 10.0 and msg.count("last grid point") == 1


@pytest.mark.fast
def test_r04_a_chamber_that_does_not_meet_the_plateau_in_the_range_is_classified_by_its_side():
    gas = Gas(T_m=T_M, report_gamma_s=1.4)
    s_0, h_0 = _stagnation(gas, 3000.0, 5.0e6)
    above = Expansion(gas, B_E, s_0, h_0, 5.0e6, 2.0e6, (100.0, 5000.0), T_m=T_M, d_relative=D_REL)
    assert above.plateau == "absent_above" and find_throat(above, n_grid=N_GRID, kink_dlnp=KINK, ln_p_xtol=XTOL).classification == "single_phase_above_plateau"
    s_1, h_1 = _stagnation(gas, 1500.0, 5.0e6)
    below = Expansion(gas, B_E, s_1, h_1, 5.0e6, 5.0e2, (100.0, 5000.0), T_m=T_M, d_relative=D_REL)
    assert below.plateau == "absent_below" and find_throat(below, n_grid=N_GRID, kink_dlnp=KINK, ln_p_xtol=XTOL).classification == "single_phase_below_plateau"


@pytest.mark.fast
def test_r05_the_settings_have_no_default_and_invalid_inputs_raise():
    gas, exp = ideal_expansion(1.4)
    with pytest.raises(TypeError):
        find_throat(exp, kink_dlnp=KINK, ln_p_xtol=XTOL)
    with pytest.raises(TypeError):
        find_throat(exp, n_grid=N_GRID, ln_p_xtol=XTOL)
    with pytest.raises(TypeError):
        find_throat(exp, n_grid=N_GRID, kink_dlnp=KINK)
    for kw in (dict(n_grid=2), dict(n_grid=10.5), dict(kink_dlnp=0.0), dict(ln_p_xtol=-1.0), dict(ln_p_xtol=math.nan)):
        args = dict(n_grid=N_GRID, kink_dlnp=KINK, ln_p_xtol=XTOL)
        args.update(kw)
        with pytest.raises(StateRecoveryError):
            find_throat(exp, **args)
    th = find_throat(exp, n_grid=N_GRID, kink_dlnp=KINK, ln_p_xtol=XTOL)
    with pytest.raises(TypeError):
        exit_state(exp, th, eps=2.0, p_a=1.0e5, eta=1.0, A_t=1.0e-3)
    for kw in (dict(eps=math.nan), dict(p_a=-1.0), dict(eta=0.0), dict(A_t=0.0), dict(n_grid=1)):
        args = dict(eps=2.0, p_a=1.0e5, eta=1.0, A_t=1.0e-3, n_grid=N_GRID)
        args.update(kw)
        with pytest.raises(StateRecoveryError):
            exit_state(exp, th, **args)
    s_0, h_0 = _stagnation(gas)
    for args in ((s_0, h_0, 5.0e6, 5.0e6, (100.0, 5000.0)), (s_0, h_0, 5.0e6, 6.0e6, (100.0, 5000.0)),
                 (math.nan, h_0, 5.0e6, 5.0e2, (100.0, 5000.0)), (s_0, h_0, 5.0e6, 5.0e2, (500.0, 400.0))):
        with pytest.raises(StateRecoveryError):
            Expansion(gas, B_E, *args)
    with pytest.raises(StateRecoveryError):
        Expansion(gas, B_E, s_0, h_0, 5.0e6, 5.0e2, (100.0, 5000.0), T_m=T_M)                # d_relative is required with T_m
    with pytest.raises(StateRecoveryError):
        exp.state_at(-1.0)


@pytest.mark.fast
def test_r06_a_probe_step_that_leaves_the_bracketing_interval_raises():
    _, exp = _kink_expansion()
    with pytest.raises(StateRecoveryError):
        find_throat(exp, n_grid=N_GRID, kink_dlnp=0.2, ln_p_xtol=XTOL)
