"""A-6 stage 2, U-H25-5 group H-I (stages H, I; reports 283, 284; decision 117; PRE_REGISTRATION_R551): the numerical nozzle flow in
an energy run (gribs.energy.numerical_path, gribs.energy.isentrope_flow.find_throat_near / exit_state_near) with SYNTHETIC providers
(every number is an ASSUMED INPUT, not data; no CEA).  No assignment and no monkeypatch.

Tolerances (decision 117 item 2): continuation against the cold search: p_t 1e-6, G* and p_e 1e-9 (relative); the closed form against
the numerical flow in the limit of a single-phase ideal gas without a plateau: mdot, F and p_e 1e-8 (T1); the entropy added to the
two plateau providers: I2 the observed order of the trapezoid for dh = v dp (at least 2 - 0.3, T2), I7 h_L - h_S = T_m (s_L - s_S) 1e-8.
The constants of the numerical settings are asserted as written in decision 117.  Any check that does not pass is not loosened.
"""
import dataclasses
import functools
import math
from dataclasses import asdict
from types import SimpleNamespace

import pytest

import gribs.config.loader
from gribs.config import Config
from gribs.energy import recovery
from gribs.energy.exceptions import ExitOutOfRangeError, StateRecoveryError
from gribs.energy.isentrope import default_d_relative
from gribs.energy.isentrope_flow import Expansion, find_throat
from gribs.energy.numerical_path import (DEFAULT_SETTINGS, KINK_WARNING_CODE, SEPARATION_WARNING_CODE, NumericalPath,
                                         NumericalPathSettings)
from gribs.energy.nozzle_from_state import nozzle_state_from_chamber
import test_energy_integrator as tei
import test_energy_plateau_integration as tpl
import test_energy_v0_convention as v0
from test_energy_isentrope import B_E, RU, T1_TOL, T_M, Gas, orders, rel

R_GAS = RU / 24.0
N_GRID, KINK, XTOL = 200, 1.0e-3, 1.0e-10
#: the ladder of the lower end of T_bracket with 300 K added: the synthetic chambers below 2300 K have s_0 below s(500 K, p_floor)
LADDER_300 = dataclasses.replace(DEFAULT_SETTINGS, T_ladder_K=(300.0,) + DEFAULT_SETTINGS.T_ladder_K)


def make_config(**over):
    c0 = gribs.config.loader.configuration_from_document(
        {**tpl.G1, "two_phase_model": {**tpl.G1["two_phase_model"], "mode": "homogeneous_equilibrium"}, "energy_model": tpl.ENERGY}).config
    sub = type("ConfigWithEnergy", (Config,), {"energy": c0.energy})
    return sub(**{**asdict(c0), "rtol": 1.0e-8, "R_p": 0.0045, **over})


def numeric_path(provider, settings=DEFAULT_SETTINGS, **over):
    c = make_config(**over)
    return c, NumericalPath(c, SimpleNamespace(provider=provider), settings)


NOT_PLATEAU = SimpleNamespace(on_plateau=False, phi=None)
At0 = 1.0e-3


def call(path, c, T, p, eps, rec=NOT_PLATEAU, gamma=1.4, key=None):
    """One call of the dispatcher for an ideal-gas chamber (Yg = 1, R_eff = Ru / 24, gamma_s = gamma)."""
    return path.nozzle(rec=rec, p_stag=p, T_stag=T, R_eff=R_GAS, gamma_s=gamma, Yg=1.0, At=At0, Ae=eps * At0, c=c, B_e=B_E,
                       key=key if key is not None else (T, p, eps), t=0.0)


def closed(c, T, p, eps, gamma=1.4):
    return nozzle_state_from_chamber(p, T, R_GAS, gamma, 1.0, At0, eps * At0, c)


@functools.lru_cache(maxsize=None)
def case_b():
    return tpl.run_case()


# ------------------------------------------------------------------------------------------------------- the settings
@pytest.mark.fast
def test_n01_the_default_settings_are_those_of_decision_117():
    s = NumericalPathSettings()
    assert (s.n_grid, s.p_floor_ratio, s.ln_p_xtol, s.kink_dlnp) == (400, 1.0e-3, 1.0e-10, 1.0e-4)
    assert (s.plateau_dlnp, s.gate_band, s.bracket_width) == (0.025, 0.1, 0.05)
    assert s.T_ladder_K == (500.0, 800.0, 1000.0, 1200.0, 1500.0, 1800.0) and s.continuation is True
    assert DEFAULT_SETTINGS == s


@pytest.mark.fast
def test_n02_d_relative_is_the_one_of_the_plateau_bounds_of_the_recovery():
    """The temperatures at which recovery.plateau_bounds calls the provider are T_m -+ d_relative T_m with the d_relative of the run."""
    c, path = numeric_path(Gas(T_m=T_M))
    seen = []

    class Recorder:
        def tp(self, T, p, B_e):
            seen.append(T)
            return tpl.PlateauProvider().tp(T, p, B_e)

    tr = c.energy.condensed.transitions[0]
    recovery.plateau_bounds(1.0e-3, {"C": 0.3, "H": 0.7}, tr, Recorder(), c.energy.condensed, c.energy.recovery)
    d = path.d_relative
    assert d == default_d_relative(c.energy.recovery.root_tolerance)
    assert (tr.melting_temperature_K - d * tr.melting_temperature_K) in seen and (tr.melting_temperature_K + d * tr.melting_temperature_K) in seen


# ------------------------------------------------------------------------ I2 and I7 of the entropy of the plateau providers
PROVIDERS = [("plateau_integration", tpl.PlateauProvider, tpl.T_M, tpl.LATENT), ("v0_convention", v0.PlateauProvider, v0.T_M, v0.H_L_REF - v0.H_S_REF)]


def _h(st, T):
    return st.u_J_per_kg + st.n_gas_per_kg * RU * T


@pytest.mark.fast
@pytest.mark.parametrize("name,cls,T_m,L", PROVIDERS)
def test_n03_i7_the_enthalpy_jump_is_t_m_times_the_entropy_jump(name, cls, T_m, L):
    prov, delta = cls(), 1.0e-9 * T_m
    for p in (1.0e6, 3.0e6):
        lo, hi = prov.tp(T_m - delta, p, B_E), prov.tp(T_m + delta, p, B_E)
        assert rel(_h(hi, T_m + delta) - _h(lo, T_m - delta), T_m * (hi.entropy_J_per_kgK - lo.entropy_J_per_kgK)) <= T1_TOL
        assert rel(T_m * (hi.entropy_J_per_kgK - lo.entropy_J_per_kgK), 0.3 * L) <= 1.0e-6     # the jump per kg of products is Yc L (plus the 2 delta of sensible heat)


@pytest.mark.fast
@pytest.mark.parametrize("name,cls,T_m,L", PROVIDERS)
def test_n04_i2_dh_equals_v_dp_above_the_plateau_with_the_observed_order_of_the_trapezoid(name, cls, T_m, L):
    prov, p0 = cls(), 5.0e6
    st = prov.tp(3000.0, p0, B_E)
    exp = Expansion(prov, B_E, st.entropy_J_per_kgK, _h(st, 3000.0), p0, 1.0e3, (100.0, 5000.0), T_m=T_m, d_relative=1.0e-9)
    assert exp.has_plateau and 1.0e3 < exp.p_low < exp.p_high < p0
    errors = []
    for k in range(5):
        n = 4 * 2 ** k
        xs = [math.log(p0) + (math.log(exp.p_high) - math.log(p0)) * i / n for i in range(n + 1)]
        pts = [exp.state_at(math.exp(x)) for x in xs]
        assert all(q.status == "ok" for q in pts)
        f = [q.p_Pa * q.v_m3_per_kg for q in pts]
        dh = pts[-1].h_J_per_kg - pts[0].h_J_per_kg
        errors.append(abs((xs[1] - xs[0]) * (0.5 * f[0] + sum(f[1:-1]) + 0.5 * f[-1]) - dh) / abs(dh))
    orders(errors, 2.0)


# ------------------------------------------------------------------------------------------------------------ the gate
@pytest.mark.fast
def test_n05_without_a_declared_transition_the_closed_form_is_returned_unchanged_and_no_provider_call_is_made():
    run = tei.base_run()
    assert run["provenance"]["numerical_path_evaluations"] == 0 and run["provenance"]["numerical_path_switches"] == 0
    assert run["provenance"]["throat_kink_evaluations"] == 0 and run["provenance"]["closed_form_evaluations"] > 0
    assert run["warnings"] == []
    calls = []

    class Counting:
        def tp(self, T, p, B_e):
            calls.append((T, p))
            return Gas().tp(T, p, B_e)

    c = make_config()
    plain = dataclasses.replace(c.energy, condensed=dataclasses.replace(c.energy.condensed, transitions=(), latent_heat_check_tolerance_rel=None))
    c0 = type("C", (), {"energy": plain})()
    path = NumericalPath(c0, SimpleNamespace(provider=Counting()))
    nz, margin = call(path, c, 3000.0, 5.0e6, 2.0)
    assert nz == closed(c, 3000.0, 5.0e6, 2.0) and margin is None and calls == [] and not path.declared


@pytest.mark.fast
def test_n05b_an_exit_below_p_floor_is_not_found_and_the_message_names_p_floor():
    """The supersonic scan ends at p_floor = p_c * p_floor_ratio; the failure keeps the leading sentence of the not_found and names
    p_floor, the last grid point (p, T, area ratio) and the requested area ratio.  Synthetic ideal gas (ASSUMED INPUT)."""
    c, path = numeric_path(Gas(T_m=T_M), p_a=1.0e5)
    p, eps = 5.0e6, 1.0e4
    p_floor = p * DEFAULT_SETTINGS.p_floor_ratio
    with pytest.raises(StateRecoveryError) as caught:
        call(path, c, 3500.0, p, eps)
    msg = str(caught.value)
    assert "not_found: the supersonic branch does not reach the area ratio on the grid" in msg
    assert f"p_floor = {p_floor!r} Pa" in msg and f"the requested area ratio is {eps!r}" in msg
    assert "last grid point: p = " in msg and " K, area ratio = " in msg


@pytest.mark.fast
def test_n06_the_gate_follows_the_rule_in_each_region():
    c, path = numeric_path(Gas(T_m=T_M), p_a=1.0e5)
    p = 5.0e6
    nz, m = call(path, c, 1500.0, p, 1.7)
    assert path.last.gate_reason == "T_stag <= T_m" and not path.last.numeric and nz == closed(c, 1500.0, p, 1.7) and m is None
    call(path, c, 3500.0, p, 1.5)
    assert path.last.gate_reason == "the closed-form exit is above the band" and not path.last.numeric and path.last.gate_x > 0.1
    call(path, c, 3500.0, p, 1.55)
    assert path.last.gate_reason == "the numerical path does not enter the plateau" and not path.last.numeric and 0.0 < path.last.gate_x <= 0.1
    nz, m = call(path, c, 3500.0, p, 1.65)
    assert path.last.gate_reason == "band" and path.last.numeric and -0.1 <= path.last.gate_x < 0.0 and m is not None
    nz, m = call(path, c, 3500.0, p, 1.8)
    assert path.last.numeric and path.last.gate_x < -0.1 and path.last.gate_reason.startswith("the closed-form exit is below")
    c, path = numeric_path(Gas(T_m=T_M), LADDER_300, p_a=1.0e5)
    rec = SimpleNamespace(on_plateau=True, phi=0.4)
    nz, m = call(path, c, T_M, p, 1.0, rec=rec)
    assert path.last.gate_reason == "on_plateau" and path.last.numeric and path.last.throat_class == "inside_plateau"
    # no flow: the closed form, no gate
    nz, m = call(path, c, 3500.0, 5.0e4, 1.7)
    assert nz["regime"] == "no-flow" and m is None


@pytest.mark.fast
def test_n07_counters_and_switches_count_the_calls():
    c, path = numeric_path(Gas(T_m=T_M), p_a=1.0e5)
    for eps in (1.5, 1.7, 1.7, 1.5):                 # closed, numerical, numerical (cache hit counts), closed
        call(path, c, 3500.0, 5.0e6, eps)
    assert path.provenance() == {"numerical_path_evaluations": 2, "closed_form_evaluations": 2, "numerical_path_switches": 2,
                                 "throat_kink_evaluations": 0}


# ------------------------------------------------------------------ the numerical flow against the closed form (the limit)
@pytest.mark.fast
@pytest.mark.parametrize("eps", [1.0, 1.7, 3.0])
def test_n08_in_the_limit_of_a_single_phase_ideal_gas_the_numerical_flow_is_the_closed_form(eps):
    c, path = numeric_path(Gas(T_m=T_M), p_a=1.0e5, Cd=0.97, eta_thrust=0.95)
    T, p = 6000.0, 5.0e6                                            # the exit stays above T_m: no plateau on the path
    ref = closed(c, T, p, eps)
    r = path._numerical(NOT_PLATEAU, p, T, 1.0, At0, eps * At0, c, B_E, (T, p), None, "test")
    nz = r.nz
    assert ref["regime"] == "choked" and nz["regime"] == "choked"
    for key in ("mdot_total", "F", "pe", "ve"):
        assert rel(nz[key], ref[key]) <= T1_TOL, key
    assert r.margin > 0.0 and nz["Me"] is None and nz["A_eff"] == eps * At0


@pytest.mark.fast
def test_n09_the_unchoked_flow_has_the_ambient_exit_pressure_and_a_negative_margin_and_the_margin_changes_sign_at_the_threshold():
    T, p, eps = 3500.0, 5.0e6, 3.0
    c, path = numeric_path(Gas(T_m=T_M), p_a=4.9e6, Cd=1.0, eta_thrust=1.0)
    ref = closed(c, T, p, eps)
    r = path._numerical(NOT_PLATEAU, p, T, 1.0, At0, eps * At0, c, B_E, (T, p), None, "test")
    assert ref["regime"] == "subsonic" and r.nz["regime"] == "subsonic" and r.margin < 0.0 and r.nz["pe"] == c.p_a
    for key in ("mdot_total", "F", "ve"):
        assert rel(r.nz[key], ref[key]) <= 1.0e-6, key             # the closed form of the subsonic branch has the same gamma
    # the sign follows p_a: bracket the closed-form threshold (p_a / p_0 at the threshold)
    from gribs.energy.nozzle_from_state import _choke_threshold_pressure
    margins = []
    p_a_threshold = p / _choke_threshold_pressure(1.0, eps, 1.4)      # the chamber threshold for p_a = 1 is p_thr(1)
    for factor in (0.995, 1.005):
        c2, path2 = numeric_path(Gas(T_m=T_M), p_a=factor * p_a_threshold, Cd=1.0, eta_thrust=1.0)
        margins.append(path2._numerical(NOT_PLATEAU, p, T, 1.0, At0, eps * At0, c2, B_E, (T, p), None, "test").margin)
    assert margins[0] > 0.0 > margins[1]


# --------------------------------------------------------------- the continuation against the cold search (decision 117)
@pytest.mark.fast
def test_n10_the_continuation_agrees_with_the_cold_search_along_a_sequence_of_chambers():
    cold_settings = dataclasses.replace(DEFAULT_SETTINGS, continuation=False)
    c, warm = numeric_path(Gas(T_m=T_M), p_a=1.0e5, Cd=1.0, eta_thrust=1.0)
    _, cold = numeric_path(Gas(T_m=T_M), cold_settings, p_a=1.0e5, Cd=1.0, eta_thrust=1.0)
    eps = 1.7                                                       # a path of a hot chamber that enters the plateau
    for k in range(14):
        T, p = 3500.0 - 40.0 * k, 5.0e6 * (1.0 + 0.01 * k)
        call(warm, c, T, p, eps)
        call(cold, c, T, p, eps)
        w, k0 = warm.last, cold.last
        assert w.numeric and k0.numeric
        assert rel(w.p_t, k0.p_t) <= 1.0e-6 and rel(w.nz["mdot_total"], k0.nz["mdot_total"]) <= 1.0e-9 and rel(w.p_e, k0.p_e) <= 1.0e-9
        assert rel(w.nz["F"], k0.nz["F"]) <= 1.0e-9 and w.throat_class == k0.throat_class and w.exit_class == k0.exit_class
    assert warm.stats["throat_warm"] > 0 and warm.stats["exit_warm"] > 0 and cold.stats["throat_warm"] == 0 and cold.stats["exit_warm"] == 0


@pytest.mark.integration
def test_n11_the_continuation_agrees_with_the_cold_search_at_the_states_of_the_run_through_the_plateau():
    c, run = case_b()
    model, sol = run["model"], run["sol_burn"]
    cold_settings = dataclasses.replace(DEFAULT_SETTINGS, continuation=False)
    cold = NumericalPath(c, SimpleNamespace(provider=tpl.PlateauProvider()), cold_settings)
    warm = NumericalPath(c, SimpleNamespace(provider=tpl.PlateauProvider()))
    idx = [j for j in range(sol.t.size)]
    used = 0
    for j in idx[:: max(1, len(idx) // 30)]:
        y = sol.y[:, j]
        k = model.terms(y, True)
        args = dict(rec=k["rec"], p_stag=k["p"], T_stag=k["T"], R_eff=tpl.RU * k["rec"].n_gas_kmol / k["rec"].m_total_kg,
                    gamma_s=k["gamma_s"], Yg=k["rec"].Yg, At=k["At"], Ae=model.Ae, c=c, B_e=model.B_e(float(y[0])),
                    key=(float(y[0]), float(y[1]), k["V_free"]))
        warm.nozzle(**args)
        cold.nozzle(**args)
        w, k0 = warm.last, cold.last
        assert w.numeric == k0.numeric
        if w.numeric:
            used += 1
            assert rel(w.p_t, k0.p_t) <= 1.0e-6 and rel(w.nz["mdot_total"], k0.nz["mdot_total"]) <= 1.0e-9
            assert rel(w.p_e, k0.p_e) <= 1.0e-9 and w.throat_class == k0.throat_class
    assert used >= 5 and warm.stats["throat_warm"] > 0


# --------------------------------------------------------------------------------------------- exceptions and warnings
@pytest.mark.fast
def test_n12_an_area_ratio_beyond_the_fold_raises_the_exit_out_of_range_error_and_an_area_ratio_before_it_does_not():
    c, path = numeric_path(Gas(T_m=T_M), LADDER_300, p_a=1.0e3, Cd=1.0, eta_thrust=1.0)
    p, rec = 5.0e6, SimpleNamespace(on_plateau=True, phi=0.45)
    exp = path._expansion(rec, p, T_M, B_E, "k")
    th = find_throat(exp, n_grid=N_GRID, kink_dlnp=KINK, ln_p_xtol=XTOL)
    from gribs.energy.isentrope_flow import supersonic_limit
    lim = supersonic_limit(exp, th, n_grid=N_GRID, kink_dlnp=KINK, ln_p_xtol=XTOL)
    ok = path._numerical(rec, p, T_M, 1.0, At0, (1.0 + 0.5 * (lim.eps_max - 1.0)) * At0, c, B_E, "k", None, "test")
    assert ok.exit_class == "choked_supersonic"
    with pytest.raises(ExitOutOfRangeError) as caught:
        path._numerical(rec, p, T_M, 1.0, At0, (1.0 + 1.5 * (lim.eps_max - 1.0)) * At0, c, B_E, "k", None, "test")
    assert ExitOutOfRangeError.code == "E_ENERGY_EXIT_OUT_OF_RANGE" and "E_ENERGY_EXIT_OUT_OF_RANGE" in str(caught.value)
    assert isinstance(caught.value, StateRecoveryError)


@pytest.mark.fast
def test_n13_a_provider_without_entropy_stops_where_the_flow_needs_it_and_not_elsewhere():
    c, path = numeric_path(Gas(T_m=T_M, with_entropy=False), p_a=1.0e5)
    nz, m = call(path, c, 1500.0, 5.0e6, 1.7)                                            # T_stag <= T_m: the closed form, no entropy needed
    assert m is None and not path.last.numeric
    nz, m = call(path, c, 3500.0, 5.0e6, 1.5)                                            # above the band: the closed form
    assert m is None and not path.last.numeric
    with pytest.raises(StateRecoveryError, match="entropy"):
        call(path, c, 3500.0, 5.0e6, 1.7)                                                # below the band: numerical, no entropy
    with pytest.raises(StateRecoveryError, match="entropy"):
        call(path, c, 3500.0, 5.0e6, 1.6)                                                # in the band: the exact judgement needs it


@pytest.mark.fast
def test_n14_the_separation_warning_is_recorded_when_the_condition_holds_and_not_otherwise():
    # choked, eps > 1, use_separation, p_e < sep_ratio p_a  -> recorded once per code, regime stays choked
    c, path = numeric_path(Gas(T_m=T_M), p_a=2.0e6, use_separation=True, sep_ratio=0.4)
    nz, m = call(path, c, 3500.0, 5.0e6, 1.7)
    assert path.last.numeric and path.last.p_e < 0.4 * c.p_a and nz["regime"] == "choked"
    w = path.warnings()
    assert [x["code"] for x in w] == [SEPARATION_WARNING_CODE] and w[0]["severity"] == "warning" and w[0]["context"]["evaluations"] == 1
    assert 0.0 < w[0]["context"]["min_pe_over_pa"] < 0.4 and w[0]["context"]["first_time_s"] == 0.0
    for kw, eps in ((dict(use_separation=False), 1.7), (dict(use_separation=True), 1.0)):          # no separation model / eps = 1
        c2, p2 = numeric_path(Gas(T_m=T_M), p_a=2.0e6, sep_ratio=0.4, **kw)
        call(p2, c2, 3500.0, 5.0e6, eps)
        assert p2.last.numeric is (eps > 1.0) and p2.warnings() == []
    c3, p3 = numeric_path(Gas(T_m=T_M), p_a=1.0e5, use_separation=True, sep_ratio=0.4)       # p_e above 0.4 p_a: no warning
    call(p3, c3, 3500.0, 5.0e6, 1.7)
    assert p3.last.numeric and p3.last.p_e >= 0.4 * c3.p_a and p3.warnings() == []


@pytest.mark.fast
def test_n15_a_throat_at_the_plateau_entry_is_recorded_as_one_warning_with_the_count_and_the_first_time():
    c, path = numeric_path(Gas(T_m=T_M), LADDER_300, p_a=1.0e5, Cd=1.0, eta_thrust=1.0)
    for k in range(3):
        path.nozzle(rec=NOT_PLATEAU, p_stag=5.0e6, T_stag=2340.0 + 0.001 * k, R_eff=R_GAS, gamma_s=1.4, Yg=1.0, At=At0, Ae=1.0 * At0, c=c,
                    B_e=B_E, key=("kink", k), t=0.5 + k)
        assert path.last.numeric and path.last.throat_class == "kink_at_plateau_entry"
    w = path.warnings()
    assert [x["code"] for x in w] == [KINK_WARNING_CODE] and w[0]["context"] == {"evaluations": 3, "first_time_s": 0.5}
    assert "independent confirmation" in w[0]["message"] and "report 282" in w[0]["message"]
    assert path.provenance()["throat_kink_evaluations"] == 3


# -------------------------------------------------------------------------------- the runs of the existing synthetic tests
@pytest.mark.integration
def test_n16_the_run_through_the_plateau_records_the_numerical_path_in_the_provenance_and_the_warning():
    c, run = case_b()
    prov = run["provenance"]
    assert prov["numerical_path_evaluations"] > 0 and prov["closed_form_evaluations"] > 0 and prov["numerical_path_switches"] > 0
    assert prov["gamma_s_min"] == prov["gamma_s_max"] == 1.2                             # the states off the plateau
    assert prov["gamma_s_plateau_interpolated_min"] is not None and prov["gamma_s_plateau_interpolated_max"] is not None
    assert abs(prov["gamma_s_plateau_interpolated_min"] - 1.2) < 1.0e-9 and abs(prov["gamma_s_plateau_interpolated_max"] - 1.2) < 1.0e-9
    assert run["stop_burn"] == "burnout"
    codes = [w["code"] for w in run["warnings"]]
    assert codes == ([KINK_WARNING_CODE] if prov["throat_kink_evaluations"] else [])
    assert run["integration_budget"]["rhs_calls"] < 2000


@pytest.mark.integration
def test_n17_the_gamma_range_is_none_when_no_state_is_off_a_plateau():
    # the provenance of a run whose declared states never left a plateau cannot be made cheaply; the helper is checked directly
    import gribs.simulation.energy_integrator as ei
    assert ei._finite_or_none(math.inf) is None and ei._finite_or_none(-math.inf) is None and ei._finite_or_none(1.2) == 1.2


# ------------------------------------------------------------------------------------- Yg rounded at the entry (H-26)
@pytest.mark.fast
def test_round_gas_mass_fraction():
    from gribs.energy.numerical_path import YG_ROUNDOFF_TOL, round_gas_mass_fraction
    assert YG_ROUNDOFF_TOL == 1.0e-12
    assert round_gas_mass_fraction(1.0000000000000002) == 1.0
    assert round_gas_mass_fraction(1.0 + 1.0e-9) == 1.0 + 1.0e-9
    assert round_gas_mass_fraction(0.85) == 0.85
    assert round_gas_mass_fraction(0.0) == 0.0


@pytest.mark.fast
def test_nozzle_entry_rounds_yg_above_one_by_roundoff():
    """Yg = 1 + 1 ulp reaches nozzle_state_from_chamber (p_stag > p_a, closed form, T_stag <= T_m) without ThermochemistryError."""
    yg = 1.0000000000000002
    c, path = numeric_path(None)
    nz, margin = path.nozzle(rec=NOT_PLATEAU, p_stag=2.0e6, T_stag=1500.0, R_eff=R_GAS, gamma_s=1.4, Yg=yg, At=At0, Ae=4.0 * At0,
                             c=c, B_e=B_E, key=("yg", 1), t=0.0)
    assert nz["regime"] == "choked" and nz["gas_mass_fraction"] == 1.0 and margin is None
    ref = closed(c, 1500.0, 2.0e6, 4.0)
    assert nz["mdot"] == ref["mdot"] and nz["F"] == ref["F"]
