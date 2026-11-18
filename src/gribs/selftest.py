"""Numerical self-tests (A-6 unit 8). Moved verbatim from section 9 of the entry: the synthetic analytic test
backend (_SyntheticTestBackend), the derivative and single-phase-limit tolerances, the exact self-test name
contract, self_tests (run by --selftest) and self_tests_post_run (the post-run conservation verdict printed
after a production run).

Tests patch self_tests_post_run here ("gribs.selftest.self_tests_post_run"); the entry re-exports every name
and reads self_tests_post_run as gribs.selftest.self_tests_post_run, so a patch on the entry would not reach
it. Names of other package modules that tests patch are read as <module>.<name> (ADR decision 1); this
includes the classes PressureTableBackend and _ThermoSwap of gribs.thermo, the only lines that differ from
the entry text.
"""
from __future__ import annotations

import copy
import math
import tempfile
from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, Optional, Sequence, Tuple

import numpy as np

import gribs.nozzle
import gribs.simulation.integrator
import gribs.simulation.postprocess
import gribs.thermo
from gribs.config import Config, geometry, web_thickness
from gribs.constants import R_UNIVERSAL, TWO_PHASE_HOMOGENEOUS, TWO_PHASE_SINGLE_LEGACY
from gribs.nozzle import area_mach, critical_ratio, mach_from_area
from gribs.simulation.integrator import IP, make_rhs
from gribs.simulation.postprocess import post_run_consistency_report
from gribs.thermo import theta_and_deriv


#: Documented tolerance of the analytic-vs-finite-difference derivative check.
#: The chamber-property tables use a C1-continuous PCHIP spline whose second
#: derivative jumps at the knots, so a centred difference with h = 1e-5 p is
#: limited to O(h) there (1e-5 relative tolerance).
DERIVATIVE_TOLERANCE_TABLE = 1.0e-5


class _SyntheticTestBackend(gribs.thermo.PressureTableBackend):
    """Synthetic analytic property table for numerical self-tests ONLY.

    Not selectable from any configuration; used to verify:
      * the single-phase limit (Yg = 1 everywhere) reproduces the legacy
        Theta model bit-for-bit (self-test: single-phase limit), and
      * the homogeneous-equilibrium run matches the legacy run when the table
        contains no condensed phase (compatibility test).
    Constant T0, gamma_s, M and a configurable constant condensed fraction are
    sufficient because these tests compare MODELS, not chemistry.
    """

    name = "synthetic_selftest"
    uses_cea = False
    supports_two_phase = True
    table_stem = "thermo_synthetic_selftest"
    equilibrium_formulation = ("synthetic constant-property analytic table "
                               "(numerical self-tests only, never a production "
                               "backend)")

    def __init__(self, c: Config, condensed_fraction: float = 0.0):
        if not (0.0 <= condensed_fraction < 1.0):
            raise ValueError("condensed_fraction must lie in [0, 1).")
        self._ycond = float(condensed_fraction)
        # cache must never interfere with the self-tests
        c2 = Config(**{**asdict(c), "cea_cache_enabled": False,
                       "cea_rebuild_cache": False})
        super().__init__(c2, c.reactants_spec_list, c.packing_fraction,
                         Path(tempfile.gettempdir()) / "gribs_selftest_cache")

    def _build_table(self) -> Dict[str, Any]:
        c = self.c
        pmin = max(c.p_fit_min / 50.0, 1.0e3)
        pmax = c.p_fit_max * 10.0
        pressures = np.geomspace(pmin, pmax, 41)
        T0 = 3000.0
        gamma = 1.2
        M = 25.0
        R = R_UNIVERSAL / M
        Yg = 1.0 - self._ycond
        n = pressures.size
        return {
            "schema": self.cache_schema,
            "backend": self.name,
            "cache_key": self.key,
            "source": "synthetic analytic table (self-tests only)",
            "pressure_Pa": pressures.tolist(),
            "temperature_K": [T0] * n,
            "gamma_s": [gamma] * n,
            "equilibrium_molecular_weight_kg_kmol": [M] * n,
            "effective_gas_constant_J_kgK": [R] * n,
            "gas_mass_fraction": [Yg] * n,
            "condensed_mass_fraction": [self._ycond] * n,
            "psi_J_kg": [R * T0 * c.eta_T0] * n,
        }

    def metadata(self) -> Dict[str, Any]:
        return {"backend": self.name, "synthetic_selftest_backend": True}


#: Tolerance of the single-phase-limit regression test.  With Yg = 1 the two
#: code paths evaluate bit-identical arithmetic, so any difference beyond
#: double-rounding level indicates a structural inconsistency.
SINGLE_PHASE_LIMIT_TOLERANCE = 1.0e-12


SELF_TEST_COMMON_NAMES = (
    "analytic vs numerical dTheta/dp",
    "geometry consistency dVg/dx = Ab",
    "area-Mach inversion round trip",
    "closed-form choked thrust identity",
    "finite initial derivatives",
    "chamber-property table physical",
    "phase fractions sum to 1 (table)",
    "phase fractions physical (table)",
    "Psi = R_eff*T0 (table + interpolator)",
    "analytic vs numerical dPsi/dp",
    "pressure factor 1 - p*Psi'/Psi > 0",
    "single-phase limit: model B == Theta model",
    "nozzle decomposition mdot_total = gas + condensed",
    "no-condensed compatibility (hom. vs legacy run)",
)


def expected_self_test_names(c: Config) -> Tuple[str, ...]:
    """Return the exact ordered self-test names for a configuration."""
    if math.isclose(
        c.eps_nozzle,
        1.0,
        rel_tol=0.0,
        abs_tol=1.0e-12,
    ):
        nozzle_name = "converging-nozzle continuity (mdot, F)"
    else:
        nozzle_name = "C-D nozzle transition mass-flow continuity"

    return (
        SELF_TEST_COMMON_NAMES[0],
        SELF_TEST_COMMON_NAMES[1],
        SELF_TEST_COMMON_NAMES[2],
        nozzle_name,
        *SELF_TEST_COMMON_NAMES[3:],
    )


def self_test_name_contract(
    c: Config,
    actual_names: Sequence[str],
) -> dict:
    """Compare reported self-test names with the exact expected contract."""
    expected = expected_self_test_names(c)
    actual = tuple(actual_names)

    missing = [
        name
        for name in expected
        if name not in actual
    ]
    unexpected = [
        name
        for name in actual
        if name not in expected
    ]
    duplicates = sorted({
        name
        for name in actual
        if actual.count(name) > 1
    })

    return {
        "expected_names": list(expected),
        "actual_names": list(actual),
        "expected_count": len(expected),
        "actual_count": len(actual),
        "missing_names": missing,
        "unexpected_names": unexpected,
        "duplicate_names": duplicates,
        "passed": bool(
            actual == expected
            and not missing
            and not unexpected
            and not duplicates
        ),
    }


def self_tests(c: Config, verbose=True) -> bool:
    """Numerical self-tests executed explicitly with ``--selftest``.

    All v0.4.3-alpha tests are kept unchanged (T1-T7); the two-phase model
    adds T8-T15 (see the test names).  Post-simulation conservation tests run
    in :func:`self_tests_post_run` after the ODE integration.
    """
    out, ok_all = [], True
    result_names = []

    def chk(name, ok, detail=""):
        nonlocal ok_all
        ok_all &= bool(ok)
        result_names.append(name)
        out.append(f"  [{'PASS' if ok else 'FAIL'}] {name:<46s} {detail}")

    backend = gribs.thermo.thermochemistry()
    two_phase_active = (c.two_phase_mode == TWO_PHASE_HOMOGENEOUS
                        and bool(getattr(backend, "supports_two_phase", False)))

    # T1 analytic dTheta/dp against a centred difference (derivative validation)
    tol_deriv = DERIVATIVE_TOLERANCE_TABLE
    e = 0.0
    for p in np.geomspace(max(c.p_fit_min * 1.05, 1.1e5), c.p_fit_max * 0.95, 25):
        hstep = p * 1e-5
        fd = (theta_and_deriv(p + hstep, c)[0] - theta_and_deriv(p - hstep, c)[0]) / (2 * hstep)
        e = max(e, abs(fd / theta_and_deriv(p, c)[1] - 1.0))
    chk("analytic vs numerical dTheta/dp", e < tol_deriv,
        f"max rel. err = {e:.2e} (tol {tol_deriv:.0e})")

    # T2 geometric consistency dVg/dx = Ab
    xw = web_thickness(c)
    e = 0.0
    for x in np.linspace(0.02 * xw, 0.98 * xw, 25):
        hstep = xw * 1e-6
        fd = (geometry(x + hstep, c)[3] - geometry(x - hstep, c)[3]) / (2 * hstep)
        e = max(e, abs(fd / geometry(x, c)[2] - 1.0))
    chk("geometry consistency dVg/dx = Ab", e < 1e-7, f"max rel. err = {e:.2e}")

    # T3 area-Mach inversion round trip
    e = 0.0
    for eps in (1.5, 2.5, 6.0, 20.0):
        for g in (1.10, 1.15, 1.25):
            e = max(e, abs(area_mach(mach_from_area(eps, g, True), g) / eps - 1.0),
                    abs(area_mach(mach_from_area(eps, g, False), g) / eps - 1.0))
    chk("area-Mach inversion round trip", e < 1e-10, f"max rel. err = {e:.2e}")

    # T4 choking-boundary continuity and transition classification.
    #
    # For a converging nozzle (Ae/At = 1), mass flow and thrust must both
    # remain continuous across the choking boundary.
    #
    # For a C-D nozzle (Ae/At > 1), the current jump model switches directly
    # between supersonic choked and subsonic solutions. Mass flow must remain
    # continuous. The thrust jump remains a diagnostic until an internal-shock
    # transition model is implemented.
    continuity_tolerance = 1.0e-6
    eps = c.eps_nozzle

    cc = Config(**{
        **asdict(c),
        "eps_nozzle": eps,
        "use_separation": False
    })

    At0 = math.pi * cc.R_t0 ** 2
    Ae = eps * At0
    pstar = gribs.nozzle.choke_limit_pressure(cc, At0, Ae)

    if pstar is None:
        chk("nozzle choking-boundary continuity", False,
            "choking-limit pressure could not be determined")
    else:
        lo = gribs.nozzle.nozzle_state(pstar * (1.0 - 1e-9), At0, Ae, cc)
        hi = gribs.nozzle.nozzle_state(pstar * (1.0 + 1e-9), At0, Ae, cc)
        worst_mdot = abs(lo["mdot"] / hi["mdot"] - 1.0)
        worst_thrust = abs((lo["F"] + 1e-12) / (hi["F"] + 1e-12) - 1.0)

        if math.isclose(eps, 1.0, rel_tol=0.0, abs_tol=1.0e-12):
            worst = max(worst_mdot, worst_thrust)
            chk("converging-nozzle continuity (mdot, F)",
                worst < continuity_tolerance,
                f"Ae/At = {eps:.6g}, mdot jump = {worst_mdot:.2e}, "
                f"F jump = {worst_thrust:.2e}")
        else:
            chk("C-D nozzle transition mass-flow continuity",
                worst_mdot < continuity_tolerance,
                f"Ae/At = {eps:.6g}, mdot jump = {worst_mdot:.2e}, "
                f"F jump = {worst_thrust:.2e} "
                "(diagnostic; known jump-model limitation)")

        # T5 closed-form check of the converging-nozzle choked thrust.
        # Mass basis: mdot_total comes from the isentropic formula evaluated
        # with R_eff and ve = sqrt(g*R_eff*Te), so the momentum term carries
        # the phase split through R_eff alone.  The pressure thrust does not:
        #   F = mdot_total*ve + (pe - pa)*At
        #     = At*(crit*p*(g + 1) - pa),  crit = (2/(g+1))^(g/(g-1))
        # This is the classical single-phase identity; no /Yg factor appears
        # because R_eff already carries Yg.
        cc5 = Config(**{**asdict(c), "eps_nozzle": 1.0, "Cd": 1.0,
                        "eta_thrust": 1.0, "use_separation": False})
        At0 = math.pi * cc5.R_t0 ** 2
        e = 0.0
        for p in np.geomspace(cc5.p_a / 0.5, 50e6, 20):
            nz = gribs.nozzle.nozzle_state(p, At0, At0, cc5)
            _, _, g = gribs.thermo.gas_props(p, cc5)
            F_ref = At0 * (critical_ratio(g) * p * (g + 1.0) - cc5.p_a)
            e = max(e, abs(nz["F"] / F_ref - 1.0))
        chk("closed-form choked thrust identity", e < 1e-12,
            f"max rel. err = {e:.2e}")

        # T6 finite and physical right-hand side at t = 0
        try:
            dy = make_rhs(c, c.eps_nozzle * math.pi * c.R_t0 ** 2, True)(
                0.0, np.array([
                    c.p0_init,
                    0.0,
                    c.R_t0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ]))
            chk("finite initial derivatives", np.all(np.isfinite(dy)),
                f"dp/dt = {dy[0]:+.3e} Pa/s, r = {dy[1]:.3e} m/s")
        except Exception as exc:                              # pragma: no cover
            chk("finite initial derivatives", False, str(exc))

    # T7 chamber-property table completeness and physicality
    try:
        pv = backend.property_validation
        ok = bool(pv.get("all_points_finite", False)) and \
            bool(pv.get("pressure_strictly_increasing", False)) and \
            float(pv.get("temperature_min_K", 0.0)) > 0.0 and \
            float(pv.get("gamma_min", 0.0)) > 1.0 and \
            float(pv.get("min_pressure_factor_1_minus_pTheta_over_Theta", -1.0)) > 0.0 and \
            float(pv.get("min_pressure_factor_1_minus_pPsi_over_Psi", -1.0)) > 0.0
        detail = (f"{pv.get('pressure_points', 'n/a')} points, "
                  f"T0 = {pv.get('temperature_min_K', float('nan')):.1f}-"
                  f"{pv.get('temperature_max_K', float('nan')):.1f} K, "
                  f"gamma = {pv.get('gamma_min', float('nan')):.4f}-"
                  f"{pv.get('gamma_max', float('nan')):.4f}")
        chk("chamber-property table physical", ok, detail)
    except Exception as exc:                                  # pragma: no cover
        chk("chamber-property table physical", False, str(exc))

    # ---- two-phase tests (new in v0.5.0-alpha) ------------------------------
    table = getattr(backend, "table", {})
    if "gas_mass_fraction" in table:
        Yg_tab = np.asarray(table["gas_mass_fraction"], dtype=float)
        Yc_tab = np.asarray(table["condensed_mass_fraction"], dtype=float)
    else:
        npts = len(table.get("pressure_Pa", [0] * 4))
        Yg_tab = np.ones(npts)
        Yc_tab = np.zeros(npts)

    # T8 phase-fraction completeness: Yg + Yc = 1 at every table point
    e = float(np.max(np.abs(Yg_tab + Yc_tab - 1.0)))
    chk("phase fractions sum to 1 (table)", e <= 1.0e-12,
        f"max |Yg+Yc-1| = {e:.2e} (tol 1e-12)")

    # T9 physical ranges: 0 < Yg <= 1, 0 <= Yc < 1
    ok = bool(np.all(Yg_tab > 0.0) and np.all(Yg_tab <= 1.0)
              and np.all(Yc_tab >= 0.0) and np.all(Yc_tab < 1.0))
    chk("phase fractions physical (table)", ok,
        f"Yg in [{Yg_tab.min():.6g}, {Yg_tab.max():.6g}], "
        f"Yc in [{Yc_tab.min():.6g}, {Yc_tab.max():.6g}]")

    # T10 Psi consistency: Psi = R_eff*T_used at every table point.
    # R_eff already carries Yg, so the solver-facing Psi and the stored
    # table column are the same quantity in both the two-phase mode and
    # the single-phase legacy mode (where Yg := 1 by definition).
    R_tab = np.asarray(table["effective_gas_constant_J_kgK"], dtype=float)
    T_tab = np.asarray(table["temperature_K"], dtype=float) * c.eta_T0
    psi_ref = R_tab * T_tab
    psi_solver_ref = psi_ref
    if "psi_J_kg" in table:
        psi_tab = np.asarray(table["psi_J_kg"], dtype=float)
    else:
        psi_tab = psi_ref
    e = float(np.max(np.abs(psi_tab - psi_ref) / np.maximum(psi_ref, 1e-300)))
    # also check the interpolated Psi against the mode-consistent product
    e_node = 0.0
    for p, ref in zip(table["pressure_Pa"], psi_solver_ref):
        e_node = max(e_node, abs(gribs.thermo.psi_and_deriv(float(p), c)[0] / ref - 1.0))
    chk("Psi = R_eff*T0 (table + interpolator)",
        e <= 1.0e-12 and e_node <= 1.0e-10,
        f"stored-vs-product rel. err = {e:.2e}, node interp err = {e_node:.2e}")

    # T11 analytic dPsi/dp against a centred difference
    e = 0.0
    for p in np.geomspace(max(c.p_fit_min * 1.05, 1.1e5), c.p_fit_max * 0.95, 25):
        hstep = p * 1e-5
        fd = (gribs.thermo.psi_and_deriv(p + hstep, c)[0] - gribs.thermo.psi_and_deriv(p - hstep, c)[0]) / (2 * hstep)
        dp = gribs.thermo.psi_and_deriv(p, c)[1]
        if dp == 0.0:
            e = max(e, abs(fd) / max(gribs.thermo.psi_and_deriv(p, c)[0], 1e-300) * p)
        else:
            e = max(e, abs(fd / dp - 1.0))
    chk("analytic vs numerical dPsi/dp", e < tol_deriv,
        f"max rel. err = {e:.2e} (tol {tol_deriv:.0e})")

    # T12 pressure-equation denominator D(p) = 1 - p*Psi'/Psi > 0 everywhere
    try:
        dmin = float(backend.property_validation[
            "min_pressure_factor_1_minus_pPsi_over_Psi"])
        chk("pressure factor 1 - p*Psi'/Psi > 0", dmin > 0.0,
            f"min D(p) = {dmin:.6f}")
    except Exception as exc:                                  # pragma: no cover
        chk("pressure factor 1 - p*Psi'/Psi > 0", False, str(exc))

    # T13 single-phase limit: with Yg = 1 everywhere the two-phase model B must
    # reduce EXACTLY to the legacy Theta model (regression test).
    # T14 nozzle decomposition identities (with a condensed phase present).
    # T15 compatibility: no condensed phase -> homogeneous and legacy runs agree.
    # Cleanup P4 v2 (A-1): every comparison builds its backends with an
    # explicit two-phase mode, independent of the mode of the configuration
    # the self-tests were started with.
    def _cfg_mode(mode):
        cc = copy.copy(c)                      # keeps reactants_spec_list
        cc.two_phase_mode = mode
        return cc
    synth = _SyntheticTestBackend(_cfg_mode(TWO_PHASE_HOMOGENEOUS),
                                  condensed_fraction=0.0)
    try:
        with gribs.thermo._ThermoSwap(synth):
            e_psi, e_dpsi, e_rhs = 0.0, 0.0, 0.0
            c_hom = Config(**{**asdict(c), "two_phase_mode": TWO_PHASE_HOMOGENEOUS})
            c_leg = Config(**{**asdict(c), "two_phase_mode": TWO_PHASE_SINGLE_LEGACY})
            # Cleanup P4 (A-1): the legacy side is evaluated with a backend
            # built in legacy mode, so the backend's own legacy branch is
            # exercised (the CEA/legacy recursion of A-6 3c24fbc went
            # undetected because both sides used the homogeneous backend).
            synth_leg = _SyntheticTestBackend(_cfg_mode(TWO_PHASE_SINGLE_LEGACY),
                                              condensed_fraction=0.0)
            for p in np.geomspace(synth.pmin * 1.01, synth.pmax * 0.99, 30):
                psi_h, dpsi_h = gribs.thermo.psi_and_deriv(float(p), c_hom)
                with gribs.thermo._ThermoSwap(synth_leg):
                    th_l, dth_l = theta_and_deriv(float(p), c_leg)
                e_psi = max(e_psi, abs(psi_h / th_l - 1.0))
                if dth_l != 0.0:
                    e_dpsi = max(e_dpsi, abs(dpsi_h / dth_l - 1.0))
                else:
                    e_dpsi = max(e_dpsi, abs(dpsi_h - dth_l))
                y0 = np.array([
                    float(p),
                    0.3 * web_thickness(c_hom),
                    c.R_t0,
                    1e-5,
                    1e-4,
                    2e-5,
                    0.0,
                ])
                dy_h = make_rhs(c_hom, Ae, True)(0.0, y0)
                with gribs.thermo._ThermoSwap(synth_leg):
                    dy_l = make_rhs(c_leg, Ae, True)(0.0, y0)
                e_rhs = max(e_rhs, float(np.max(np.abs(dy_h - dy_l)
                                                / np.maximum(np.abs(dy_l), 1e-300))))
            chk("single-phase limit: model B == Theta model",
                max(e_psi, e_dpsi, e_rhs) < SINGLE_PHASE_LIMIT_TOLERANCE,
                f"Psi err = {e_psi:.2e}, dPsi err = {e_dpsi:.2e}, "
                f"rhs err = {e_rhs:.2e} (tol {SINGLE_PHASE_LIMIT_TOLERANCE:.0e})")

            # T14 nozzle flow decomposition: mdot_total = mdot_gas + mdot_condensed,
            # mdot_gas = Yg*mdot_total, mdot_condensed = Yc*mdot_total.
            synth_c = _SyntheticTestBackend(_cfg_mode(TWO_PHASE_HOMOGENEOUS),
                                            condensed_fraction=0.3)
            with gribs.thermo._ThermoSwap(synth_c):
                e_dec = 0.0
                c_hom_c = Config(**{**asdict(c), "two_phase_mode": TWO_PHASE_HOMOGENEOUS})
                for p in np.geomspace(c.p_a * 1.05, c.p_fit_max, 30):
                    nz = gribs.nozzle.nozzle_state(float(p), At0, Ae, c_hom_c)
                    if nz["mdot_total"] <= 0.0:
                        continue
                    Yg_p, Yc_p = gribs.thermo.phase_fractions(float(p), c_hom_c)
                    e_dec = max(e_dec,
                                abs(nz["mdot_gas"] + nz["mdot_condensed"]
                                    - nz["mdot_total"]) / nz["mdot_total"],
                                abs(nz["mdot_gas"] / nz["mdot_total"] - Yg_p),
                                abs(nz["mdot_condensed"] / nz["mdot_total"] - Yc_p),
                                abs(nz["mdot_total"] - nz["mdot_total_reference"])
                                / nz["mdot_total_reference"])
                chk("nozzle decomposition mdot_total = gas + condensed",
                    e_dec < 1.0e-12, f"max rel. err = {e_dec:.2e}")

            # T15 compatibility regression: with a condensed-FREE table the
            # homogeneous-equilibrium run must match the legacy run (short motor).
            fast = {**asdict(c), "a_burn": 10.0 * c.a_burn, "blowdown": False}
            c_hom_f = Config(**{**fast, "two_phase_mode": TWO_PHASE_HOMOGENEOUS})
            c_leg_f = Config(**{**fast, "two_phase_mode": TWO_PHASE_SINGLE_LEGACY})
            res_h = gribs.simulation.integrator.run_model(c_hom_f)
            hist_h = gribs.simulation.postprocess.sample(res_h, c_hom_f)
            s_h = gribs.simulation.postprocess.summarize(res_h, hist_h, c_hom_f)
            with gribs.thermo._ThermoSwap(synth_leg):
                res_l = gribs.simulation.integrator.run_model(c_leg_f)
                hist_l = gribs.simulation.postprocess.sample(res_l, c_leg_f)
                s_l = gribs.simulation.postprocess.summarize(res_l, hist_l, c_leg_f)
            e_run = max(abs(s_h["p_max_Pa"] / s_l["p_max_Pa"] - 1.0),
                        abs(s_h["burn_time_s"] / s_l["burn_time_s"] - 1.0),
                        abs(s_h["impulse_total_Ns"] / s_l["impulse_total_Ns"] - 1.0),
                        abs(s_h["Isp_s"] / s_l["Isp_s"] - 1.0))
            t_cmp = np.linspace(0.0, 0.98 * min(s_h["burn_time_s"], s_l["burn_time_s"]), 200)
            p_h = res_h["sol_burn"].sol(t_cmp)[IP]
            p_l = res_l["sol_burn"].sol(t_cmp)[IP]
            e_traj = float(np.max(np.abs(p_h / p_l - 1.0)))
            e_run = max(e_run, e_traj)
            chk("no-condensed compatibility (hom. vs legacy run)", e_run < 1.0e-9,
                f"summary/trajectory rel. err = {e_run:.2e} (tol 1e-9)")
    except Exception as exc:                                  # pragma: no cover
        chk("two-phase regression tests", False, f"{type(exc).__name__}: {exc}")

    name_contract = self_test_name_contract(c, result_names)
    ok_all &= bool(name_contract["passed"])

    if verbose:
        print(f"Self-tests  (thermochemistry backend: {backend.name}, "
              f"two-phase mode: {c.two_phase_mode}"
              f"{'' if two_phase_active else ' (phase split inactive)'})")
        print("\n".join(out))
        if not name_contract["passed"]:
            print("SELF-TEST CONTRACT ERROR")
            print(
                "  expected count: "
                f"{name_contract['expected_count']}, "
                "actual count: "
                f"{name_contract['actual_count']}"
            )
            if name_contract["missing_names"]:
                print(
                    "  missing names: "
                    + ", ".join(name_contract["missing_names"])
                )
            if name_contract["unexpected_names"]:
                print(
                    "  unexpected names: "
                    + ", ".join(name_contract["unexpected_names"])
                )
            if name_contract["duplicate_names"]:
                print(
                    "  duplicate names: "
                    + ", ".join(name_contract["duplicate_names"])
                )
        print(f"  -> {'all checks passed' if ok_all else 'FAILURES DETECTED'}\n")
    return ok_all


def self_tests_post_run(
    res: dict,
    h: dict,
    c: Config,
    verbose=True,
    report: Optional[dict] = None,
) -> bool:
    """Display and return the structured post-run consistency verdict."""
    report = (
        post_run_consistency_report(res, h, rtol=c.rtol)
        if report is None
        else report
    )

    final = report["final_total_mass"]
    eos = report["eos_history"]

    out = [
        (
            "  [{status}] {name:<46s} "
            "mt0 + m_gen = {lhs:.9e} kg vs "
            "m_out + mt = {rhs:.9e} kg "
            "(rel. err = {error:.2e})"
        ).format(
            status="PASS" if final["passed"] else "FAIL",
            name="total mass conservation (final)",
            lhs=final["lhs_initial_plus_generated_kg"],
            rhs=final["rhs_discharged_plus_remaining_kg"],
            error=final["relative_error"],
        ),
        (
            "  [{status}] {name:<46s} "
            "max rel. err = {error:.2e} "
            "(fixed reference = {reference:.9e} kg, "
            "tol {tolerance:.0e})"
        ).format(
            status="PASS" if eos["passed"] else "FAIL",
            name="EOS mass vs balance mass (history)",
            error=eos["invariant_relative_error"],
            reference=eos["reference_mass_kg"],
            tolerance=report["tolerance"],
        ),
        (
            "  [INFO] {name:<46s} "
            "max rel. err = {error:.2e}"
        ).format(
            name="legacy shrinking-mass metric (deprecated)",
            error=eos[
                "legacy_instantaneous_mass_relative_error"
            ],
        ),
    ]

    if verbose:
        print("Post-run conservation tests")
        print("\n".join(out))
        print(
            f"  -> {'all checks passed' if report['passed'] else 'FAILURES DETECTED'}\n"
        )

    return bool(report["passed"])
