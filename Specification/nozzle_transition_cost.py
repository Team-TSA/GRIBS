#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase-A probe #3: what the C-D nozzle switch costs in a real run.

The pre-run self-test T4 aborts every configuration with expansion_ratio > 1,
so the only way to see the effect is to call the solver directly, bypassing the
self-tests (exactly what a user cannot do today).

Reported, for the same motor with Ae/At = 4:
  * time and impulse spent in each nozzle regime,
  * the impulse accumulated in the choked/subsonic transition band,
  * the same band evaluated with two alternative branch treatments:
      - the model's own subsonic-branch formula (continuous at the switch),
      - the classical normal-shock solution (independent reference).
"""
from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
BASE = HERE.parent / "gribs_baseline"


def load_gribs():
    spec = importlib.util.spec_from_file_location("gribs", BASE / "gribs.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["gribs"] = mod
    spec.loader.exec_module(mod)
    return mod


def trapz(y, x):
    return float(np.trapezoid(y, x)) if hasattr(np, "trapezoid") else float(np.trapz(y, x))


def masked_trapz(y, x, mask):
    """Trapezoidal integral over the *contiguous* pieces of a boolean mask.

    np.trapz on a masked slice would bridge the gaps between pieces and invent
    impulse (a regime can be entered twice in one run)."""
    m = np.asarray(mask, dtype=bool)
    pair = m[:-1] & m[1:]
    return float(np.sum(0.5 * (y[:-1] + y[1:])[pair] * np.diff(x)[pair]))


def main() -> int:
    gribs = load_gribs()
    loaded = gribs.load_json_configuration(BASE / "gribs_config.json")   # Ae/At = 4
    cfg = loaded.config
    gribs.initialize_thermochemistry(cfg, BASE / "results" / "thermo_cache")
    c = cfg
    At = math.pi * c.R_t0 ** 2
    Ae = c.eps_nozzle * At
    pa = c.p_a
    p_thr = gribs.choke_limit_pressure(c, At, Ae)

    res = gribs.run_model(c)                       # <-- self-tests bypassed
    tg, ys = [], []
    for sol in (res["sol_burn"], res["sol_blow"]):
        if sol is None:
            continue
        dt = sol.t[-1] - sol.t[0]
        n = max(int(2000.0 * dt) + 1, 2001)          # equal time resolution
        t_seg = sol.t[0] + dt * np.linspace(0.0, 1.0, n)
        tg.append(t_seg); ys.append(sol.sol(t_seg))
    tg = np.concatenate(tg); ys = np.concatenate(ys, axis=1)
    order = np.argsort(tg)
    tg = tg[order]; ys = ys[:, order]
    p0 = np.maximum(ys[gribs.IP], 1.0)
    F, mdot, regime, pe = [], [], [], []
    for p in p0:
        nz = gribs.nozzle_state(float(p), At, Ae, c)
        F.append(nz["F"]); mdot.append(nz["mdot_total"])
        regime.append(nz["regime"]); pe.append(nz["pe"])
    F = np.array(F); mdot = np.array(mdot); regime = np.array(regime); pe = np.array(pe)
    I_tot = trapz(F, tg)
    print("=" * 78)
    print("C-D nozzle (Ae/At = 4), burning phase + blowdown, self-tests bypassed")
    print("=" * 78)
    print(f"choking onset pressure           : {p_thr:.1f} Pa")
    print(f"total computed time               : {tg[-1]:.6f} s")
    print(f"regimes visited                  : {sorted(set(regime))}")
    dt_all = np.diff(np.concatenate([[tg[0]], tg]))      # time weight of each sample
    for r in sorted(set(regime)):
        m = regime == r
        print(f"   {r:<10} {dt_all[m].sum()*100.0/tg[-1]:6.2f} % of the time, "
              f"{masked_trapz(F, tg, m):10.4f} N.s")
    print(f"total computed impulse           : {I_tot:.4f} N.s (burning + blowdown)")

    # ---- transition band: from the choke switch up to where the model's own
    # separation criterion takes over (p_e,full = sep_ratio*p_a)
    Me = gribs.mach_from_area(c.eps_nozzle, 1.2, True)          # representative
    R, T0, g = gribs.gas_props(5.0e6, c)
    Me = gribs.mach_from_area(c.eps_nozzle, g, True)
    pe_over_p0 = (1 + 0.5 * (g - 1) * Me ** 2) ** (-g / (g - 1))
    p_sep = c.sep_ratio * pa / pe_over_p0
    print(f"\nmodel's separation onset (p_e = 0.4 p_a) : {p_sep:,.0f} Pa")
    m_band = (p0 > p_thr) & (p0 < p_sep)
    print(f"time inside the transition band          : "
          f"{dt_all[m_band].sum()*1000.0:.2f} ms "
          f"({dt_all[m_band].sum()*100.0/tg[-1]:.2f} % of the run)")
    I_band = masked_trapz(F, tg, m_band)
    print(f"impulse inside the band (model)          : {I_band:.4f} N.s "
          f"({100*I_band/I_tot:.3f} % of the total computed impulse)")

    # alternative treatments inside the band
    F_sub = np.zeros_like(p0); F_sh = np.full_like(p0, np.nan)
    for i, p in enumerate(p0):
        if not m_band[i]:
            continue
        R_, T0_, g_ = gribs.gas_props(float(p), c)
        sq = math.sqrt(g_ / (R_ * T0_))
        Me_s = math.sqrt(max(2.0 / (g_ - 1.0) * ((p / pa) ** ((g_ - 1.0) / g_) - 1.0), 0.0))
        md_ref = c.Cd * Ae * p * sq * Me_s \
            * (1.0 + 0.5 * (g_ - 1.0) * Me_s ** 2) ** (-(g_ + 1.0) / (2.0 * (g_ - 1.0)))
        Yg, _ = gribs.phase_fractions(float(p), c)
        md_tot = md_ref / Yg
        Te = T0_ / (1.0 + 0.5 * (g_ - 1.0) * Me_s ** 2)
        F_sub[i] = c.eta_thrust * md_tot * Me_s * math.sqrt(g_ * R_ * Te)
        F_sh[i] = shock_thrust(g_, R_, T0_, md_tot, float(p), c.eps_nozzle, pa)
    print(f"impulse inside the band (model's own subsonic branch): "
          f"{masked_trapz(F_sub, tg, m_band):.4f} N.s")
    ok = np.isfinite(F_sh)
    if ok.any():
        idx = m_band & ok
        print(f"impulse inside the band (normal-shock reference)     : "
              f"{masked_trapz(np.nan_to_num(F_sh), tg, idx):.4f} N.s  "
              f"({idx.sum()} of {m_band.sum()} samples have a shock solution)")
    print()
    return 0


def shock_thrust(g, R, T0, mdot_total, p0, eps, pa):
    Me_sub = brentq(lambda M: (1.0 / M) * ((2.0 / (g + 1.0))
                    * (1.0 + 0.5 * (g - 1.0) * M * M)) ** ((g + 1.0) / (2.0 * (g - 1.0)))
                    - eps, 1e-7, 1.0, xtol=1e-15)
    p_ratio_exit = (1.0 + 0.5 * (g - 1.0) * Me_sub ** 2) ** (-g / (g - 1.0))

    def pe_over_p0(M1):
        pt2_pt1 = (((g + 1.0) * M1 * M1) / ((g - 1.0) * M1 * M1 + 2.0)) ** (g / (g - 1.0)) \
            * ((g + 1.0) / (2.0 * g * M1 * M1 - (g - 1.0))) ** (1.0 / (g - 1.0))
        return pt2_pt1 * p_ratio_exit
    try:
        M1 = brentq(lambda m: p0 * pe_over_p0(m) - pa, 1.0 + 1e-9, 200.0, xtol=1e-12)
    except ValueError:
        return float("nan")
    A_sh = (1.0 / M1) * ((2.0 / (g + 1.0)) * (1.0 + 0.5 * (g - 1.0) * M1 * M1)) ** \
        ((g + 1.0) / (2.0 * (g - 1.0)))
    if not (1.0 < A_sh < eps):
        return float("nan")
    Te = T0 / (1.0 + 0.5 * (g - 1.0) * Me_sub ** 2)
    return mdot_total * Me_sub * math.sqrt(g * R * Te)


if __name__ == "__main__":
    raise SystemExit(main())
