#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase-A probe #2: the C-D nozzle regime switch of GRIBS v0.5.0-alpha.

Two questions, both answered with the code itself (no re-implementation of
GRIBS physics):

  Q1  How large is the discontinuity of the modelled thrust at the
      choked/subsonic switch, as a function of the expansion ratio?
      (This is what self-test T4 measures; it aborts the run when eps > 1.)

  Q2  Between the choking onset and the "design" pressure (where the fully
      expanded supersonic exit pressure equals the ambient pressure) the real
      nozzle flow is choked with a normal shock inside the diverging section
      and p_e = p_a.  GRIBS instead uses the fully expanded supersonic solution.
      How large is the resulting thrust error?

The shock solution is computed here independently (classical 1-D isentropic +
Rankine-Hugoniot relations) and is only used as a reference; it is never fed
back into GRIBS.
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


def shock_thrust(g, R, T0, mdot_total, p0, eps, pa):
    """Thrust when the flow is choked with a normal shock inside the diverging
    part and the exit is subsonic at p_e = p_a (independent reference solution).

    Returns (F, M_exit, shock_area_ratio) or None when no such solution exists
    (choking onset region or p0 below the design point of the shock-free case).
    """
    Me_sub = brentq(lambda M: (1.0 / M) * ((2.0 / (g + 1.0))
                    * (1.0 + 0.5 * (g - 1.0) * M * M)) ** ((g + 1.0) / (2.0 * (g - 1.0)))
                    - eps, 1e-7, 1.0, xtol=1e-15)
    # p_e/p_t2 for the subsonic exit
    p_ratio_exit = (1.0 + 0.5 * (g - 1.0) * Me_sub ** 2) ** (-g / (g - 1.0))

    def pe_over_p0(M1):
        # total-pressure ratio across the standing normal shock at M1
        pt2_pt1 = (((g + 1.0) * M1 * M1) / ((g - 1.0) * M1 * M1 + 2.0)) ** (g / (g - 1.0)) \
            * ((g + 1.0) / (2.0 * g * M1 * M1 - (g - 1.0))) ** (1.0 / (g - 1.0))
        return pt2_pt1 * p_ratio_exit

    f = lambda M1: p0 * pe_over_p0(M1) - pa
    # search over upstream supersonic Mach numbers whose area ratio is <= Ae
    M_lo, M_hi = 1.0 + 1e-9, 200.0
    try:
        if f(M_hi) > 0.0:            # shock would have to sit outside the nozzle
            return None
        M1 = brentq(f, M_lo, M_hi, xtol=1e-12, maxiter=200)
    except ValueError:
        return None
    A_sh = (1.0 / M1) * ((2.0 / (g + 1.0)) * (1.0 + 0.5 * (g - 1.0) * M1 * M1)) ** \
        ((g + 1.0) / (2.0 * (g - 1.0)))
    if 1.0 < A_sh < eps:
        Te = T0 / (1.0 + 0.5 * (g - 1.0) * Me_sub ** 2)
        ve = Me_sub * math.sqrt(g * R * Te)
        return mdot_total * ve, Me_sub, A_sh
    return None


def main() -> int:
    gribs = load_gribs()
    cfg = gribs.load_json_configuration(BASE / "gribs_config_eps1.json").config
    gribs.initialize_thermochemistry(cfg, BASE / "results" / "thermo_cache")

    At = math.pi * cfg.R_t0 ** 2
    pa = cfg.p_a
    print("=" * 78)
    print("Q1  thrust discontinuity at the choked/subsonic switch (GRIBS model)")
    print("=" * 78)
    print(f"{'Ae/At':>6} {'p_switch [Pa]':>14} {'mdot jump':>11} {'F_sub [N]':>10} "
          f"{'F_choked [N]':>13} {'F jump':>9}")
    rows = []
    for eps in (1.0, 1.5, 2.0, 4.0, 6.0):
        Ae = eps * At
        p_thr = gribs.choke_limit_pressure(cfg, At, Ae)
        if p_thr is None:
            print(f"{eps:6.2f}   choke limit could not be bracketed")
            continue
        lo = gribs.nozzle_state(p_thr * (1 - 1e-9), At, Ae, cfg)
        hi = gribs.nozzle_state(p_thr * (1 + 1e-9), At, Ae, cfg)
        jm = abs(lo["mdot_total"] / hi["mdot_total"] - 1.0)
        jf = abs(lo["F"] / hi["F"] - 1.0)
        rows.append((eps, p_thr, jm, lo["F"], hi["F"], jf))
        print(f"{eps:6.2f} {p_thr:14.1f} {jm:11.2e} {lo['F']:10.3f} {hi['F']:13.3f} "
              f"{jf:9.3f}")

    print()
    print("=" * 78)
    print("Q2  off-design error of the C-D nozzle model (complete-entrainment basis)")
    print("=" * 78)
    print("  'shock solution' = choked flow, normal shock in the diverging part,")
    print("  subsonic exit at p_e = p_a  (the classical 1-D solution GRIBS omits).")
    for eps in (2.0, 4.0, 6.0):
        Ae = eps * At
        # pressure at which the fully expanded supersonic exit pressure equals p_a
        R, T0, g = gribs.gas_props(1.0e6, cfg)
        Me = gribs.mach_from_area(eps, g, True)
        ratio = (1 + 0.5 * (g - 1) * Me * Me) ** (-g / (g - 1))   # p_e/p_0
        p_design = pa / ratio                                     # p_e = p_a
        p_thr = gribs.choke_limit_pressure(cfg, At, Ae)
        print(f"\n  Ae/At = {eps:g}:  choking onset p = {p_thr:,.0f} Pa, "
              f"design (p_e = p_a) p = {p_design:,.0f} Pa")
        print(f"    {'p0 [Pa]':>12} {'F model [N]':>12} {'F shock [N]':>12} "
              f"{'error':>9} {'shock at A/At':>14}")
        for p0 in (0.3e6, 0.6e6, 1.0e6, 2.0e6, 4.0e6):
            if p0 > 1.0001 * p_design:
                continue
            R, T0, g = gribs.gas_props(p0, cfg)
            nz = gribs.nozzle_state(p0, At, Ae, cfg)
            ref = shock_thrust(g, R, T0, nz["mdot_total"], p0, eps, pa)
            if ref is None:
                print(f"    {p0:12.0f} {nz['F']:12.3f} {'(no shock solution)':>12}")
                continue
            Fs, Me_s, A_sh = ref
            print(f"    {p0:12.0f} {nz['F']:12.3f} {Fs:12.3f} "
                  f"{(nz['F']/Fs - 1.0):+9.1%} {A_sh:14.3f}")
    print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
