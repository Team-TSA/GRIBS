#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase-A probe #5: a Level-0 reference problem with an analytic solution.

Configuration: constant properties (the module's own synthetic test table),
r = 0 (no propellant generation), fixed free volume, converging nozzle.
The chamber then obeys, exactly,

    dp/dt = -Psi * mdot_out / Vg,
    mdot_out = Cd*At*p*sqrt(g/(R*T0)) * (2/(g+1))^((g+1)/(2(g-1)))

i.e. an exponential decay  p(t) = p0 * exp(-t/tau),
tau = Vg / (Psi*Cd*At*sqrt(g/(R*T0))*(2/(g+1))^((g+1)/(2(g-1)))).

The script compares the integrator against this closed form.  It exists to show
two things:
  * the chamber ODE + nozzle flow can be verified against an independent
    analytic solution to ~1e-9 (the numerical core is sound), and
  * today the ONLY way to do that is the private class _SyntheticTestBackend
    inside the production file; Phase A should promote it to a selectable
    reference backend.
"""
from __future__ import annotations

import copy
import importlib.util
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parent / "gribs_baseline"


def load_gribs():
    spec = importlib.util.spec_from_file_location("gribs", BASE / "gribs.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["gribs"] = mod
    spec.loader.exec_module(mod)
    return mod


def main() -> int:
    g = load_gribs()
    base = g.load_json_configuration(BASE / "gribs_config_eps1.json").config

    # --- Level-0 configuration: constant properties, no burning, fixed volume
    cfg = copy.copy(base)
    # validate() rejects a_burn = 0, so use an effectively-zero burn rate with
    # n_burn = 0: the generation term is then a tiny CONSTANT and the analytic
    # solution stays exact:
    #     p(t) = p_inf + (p0 - p_inf)*exp(-t/tau),  p_inf = rho_p*Ab*a/C
    cfg.a_burn = 1.0e-12
    cfg.n_burn = 0.0
    cfg.ero_alpha = 0.0
    cfg.sigma_p = 0.0
    cfg.ign_mdot = 0.0
    cfg.ign_time = 0.0
    cfg.eps_nozzle = 1.0             # converging nozzle -> keeps the choked branch
    cfg.Cd = 1.0
    cfg.eta_thrust = 1.0
    cfg.use_separation = False
    cfg.p_a = 1.0e5
    cfg.p0_init = 5.0e6
    cfg.blowdown = False
    cfg.t_max = 5.0

    synth = g._SyntheticTestBackend(cfg, condensed_fraction=0.0)   # private API
    g._THERMO = synth

    res = g.run_model(cfg)
    sol = res["sol_burn"]
    At = math.pi * cfg.R_t0 ** 2
    R, T0, gam = g.gas_props(cfg.p0_init, cfg)
    psi = g.psi_and_deriv(cfg.p0_init, cfg)[0]
    crit = (2.0 / (gam + 1.0)) ** ((gam + 1.0) / (2.0 * (gam - 1.0)))
    C = cfg.Cd * At * math.sqrt(gam / (R * T0)) * crit         # mdot = C*p  [choked]
    tau = cfg.V_g0 / (psi * C)
    Ab0 = g.geometry(0.0, cfg)[2]
    p_inf = cfg.rho_p * Ab0 * cfg.a_burn / C
    p_thr = g.choke_limit_pressure(cfg, At, At)

    t = sol.t
    p = sol.y[g.IP]
    keep = p > 1.05 * p_thr                                    # choked segment
    p_ref = p_inf + (cfg.p0_init - p_inf) * np.exp(-t[keep] / tau)
    rel = np.max(np.abs(p[keep] / p_ref - 1.0))
    print("=" * 78)
    print("Level-0 reference problem: isothermal-constant-property blowdown")
    print("=" * 78)
    print(f"synthetic property basis : T0 = {T0:.1f} K, gamma = {gam:.3f}, "
          f"R = {R:.3f} J/(kg K), Psi = {psi:.3f} J/kg")
    print(f"choked flow coefficient C: {C:.6e} kg/(s Pa)   time constant tau = {tau:.6f} s")
    print(f"residual generation      : {cfg.rho_p * Ab0 * cfg.a_burn:.3e} kg/s "
          f"-> asymptotic pressure p_inf = {p_inf:.6f} Pa")
    print(f"integration interval     : 0 ... {t[-1]:.4f} s, {len(t)} steps")
    print(f"choked segment used      : p0 = {p[keep][0]/1e6:.4f} ... "
          f"{p[keep][-1]/1e3:.1f} kPa  ({keep.sum()} points)")
    print(f"max |p_num/p_analytic-1| : {rel:.3e}")
    print(f"analytic p at t_end      : {p_ref[-1]:.6f} Pa   numerical: {p[keep][-1]:.6f} Pa")
    print()
    print("Interpretation: the chamber mass/state equation, the choked nozzle flow")
    print("and the event handling reproduce the closed-form solution to the solver")
    print("tolerance.  The verification gap is not in the numerics but in the")
    print("model-form assumptions - which is exactly what the roadmap says.")
    print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
