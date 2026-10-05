#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase-A probe #4: numerical-convergence sweep prototype.

This is the kind of tool the roadmap asks for in step 2 ('数値収束試験'): it
perturbs one numerical setting at a time, re-runs the same motor and reports the
converged quantities plus the relative change against the finest run.

It also demonstrates two things about the current code that matter for Phase A:
  * the solver + thermochemistry can be driven from an external script only
    through module-level globals (``_THERMO``), so the sweep must live in the
    same process as the engine,
  * the CEA table is a *numerical* setting as well: changing pressure_points
    changes the interpolated properties and therefore the "converged" answer.

Usage:  python3 convergence_sweep.py [config.json]
Writes: convergence_sweep.csv  and  prints a Markdown table.
"""
from __future__ import annotations

import copy
import importlib.util
import json
import math
import sys
import time
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


def metrics(g, res, c, n=2001):
    """Converged quantities taken straight from the ODE solution (no sampling of
    the plotting/interpolation layer is involved)."""
    t_all, p_all, F_all = [], [], []
    for sol in (res["sol_burn"], res["sol_blow"]):
        if sol is None:
            continue
        tt = sol.t[0] + (sol.t[-1] - sol.t[0]) * np.linspace(0.0, 1.0, n)
        ys = sol.sol(tt)
        pp = np.maximum(ys[g.IP], 1.0)
        At = math.pi * np.maximum(ys[g.IRT], 1e-9) ** 2
        FF = np.array([g.nozzle_state(float(p), float(a), res["Ae"], c)["F"]
                       for p, a in zip(pp, At)])
        t_all.append(tt); p_all.append(pp); F_all.append(FF)
    t = np.concatenate(t_all); p = np.concatenate(p_all); F = np.concatenate(F_all)
    tr = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    sol_b = res["sol_burn"]
    tb = sol_b.t
    yb = sol_b.y
    p_b = np.maximum(yb[g.IP], 1.0)
    mt0 = res["m_total0"]
    bal = mt0 + yb[g.IMG] - yb[g.IMO]
    psi = np.array([g.psi_and_deriv(float(pp), c)[0] for pp in p_b])
    Vg = np.array([g.geometry(float(x), c)[3] for x in yb[g.IX]])
    eos = p_b * Vg / psi
    mass_err = float(np.max(np.abs(bal - eos) / np.maximum(np.abs(eos), 1e-15)))
    return dict(
        p_max=float(np.max(p)), t_p_max=float(t[int(np.argmax(p))]),
        F_max=float(np.max(F)),
        burn_time=float(tb[-1]),
        impulse_total=float(tr(F, t)),
        m_out_total=float(yb[g.IMO][-1]),
        p_final=float(yb[g.IP][-1]),
        max_mass_error=mass_err,
        n_steps=int(len(tb) + (len(res["sol_blow"].t) if res["sol_blow"] else 0)),
    )


def main(config_name: str) -> int:
    g = load_gribs()
    cfg_path = BASE / config_name
    base_cfg = g.load_json_configuration(cfg_path).config
    cache_dir = BASE / "results" / "thermo_cache"

    cases = []
    for rtol in (1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11):
        cases.append(("rtol", f"{rtol:.0e}", dict(rtol=rtol)))
    for atol in (1e-1, 1e-2, 1e-3, 1e-4):
        cases.append(("atol_p", f"{atol:.0e}", dict(atol_p=atol)))
    for ms in (0.2, 0.05, 0.02, 0.005):
        cases.append(("max_step_burn", f"{ms:g}", dict(max_step_burn=ms)))
    for np_ in (21, 41, 81, 161):
        cases.append(("cea_pressure_points", f"{np_}", dict(cea_pressure_points=np_)))

    rows = []
    print(f"{'setting':>20} {'value':>6} {'p_max [MPa]':>12} {'t_bo [s]':>10} "
          f"{'F_max [N]':>10} {'I [N.s]':>10} {'m_out [kg]':>11} "
          f"{'eps_m':>10} {'steps':>7} {'wall [s]':>8}")
    last_key, last_val = None, None
    for key, val, over in cases:
        cfg = copy.copy(base_cfg)          # keeps reactants_spec_list etc.
        for k, v in over.items():
            setattr(cfg, k, v)
        if key != last_key:
            # rebuild the backend only when the thermochemistry options change
            g.initialize_thermochemistry(cfg, cache_dir)
        t0 = time.time()
        try:
            res = g.run_model(cfg)
            m = metrics(g, res, cfg)
        except Exception as exc:                       # keep the sweep running
            print(f"{key:>20} {val:>6}   FAILED: {type(exc).__name__}: {exc}")
            rows.append(dict(setting=key, value=val, error=f"{type(exc).__name__}: {exc}"))
            last_key, last_val = key, val
            continue
        wall = time.time() - t0
        rows.append(dict(setting=key, value=val, wall_s=round(wall, 2), **m))
        print(f"{key:>20} {val:>6} {m['p_max']/1e6:12.6f} {m['burn_time']:10.6f} "
              f"{m['F_max']:10.4f} {m['impulse_total']:10.4f} {m['m_out_total']:11.8f} "
              f"{m['max_mass_error']:10.2e} {m['n_steps']:7d} {wall:8.2f}")
        last_key, last_val = key, val

    # relative change against the last (finest) case of each setting
    print("\nrelative change of each setting against its last (finest) case")
    from collections import defaultdict
    per = defaultdict(list)
    for r in rows:
        if "error" not in r:
            per[r["setting"]].append(r)
    for key, rs in per.items():
        ref = rs[-1]
        for r in rs[:-1]:
            d = {k: f"{(r[k]/ref[k]-1.0)*100:+.4f}%" for k in
                 ("p_max", "burn_time", "F_max", "impulse_total", "m_out_total")}
            print(f"  {key:>20} {r['value']:>6} vs {ref['value']:>6} : "
                  f"p_max {d['p_max']:>10}  t_bo {d['burn_time']:>10}  "
                  f"F_max {d['F_max']:>10}  I {d['impulse_total']:>10}  "
                  f"m_out {d['m_out_total']:>10}")
    out = HERE / "convergence_sweep.json"
    out.write_text(json.dumps(rows, indent=2), encoding="utf-8")
    print(f"\nwritten: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1] if len(sys.argv) > 1 else "gribs_config_wide.json"))
