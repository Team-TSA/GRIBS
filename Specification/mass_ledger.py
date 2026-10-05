#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase-A probe #1: itemised mass ledger rebuilt from a GRIBS history CSV.

The point of this script is to show what the *existing* outputs already allow
and what they do NOT allow, using only the files GRIBS v0.5.0-alpha writes:

  results/*.csv   (time history, total-product-mass basis)
  results/*.json  (summary, resolved configuration)

Usage:
    python3 mass_ledger.py results/conv_nozzle_history.csv results/conv_nozzle_summary.json
"""
from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import numpy as np


def trapz(y, x):
    return float(np.trapezoid(y, x)) if hasattr(np, "trapezoid") else float(np.trapz(y, x))


def load(csv_path: Path):
    with csv_path.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    cols = {}
    for k in rows[0]:
        vals = [r[k] for r in rows]
        try:
            cols[k] = np.array([float(v) for v in vals])
        except ValueError:
            cols[k] = np.array(vals, dtype=object)
    return cols


def main(csv_path: Path, json_path: Path) -> int:
    h = load(csv_path)
    s = json.loads(json_path.read_text(encoding="utf-8"))
    cfg = s["resolved_configuration"]
    summ = s["results"]
    t = h["t[s]"]
    burn = h["phase"] == "burn"

    # ------------------------------------------------------------------ inputs
    rho_p = cfg["rho_p"]
    R_i0, R_p, L_p0, n_end = cfg["R_i0"], cfg["R_p"], cfg["L_p0"], cfg["n_end"]
    V_p0 = math.pi * (R_p**2 - R_i0**2) * L_p0
    m_prop = rho_p * V_p0
    ign_total = cfg["ign_mdot"] * cfg["ign_time"]

    # ------------------------------------------------------ ODE state at t_end
    mt0 = float(h["m_total_eos[kg]"][0])            # initial chamber inventory
    m_gen = float(h["m_gen[kg]"][-1])               # integrator state (total products)
    m_out = float(h["m_out[kg]"][-1])               # integrator state (discharge)
    mt_end = float(h["m_total_eos[kg]"][-1])        # chamber content from EOS p*Vg/Psi

    # geometry-based residual propellant (exact, from the resolved geometry)
    x_end = float(h["x[m]"][-1])
    xw = min(R_p - R_i0, (L_p0 / n_end) if n_end > 0 else math.inf)
    xe = min(max(x_end, 0.0), xw)
    Ri = min(R_i0 + xe, R_p)
    Lp = max(L_p0 - n_end * xe, 0.0)
    m_prop_remaining = rho_p * math.pi * max(R_p**2 - Ri**2, 0.0) * Lp
    m_prop_burned = m_prop - m_prop_remaining

    # ------------------------------------------------------ discharge quadrature
    md_out_tot = h["mdot_out_total[kg/s]"]
    md_out_gas = h["mdot_out_gas[kg/s]"]
    md_out_cnd = h["mdot_out_condensed[kg/s]"]
    out_tot_q = trapz(md_out_tot, t)
    out_gas_q = trapz(md_out_gas, t)
    out_cnd_q = trapz(md_out_cnd, t)

    lines = []
    add = lines.append
    add("=" * 78)
    add("GRIBS v0.5.0-alpha  --  itemised mass ledger (rebuilt from the CSV)")
    add("=" * 78)
    add(f"configuration          : {s['configuration_source']}")
    add(f"propellant (geometry)  : {m_prop:.9e} kg   (rho_p = {rho_p:.3f} kg/m3)")
    add(f"igniter charge         : {ign_total:.9e} kg   (mdot = {cfg['ign_mdot']} kg/s, "
        f"t = {cfg['ign_time']} s)")
    add("")
    add("[ ledger items ]")
    add(f"  A  initial chamber inventory mt0          : {mt0:.9e} kg")
    add(f"  B  generated products  m_gen (state)      : {m_gen:.9e} kg")
    add(f"     B1  of which propellant-derived        : {m_gen - ign_total:.9e} kg")
    add(f"     B2  of which igniter-derived           : {ign_total:.9e} kg")
    add(f"  C  discharged          m_out (state)      : {m_out:.9e} kg")
    add(f"  D  chamber content     p0*Vg/Psi (EOS)    : {mt_end:.9e} kg")
    add(f"  E  unburned solid propellant (geometry)   : {m_prop_remaining:.9e} kg")
    add(f"  F  burned solid propellant   (geometry)   : {m_prop_burned:.9e} kg")
    add("")
    add("[ closure residuals ]")
    lhs = mt0 + m_gen
    rhs = m_out + mt_end
    add(f"  A + B        = {lhs:.9e} kg")
    add(f"  C + D        = {rhs:.9e} kg")
    add(f"  eps_m (F.1)  = |A+B-C-D| / (A+B)          : {abs(lhs-rhs)/abs(lhs):.3e}")
    e_prop = (m_gen - ign_total) - m_prop_burned
    add(f"  F.2  propellant ledger: (B - igniter) - F : {e_prop:+.3e} kg "
        f"({e_prop/max(m_prop,1e-30):+.3e} relative to m_prop)")
    add(f"  F.3  the number GRIBS itself prints as 'propellant mass balance error'")
    add(f"       is  |B / m_prop - 1| = {abs(m_gen/m_prop - 1.0):.3e}  <-- igniter mass")
    add(f"       ({ign_total:.3e} kg) is folded into m_gen, so F.3 reports the igniter")
    add("       charge as a propellant error.  The itemised ledger does not.")
    add("")
    add("[ what the existing outputs do NOT close: the phase split ]")
    add("  The ODE carries ONE discharge state m_out (gas + condensed together).")
    add("  The gas/condensed split of the discharge exists only as sampled rates,")
    add("  so its accumulated value depends on the post-processing grid:")
    add(f"    integral(mdot_out_total)     = {out_tot_q:.9e} kg   (grid quadrature)")
    add(f"    ODE state m_out                                 = {m_out:.9e} kg")
    add(f"    grid quadrature error                           = "
        f"{abs(out_tot_q - m_out)/m_out:.3e}")
    add(f"    integral(mdot_out_gas)       = {out_gas_q:.9e} kg")
    add(f"    integral(mdot_out_condensed) = {out_cnd_q:.9e} kg")
    add(f"    gas + condensed             = {out_gas_q + out_cnd_q:.9e} kg")
    add(f"    gas/condensed split of the DISCHARGE therefore carries a quadrature")
    add(f"    error evaluated over the same sampled grid.")
    add("")
    add("[ chamber residual split (post-processed, not an ODE state) ]")
    add(f"  m_gas_equilibrium (final) : {float(h['m_gas_equilibrium[kg]'][-1]):.9e} kg")
    add(f"  m_condensed       (final) : {float(h['m_condensed[kg]'][-1]):.9e} kg")
    add(f"  Yg final = {float(h['Y_gas[-]'][-1]):.6f}   Yc final = "
        f"{float(h['Y_condensed[-]'][-1]):.6f}")
    add("")
    add("[ energy ledger ]  : NOT AVAILABLE -- no energy state is integrated and no")
    add("  energy residual is written.  (Only the model-internal consistency of p, V,")
    add("  m and Psi is checked.)")
    text = "\n".join(lines)
    print(text)
    out = Path(csv_path).with_name(Path(csv_path).stem + "_ledger.txt")
    out.write_text(text + "\n", encoding="utf-8")
    print(f"\nwritten: {out}")
    return 0


if __name__ == "__main__":
    a = sys.argv[1:]
    default_dir = Path(__file__).resolve().parent.parent / "gribs_baseline" / "results"
    csv_p = Path(a[0]) if a else default_dir / "conv_nozzle_history.csv"
    js_p = Path(a[1]) if len(a) > 1 else default_dir / "conv_nozzle_summary.json"
    raise SystemExit(main(csv_p, js_p))
