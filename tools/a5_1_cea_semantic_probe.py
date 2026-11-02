"""A-5.1 G1: independent CEA mass-basis semantic probe.

Uses the official ``cea`` package directly (no GRIBS code) to measure, for the
reactant sets of the existing golden configs, what EqSolution.M / n / MW /
density / mass_fractions mean, and where a sonic velocity is available.

Outputs <outdir>/cea_semantics.json and <outdir>/cea_semantics.csv.
Usage: python tools/a5_1_cea_semantic_probe.py <outdir> <config.json>...
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
import cea

PRESSURES_PA = [1.0e5, 1.0e6, 3.0e6, 8.0e6]
RU = float(cea.R)  # J/(kmol K); CEA 3.3.4 reports 8314.51 (same as GRIBS R_UNIVERSAL)


def probe_case(cfg_path):
    doc = json.loads(Path(cfg_path).read_text())
    reac = doc["propellant"]["reactants"]
    names = [r["name"] for r in reac]
    w = np.array([r["wt_percent"] for r in reac], float)
    t = np.array([r["temperature_K"] for r in reac], float)
    rmix = cea.Mixture(names)
    pmix = cea.Mixture(names, products_from_reactants=True, omit=[], ions=False)
    solver = cea.EqSolver(pmix, reactants=rmix, transport=False, ions=False, trace=1e-10, insert=[])
    sol = cea.EqSolution(solver)
    h = float(rmix.calc_property(cea.ENTHALPY, w, t))
    hR = h / float(cea.R)
    species = list(pmix.species_names)
    ng = int(solver.num_gas)
    gas = set(species[:ng])
    rows = []
    for p in PRESSURES_PA:
        solver.solve(sol, cea.HP, hR, p / 1e5, w)
        mf = sol.mass_fractions
        Yg = sum(float(mf[s]) for s in gas)
        T, M, n, MW, rho = float(sol.T), float(sol.M), float(sol.n), float(sol.MW), float(sol.density)
        nj = np.asarray(sol.nj, float)[:len(species)]
        sum_gas_nj = float(sum(nj[i] for i, s in enumerate(species) if s in gas))
        sum_all_nj = float(np.sum(nj))
        # independent gas-phase mean molar mass: gas mass per kg-total / gas kmol per kg-total
        M_gas = Yg / sum_gas_nj
        rows.append({
            "case": Path(cfg_path).stem, "p_Pa": p, "converged": bool(sol.converged),
            "T_K": T, "M": M, "n": n, "MW": MW, "density": rho, "Yg": Yg,
            "M_times_n_minus_1": M * n - 1.0,
            "density_over_pM_CEA_RuT_minus_1": rho / (p * M / (RU * T)) - 1.0,
            "density_over_pMgas_RuT_minus_1": rho / (p * M_gas / (RU * T)) - 1.0,
            "M_over_MW": M / MW,
            "sum_gas_nj": sum_gas_nj, "sum_all_nj": sum_all_nj,
            "n_over_sum_gas_nj_minus_1": n / sum_gas_nj - 1.0,
            "MW_times_sum_all_nj_minus_1": MW * sum_all_nj - 1.0,
            "M_gas_independent": M_gas,
            "M_over_Mgas_div_Yg_minus_1": M / (M_gas / Yg) - 1.0,
        })
    # sonic API location (attribute presence only on EqSolution; RocketSolution values)
    eq_has_sonic = hasattr(sol, "sonic_velocity")
    rocket = {}
    try:
        rs = cea.RocketSolver(pmix, reactants=rmix, transport=False, trace=1e-10)
        rsol = cea.RocketSolution(rs)
        rs.solve(rsol, w, 30.0, supar=[2.0], hc=hR, iac=True)
        k = int(rsol.num_pts)
        rocket = {"converged": bool(rsol.converged), "num_pts": k,
                  "Mach": [float(x) for x in np.asarray(rsol.Mach)[:k]],
                  "sonic_velocity_m_s": [float(x) for x in np.asarray(rsol.sonic_velocity)[:k]],
                  "T_K": [float(x) for x in np.asarray(rsol.T)[:k]]}
    except Exception as exc:  # recorded, not hidden
        rocket = {"error": f"{type(exc).__name__}: {exc}"}
    return rows, {"EqSolution_has_sonic_velocity": eq_has_sonic,
                  "EqSolution_public_attrs": sorted(a for a in dir(sol) if not a.startswith("_")),
                  "RocketSolution_30bar_supar2": rocket}


def main():
    out = Path(sys.argv[1]); out.mkdir(parents=True, exist_ok=True)
    allrows, api = [], {}
    for cfg in sys.argv[2:]:
        rows, a = probe_case(cfg)
        allrows += rows; api[Path(cfg).stem] = a
    keys = list(allrows[0])
    with open(out / "cea_semantics.csv", "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=keys); wr.writeheader(); wr.writerows(allrows)
    summary = {
        "cea_version": getattr(cea, "__version__", "unknown"),
        "cea_lib_version": str(cea.lib_version()),
        "Ru_J_kmolK": RU,
        "max_abs_M_times_n_minus_1": max(abs(r["M_times_n_minus_1"]) for r in allrows),
        "max_abs_density_vs_pM_CEA_RuT": max(abs(r["density_over_pM_CEA_RuT_minus_1"]) for r in allrows),
        "max_abs_n_vs_sum_gas_nj": max(abs(r["n_over_sum_gas_nj_minus_1"]) for r in allrows),
        "max_abs_MW_vs_sum_all_nj": max(abs(r["MW_times_sum_all_nj_minus_1"]) for r in allrows),
        "max_abs_M_vs_Mgas_over_Yg": max(abs(r["M_over_Mgas_div_Yg_minus_1"]) for r in allrows),
        "max_abs_density_vs_pMgas_RuT": max(abs(r["density_over_pMgas_RuT_minus_1"]) for r in allrows),
        "rows": allrows, "api": api,
    }
    (out / "cea_semantics.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps({k: v for k, v in summary.items() if k not in ("rows", "api")}, indent=1))
    for c, a in api.items():
        r = a["RocketSolution_30bar_supar2"]
        print(c, "EqSolution.sonic_velocity:", a["EqSolution_has_sonic_velocity"],
              "| Rocket:", r.get("error") or (r["Mach"], [round(x, 2) for x in r["sonic_velocity_m_s"]]))
    for r in allrows:
        print(f'{r["case"]:<26} p={r["p_Pa"]:.0e} Yg={r["Yg"]:.4f} M={r["M"]:.4f} 1/n={1/r["n"]:.4f} MW={r["MW"]:.4f} rho/(pM/RuT)-1={r["density_over_pM_CEA_RuT_minus_1"]:.2e}')


if __name__ == "__main__":
    main()
