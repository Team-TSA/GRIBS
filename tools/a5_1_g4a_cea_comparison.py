"""A-5.1 G4a candidate evidence: CEA-referenced semantic/numerical comparison.

GRIBS side : the GRIBS CEA backend table (PCHIP Psi, gamma_s, R_eff, Yg) and
             the GRIBS ideal choked-flow c* formula sqrt(Psi)/Gamma(gamma_s).
CEA side   : official ``cea`` RocketSolver/RocketSolution (IAC), computed here
             directly, not via GRIBS code.
Tolerances follow the 2026-10-09 guide (N-1, N-4 section 2):
  density identity 1e-14, sum of phase fractions 1e-13,
  sonic a = sqrt(gamma_s*Ru*T/M_CEA) vs RocketSolution 1e-12,
  GRIBS c* vs CEA IAC c* 2e-3 (10/09 N-1; different models: validation comparison).
CEA throat validity (added 2026-11-02): a point's IAC c* is a usable reference only if the
RocketSolution throat station is valid: converged, |Mach_t - 1| <= 1e-3, Ae/At at that
station == 1 (1e-9), and c*_CEA == p0/(rho_t*u_t) within 1e-6 (converged=True alone is not
enough: Mach_t ~0.90 with converged=True was observed when Al2O3 is near its melting point at
the throat).  Invalid points are listed as 'reference_invalid'; their c* is never replaced.
The overall c* check is True only if there is at least one point, every point has a valid
reference, and every point is within tolerance -- an invalid or empty reference set is never
reported as an overall pass based on the remaining points.
This produces EVIDENCE for G4a; it does not approve G4a.
Usage: python tools/a5_1_g4a_cea_comparison.py <outdir> <case.json>...
"""
import copy, hashlib, importlib.util, json, math, sys, tempfile
from pathlib import Path
import numpy as np
import cea

ROOT = Path(__file__).resolve().parents[1]
PRIVATE_MODULE_NAME = "gribs_a5_1_g4a_tool_body"


def load_gribs_body(root):
    """Load the GRIBS body of THIS repository privately.

    Never reuses ``sys.modules["gribs"]`` (it may come from another tree or an
    older file state) and never replaces it.  The module is executed from the
    exact bytes whose sha256 is recorded, so the reported hash is the code that ran.
    """
    candidates = sorted(Path(root).glob("GRIBS_v*.py"))
    if len(candidates) != 1:
        raise RuntimeError(f"expected exactly one GRIBS_v*.py in {root}, found {[c.name for c in candidates]}")
    path = candidates[0].resolve()
    source = path.read_bytes()
    spec = importlib.util.spec_from_file_location(PRIVATE_MODULE_NAME, path)
    module = importlib.util.module_from_spec(spec)
    previous = sys.modules.get(PRIVATE_MODULE_NAME)
    sys.modules[PRIVATE_MODULE_NAME] = module  # needed while executing (dataclasses look up their module)
    try:
        exec(compile(source, str(path), "exec"), module.__dict__)
    except BaseException:
        if previous is None:
            sys.modules.pop(PRIVATE_MODULE_NAME, None)
        else:
            sys.modules[PRIVATE_MODULE_NAME] = previous
        raise
    info = {"path": str(path), "sha256": hashlib.sha256(source).hexdigest(),
            "module_name": PRIVATE_MODULE_NAME, "reused_existing_module": False}
    return module, info


gribs, GRIBS_MODULE_INFO = load_gribs_body(ROOT)

PC_BAR = [10.0, 30.0, 70.0]
SUPAR = 4.0
TOL = {"density": 1e-14, "phase_sum": 1e-13, "sonic": 1e-12, "cstar": 2e-3}  # c* per 10/09 N-1 table
THROAT_TOL = {"mach": 1e-3, "own_throat_cstar": 1e-6, "area_ratio": 1e-9}


def gamma_fn(g):
    return math.sqrt(g) * (2.0 / (g + 1.0)) ** ((g + 1.0) / (2.0 * (g - 1.0)))


def throat_validity(converged, mach, area_ratio, cstar_cea, p0_pa, rho_t, u_t):
    """Classify the CEA throat station of one RocketSolution as a c* reference."""
    vals = [mach, area_ratio, cstar_cea, p0_pa, rho_t, u_t]
    finite = all(isinstance(v, (int, float)) and math.isfinite(v) for v in vals)
    own = (p0_pa / (rho_t * u_t) / cstar_cea - 1.0) if finite and rho_t > 0 and u_t > 0 and cstar_cea > 0 else float("nan")
    reasons = []
    if not converged:
        reasons.append("not_converged")
    if not finite:
        reasons.append("non_finite_station_value")
    else:
        if not abs(mach - 1.0) <= THROAT_TOL["mach"]:
            reasons.append("throat_mach_not_sonic")
        if not abs(area_ratio - 1.0) <= THROAT_TOL["area_ratio"]:
            reasons.append("throat_area_ratio_not_one")
        if not abs(own) <= THROAT_TOL["own_throat_cstar"]:
            reasons.append("cstar_inconsistent_with_own_throat")
    return {"converged": bool(converged), "mach": float(mach), "area_ratio": float(area_ratio),
            "own_throat_cstar_rel": float(own), "valid": not reasons, "invalid_reasons": reasons}


def aggregate_cstar(rows, tol=TOL["cstar"]):
    """Overall c* verdict.  Never passes on the remaining points when a reference is
    invalid, and never passes when there is no valid point at all."""
    valid = [r for r in rows if r["cea_throat"]["valid"]]
    invalid = [f'{r["case"]}@{r["pc_bar"]}' for r in rows if not r["cea_throat"]["valid"]]
    within = [r for r in valid if abs(r["cstar_rel_diff"]) <= tol]
    if not rows:
        status = "no_points"
    elif not valid:
        status = "no_valid_reference"
    elif invalid:
        status = "reference_invalid_points_present"
    elif len(within) == len(valid):
        status = "all_within"
    else:
        status = "exceeds"
    return {"passed": status == "all_within", "status": status, "n_points": len(rows),
            "n_valid_reference": len(valid), "n_within_among_valid": len(within),
            "reference_invalid_points": invalid,
            "max_abs_rel_diff_valid": max((abs(r["cstar_rel_diff"]) for r in valid), default=None)}


def case(cfg_path, work):
    doc = copy.deepcopy(json.loads(Path(cfg_path).read_text()))
    doc["output"]["directory"] = str(work / "out")
    doc["thermochemistry"]["cea_python"]["cache_directory"] = str(work / "cache")
    conf = gribs.configuration_from_document(doc).config
    b = gribs.initialize_thermochemistry(conf, work / "cache")
    reac = conf.reactants_spec_list
    names = [r.name for r in reac]
    w = np.array([r.wt_percent for r in reac], float)
    t = np.array([r.temperature_K for r in reac], float)
    rmix = cea.Mixture(names)
    pmix = cea.Mixture(names, products_from_reactants=True, omit=[], ions=False)
    hR = float(rmix.calc_property(cea.ENTHALPY, w, t)) / float(cea.R)
    Ru = float(cea.R)
    rows = []
    for pc in PC_BAR:
        rs = cea.RocketSolver(pmix, reactants=rmix, transport=False, trace=1e-10)
        sol = cea.RocketSolution(rs)
        rs.solve(sol, w, pc, supar=[SUPAR], hc=hR, iac=True)
        k = int(sol.num_pts)
        T = np.asarray(sol.T, float)[:k]; Mc = np.asarray(sol.M, float)[:k]
        gs = np.asarray(sol.gamma_s, float)[:k]; a = np.asarray(sol.sonic_velocity, float)[:k]
        P = np.asarray(sol.P, float)[:k]; rho = np.asarray(sol.density, float)[:k]
        cstar_cea = float(np.asarray(sol.c_star, float)[0])
        Ma = np.asarray(sol.Mach, float)[:k]; ae = np.asarray(sol.ae_at, float)[:k]
        it = int(np.argmin(np.abs(Ma - 1.0)))
        throat = throat_validity(bool(sol.converged), float(Ma[it]), float(ae[it]), cstar_cea,
                                 float(P[0]) * 1e5, float(rho[it]), float(Ma[it] * a[it]))
        throat["station"] = it
        sonic_err = float(np.max(np.abs(np.sqrt(gs * Ru * T / Mc) / a - 1.0)))
        dens_err = float(np.max(np.abs(rho / (P * 1e5 * Mc / (Ru * T)) - 1.0)))
        p = pc * 1e5
        psi, _ = b.psi_and_derivative(p)
        g_tab = float(b.gi(math.log(p)))
        yg, yc = b.phase_fractions(p)
        cstar_gribs = math.sqrt(psi) / gamma_fn(g_tab)
        rows.append({
            "case": Path(cfg_path).stem, "pc_bar": pc, "Yg": yg,
            "phase_sum_err": abs(yg + yc - 1.0),
            "sonic_identity_max_rel_err_stations": sonic_err,
            "density_identity_max_rel_err_stations": dens_err,
            "cea_throat": throat,
            "cstar_cea_iac_m_s": cstar_cea, "cstar_gribs_formula_m_s": cstar_gribs,
            "cstar_rel_diff": cstar_gribs / cstar_cea - 1.0,
            "cstar_if_extra_div_yg_rel_diff": cstar_gribs * yg / cstar_cea - 1.0,
            "cstar_if_gas_R_then_div_yg_rel_diff": cstar_gribs * math.sqrt(yg) / cstar_cea - 1.0,
        })
    return rows


def main():
    out = Path(sys.argv[1]); out.mkdir(parents=True, exist_ok=True)
    rows = []
    for c in sys.argv[2:]:
        rows += case(c, Path(tempfile.mkdtemp(dir=out)))
    res = build_result(rows)
    (out / "g4a_cea_comparison.json").write_text(json.dumps(res, indent=1))
    print(f'gribs body: {GRIBS_MODULE_INFO["path"]} sha256={GRIBS_MODULE_INFO["sha256"]}')
    print(json.dumps(res["checks_within_tolerance"]))
    agg = res["cstar_aggregate"]
    print(f'cstar aggregate: status={agg["status"]} valid_reference={agg["n_valid_reference"]}/{agg["n_points"]}')
    if agg["reference_invalid_points"]:
        print("reference_invalid (overall c* check cannot pass):", ", ".join(agg["reference_invalid_points"]))
    for r in rows:
        print(f'{r["case"]:<26} pc={r["pc_bar"]:>4} Yg={r["Yg"]:.4f} c*_rel={r["cstar_rel_diff"]:+.2e} '
              f'(/Yg bug {r["cstar_if_extra_div_yg_rel_diff"]:+.3f}, sqrtYg {r["cstar_if_gas_R_then_div_yg_rel_diff"]:+.3f}) '
              f'sonic={r["sonic_identity_max_rel_err_stations"]:.1e} dens={r["density_identity_max_rel_err_stations"]:.1e} '
              f'throat={"valid" if r["cea_throat"]["valid"] else "INVALID"} Mt={r["cea_throat"]["mach"]:.6f}')


def build_result(rows):
    """Assemble the JSON result.  Every check is False for an empty row set."""
    def mx(key):
        return max((r[key] for r in rows), default=float("inf"))
    agg = aggregate_cstar(rows)
    checks = {
        "phase_sum": bool(rows) and mx("phase_sum_err") <= TOL["phase_sum"],
        "sonic": bool(rows) and mx("sonic_identity_max_rel_err_stations") <= TOL["sonic"],
        "density": bool(rows) and mx("density_identity_max_rel_err_stations") <= TOL["density"],
        "cstar_reference_valid": agg["n_points"] > 0 and agg["n_valid_reference"] == agg["n_points"],
        "cstar": agg["passed"],
    }
    return {"cea_version": cea.__version__, "gribs_module": dict(GRIBS_MODULE_INFO),
            "tolerances": TOL, "throat_validity_tolerances": THROAT_TOL,
            "supar": SUPAR, "checks_within_tolerance": checks, "cstar_aggregate": agg, "rows": rows,
            "note": "evidence only; G4a approval requires the owner's review"}


if __name__ == "__main__":
    main()
