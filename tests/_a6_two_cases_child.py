"""Test helper (not collected): configuration cases one after another in ONE process, WITHOUT restoring the
thermochemistry backend between them (A-6 unit 3a; ADR_A6_PACKAGE_SPLIT, report 112 addendum).
argv: <body.py> <work dir> <case> [<case> ...]; cases are keys of PROPERTIES. Prints one JSON line per case."""
import hashlib, importlib.util, json, sys
from pathlib import Path
body, work, cases = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3:]
spec = importlib.util.spec_from_file_location("gribs_two_cases", body); m = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = m; spec.loader.exec_module(m)
# reference_constant properties (temperature K, gamma, molecular weight kg/kmol): no CEA, deterministic
PROPERTIES = {"A": (3000.0, 1.2, 25.0), "B": (2400.0, 1.25, 21.0)}


def document(case):
    doc = json.loads((body.parent / "gribs_config.json").read_text(encoding="utf-8"))
    th = doc["thermochemistry"]
    th["backend"] = "reference_constant"
    th.pop("cea_python", None); th.pop("cea_legacy_executable", None)
    t0, gamma, mw = PROPERTIES[case]
    th["reference_constant"] = {"temperature_K": t0, "gamma": gamma, "molecular_weight_kg_kmol": mw}
    return doc


installed = []  # (backend, case that installed it); the objects are kept alive, so identity is unambiguous


def installer(backend):
    return next((case for b, case in installed if b is backend), "unknown")


for position, case in enumerate(cases):
    c = m.configuration_from_document(document(case)).config
    before = m._thermo_state()
    before_case = None if before is None else installer(before)
    backend = m.initialize_thermochemistry(c, work / f"{position}_{case}")
    if installer(backend) == "unknown":
        installed.append((backend, case))
    result = m.run_model(c)
    digest = hashlib.sha256()
    for key in ("sol_burn", "sol_blow"):
        sol = result[key]
        if sol is not None:
            digest.update(key.encode() + b"\0" + sol.t.tobytes() + sol.y.tobytes())
    print(json.dumps({
        "case": case, "position": position,
        "state_before": before_case,
        "installed": m._thermo_state() is backend and m.thermochemistry() is backend,
        "rho_p": repr(c.rho_p), "stop_burn": str(result["stop_burn"]), "n_burn": int(result["sol_burn"].t.size),
        "digest": digest.hexdigest(),
    }))
