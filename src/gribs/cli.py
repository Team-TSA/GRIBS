"""Command line and main program (A-6 unit 9). Moved verbatim from section 13 of the entry: environment_info,
build_parser, _print_backend_banner, _cea_theoretical_diagnostics and main. The entry re-exports every name and
calls main when it is run as a script.

Tests patch _print_backend_banner, environment_info and _cea_theoretical_diagnostics here, e.g.
"gribs.cli.environment_info"; main reads them as module globals of this module, so a patch on the entry would
not reach them. Names of other package modules that tests patch are read as <module>.<name> (ADR decision 1);
the two reads of load_json_configuration (gribs.config.loader) are two of the three lines that differ from the
entry text. The third is the default of --config: gribs_config.json in _PROGRAM_DIR (gribs.thermo: the
directory that holds src/, i.e. the entry's directory) instead of the directory of the entry file itself.
"""
from __future__ import annotations

import argparse
import json
import math
import platform
import sys
from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, Optional, Sequence

import numpy as np
import matplotlib

import gribs.config
import gribs.config.loader
import gribs.plotting
import gribs.reporting
import gribs.selftest
import gribs.simulation.integrator
import gribs.simulation.postprocess
import gribs.thermo
from gribs.config import Config
from gribs.constants import (
    BACKEND_CEA_LEGACY_EXECUTABLE,
    BACKEND_CEA_PYTHON,
    BACKEND_REFERENCE_CONSTANT,
    KNOWN_BACKENDS,
    PROGRAM_NAME,
    PROGRAM_VERSION,
    SCHEMA_VERSION,
    TWO_PHASE_HOMOGENEOUS,
)
from gribs.errors import ConfigurationError, GribsError, ThermochemistryError
from gribs.reporting import cea_diagnostics_status
from gribs.selftest import self_tests
from gribs.simulation.postprocess import (
    EXIT_CODE_SUCCESS,
    append_run_verdict_warnings,
    run_verdict_passed,
    run_verification_exit_code,
)
from gribs.thermo import _PROGRAM_DIR, ThermochemistryBackend


def environment_info() -> Dict[str, Any]:
    info: Dict[str, Any] = {
        "program": PROGRAM_NAME,
        "program_version": PROGRAM_VERSION,
        "schema_version": SCHEMA_VERSION,
        "python_version": sys.version.split()[0],
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "system": platform.system(),
        "machine": platform.machine(),
        "numpy_version": np.__version__,
        "scipy_version": __import__("scipy").__version__,
        "matplotlib_version": matplotlib.__version__,
    }
    try:
        # only for provenance; never lets the import load a database (P5-1)
        cea = gribs.thermo.import_cea_without_auto_init()
        info["cea_python_version"] = getattr(cea, "__version__", None)
        info["cea_python_module_path"] = getattr(cea, "__file__", None)
        info["cea_library_version"] = str(cea.lib_version())
    except Exception as exc:
        info["cea_python_version"] = None
        info["cea_python_import_error"] = f"{type(exc).__name__}: {exc}"
    return info


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=f"GRIBS v{PROGRAM_VERSION} JSON-configured internal ballistics")
    p.add_argument("--config", type=Path,
                   default=_PROGRAM_DIR / "gribs_config.json",
                   help="complete JSON configuration file")
    p.add_argument("--backend", choices=list(KNOWN_BACKENDS), default=None,
                   help="override thermochemistry.backend for this run")
    p.add_argument("--selftest", action="store_true", help="run numerical self-tests only")
    p.add_argument("--validate-config", action="store_true",
                   help="validate the JSON configuration and exit")
    p.add_argument("--dump-thermo", type=Path, default=None, metavar="FILE",
                   help="write the chamber-property table to FILE and exit")
    p.add_argument("--no-cache", action="store_true",
                   help="force a rebuild of the chamber-property table")
    p.add_argument(
        "--allow-verification-failure",
        action="store_true",
        help=(
            "return exit 0 only for a mass-ledger-only verification failure "
            "(diagnostics remain marked failed)"
        ),
    )
    p.add_argument("--version", action="version",
                   version=f"{PROGRAM_NAME} {PROGRAM_VERSION}")
    return p


def _print_backend_banner(cfg: Config, backend: ThermochemistryBackend) -> None:
    print(f"{PROGRAM_NAME} v{PROGRAM_VERSION}  (schema {SCHEMA_VERSION})")
    print(f"thermochemistry backend : {backend.name}")
    two_phase_active = (cfg.two_phase_mode == TWO_PHASE_HOMOGENEOUS
                        and bool(getattr(backend, "supports_two_phase", False)))
    print(f"  two-phase model       : {cfg.two_phase_mode}"
          + ("  (ACTIVE: total product mass conserved, EOS on the gas phase)"
             if two_phase_active else
             ("  (phase split inactive: Yg := 1)"
              if cfg.two_phase_mode == TWO_PHASE_HOMOGENEOUS else
              "  (v0.4.3 regression pathway)")))
    if two_phase_active:
        print(f"  nozzle entrainment    : {cfg.two_phase_nozzle_entrainment} "
              f"(no particle slip)")
        print(f"  condensed volume      : {cfg.two_phase_condensed_volume}")
    if backend.name == BACKEND_CEA_PYTHON:
        print(f"  official CEA package  : {getattr(backend.cea, '__version__', '?')} "
              f"at {getattr(backend.cea, '__file__', '?')}")
        print(f"  CEA library version   : {backend.cea.lib_version()}")
        print(f"  equilibrium           : HP (assigned enthalpy and pressure)")
    elif backend.name == BACKEND_CEA_LEGACY_EXECUTABLE:
        print("  COMPATIBILITY BACKEND (transitional, v0.4-alpha only).")
        print("  Prefer 'cea_python' (official package) for new work.")
    elif backend.name == BACKEND_REFERENCE_CONSTANT:
        print("  property model        : constant analytic properties")
        print(
            "  reference properties : "
            f"T0 = {backend._temperature:g} K, "
            f"gamma = {backend._gamma:g}, "
            f"M = {backend._molecular_weight:g} kg/kmol"
        )
        print("  CEA used              : no")
        print("  property cache        : disabled")
        print(
            f"  validity range        : {backend.pmin/1e5:g} - "
            f"{backend.pmax/1e5:g} bar, "
            f"policy = {backend.outside_range_policy}"
        )
    else:
        print(f"  property table        : {backend.c.p_fit_min/1e5:g} - "
              f"{backend.c.p_fit_max/1e5:g} bar, "
              f"{backend.c.cea_pressure_points} points, "
              f"eta_T0 = {backend.c.eta_T0:g}, "
              f"policy = {backend.c.property_policy}")
    if getattr(backend, "cache_path", None):
        print(f"  cache file            : {backend.cache_path}")
        print(f"  cache reused          : {backend.cache_used}")
    print(f"  ideal mixture density : {backend.ideal_mixture_density:.6f} kg/m3")
    print(f"  packing fraction      : {backend.packing_fraction:.6f}")
    print(f"  bulk rho_p used       : {cfg.rho_p:.6f} kg/m3")


def _cea_theoretical_diagnostics(cfg: Config, backend: ThermochemistryBackend,
                                 summ: dict) -> Optional[dict]:
    """Theoretical CEA rocket comparison values (never fed back into GRIBS).

    A failed, incomplete or invalid diagnostic is returned with an explicit
    ``cea_theoretical_status`` (plus the error or the missing / invalid field
    names) so that it is preserved in summary.json and shown in the text
    summary and the figure; a successful result is returned unchanged.
    """
    if not getattr(cfg, "cea_py_rocket_diagnostics", False):
        return None
    if not hasattr(backend, "theoretical_rocket_performance"):
        return None
    p_ref = summ.get("p_mean_burn_Pa")
    if not (isinstance(p_ref, float) and math.isfinite(p_ref) and p_ref > 0.0):
        p_ref = summ.get("p_max_Pa")
    try:
        result = backend.theoretical_rocket_performance(float(p_ref), float(cfg.eps_nozzle))
    except Exception as exc:                                   # diagnostic only
        print(f"WARNING: CEA theoretical rocket diagnostics unavailable: "
              f"{type(exc).__name__}: {exc}")
        return {"cea_theoretical_status": "failed",
                "cea_theoretical_error": f"{type(exc).__name__}: {exc}"}
    status, missing, invalid = cea_diagnostics_status(result)
    if status == "ok":
        return result
    annotated = dict(result) if isinstance(result, dict) else {"cea_theoretical_raw": repr(result)[:200]}
    annotated["cea_theoretical_status"] = status
    if missing:
        annotated["cea_theoretical_missing_fields"] = list(missing)
    if invalid:
        annotated["cea_theoretical_invalid_fields"] = list(invalid)
    print(f"WARNING: CEA theoretical rocket diagnostics {status}: "
          f"missing={list(missing)} invalid={list(invalid)}")
    return annotated


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        if args.backend:
            loaded = gribs.config.loader.load_json_configuration(args.config,
                                             backend_override=args.backend)
        else:
            loaded = gribs.config.loader.load_json_configuration(args.config)
    except ConfigurationError as exc:
        print(f"Configuration error:\n{exc}", file=sys.stderr)
        return 2

    cfg = loaded.config
    if args.no_cache:
        cfg.cea_rebuild_cache = True
        loaded.migration_notes.append("Chamber-property cache rebuild forced by --no-cache.")

    for note in loaded.migration_notes:
        print(f"WARNING: {note}", file=sys.stderr)

    try:
        gribs.config.validate(cfg)
    except (ConfigurationError, ValueError) as exc:
        print(f"Configuration error:\n{exc}", file=sys.stderr)
        return 2

    if args.validate_config:
        print(f"Configuration {loaded.config.configuration_source} is valid "
              f"(schema {SCHEMA_VERSION}, backend {cfg.thermo_backend}).")
        return 0

    cache_dir = loaded.output_dir / "thermo_cache"
    try:
        backend = gribs.thermo.initialize_thermochemistry(cfg, cache_dir)
    except (ThermochemistryError, ConfigurationError) as exc:
        error_code = getattr(exc, "code", None)
        if error_code:
            print(
                f"Thermochemistry backend error [{error_code}]:\n{exc}",
                file=sys.stderr,
            )
            error_record = {
                "code": error_code,
                "severity": "error",
                "message": str(exc),
                "context": getattr(exc, "context", {}),
            }
            print(
                json.dumps(error_record, ensure_ascii=False, sort_keys=True),
                file=sys.stderr,
            )
        else:
            print(f"Thermochemistry backend error:\n{exc}", file=sys.stderr)
        return 3

    _print_backend_banner(cfg, backend)

    if args.dump_thermo is not None:
        payload = {
            "backend": backend.name,
            "provenance": backend.metadata(),
            **backend.thermo_dump_table(),
        }
        args.dump_thermo.write_text(
            json.dumps(payload, indent=2) + "\n",
            encoding="utf-8",
        )
        print(f"Chamber-property table written to {args.dump_thermo}")
        return 0

    if args.selftest:
        return 0 if self_tests(cfg) else 1

    applicability_error = gribs.simulation.postprocess.model_applicability_error(cfg, backend)
    if applicability_error is not None:
        print(
            f"Model applicability error [{applicability_error['code']}]: "
            f"{applicability_error['message']}",
            file=sys.stderr,
        )
        error_record = dict(applicability_error)
        error_record["context"] = {
            **applicability_error["context"],
            "simulation_started": False,
        }
        print(
            json.dumps(error_record, ensure_ascii=False, sort_keys=True),
            file=sys.stderr,
        )
        return 1

    print("Running the internal-ballistics simulation ...")
    try:
        res = gribs.simulation.integrator.run_model(cfg)
        hist = gribs.simulation.postprocess.sample(res, cfg)
    except (GribsError, ValueError, RuntimeError) as exc:
        print(f"Simulation failed: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 4
    _ib = (res.get("integration_budget") if isinstance(res, dict) else None) or {}
    if _ib.get("table_restarts_burn") or _ib.get("table_restarts_blowdown"):
        # Diagnostic record only (stderr, not JSON-first so that structured error parsers ignore it); the
        # summary / CSV schema is unchanged.
        print("INTEGRATION_TABLE_RESTARTS " + json.dumps(_ib, sort_keys=True), file=sys.stderr)
    summ = gribs.simulation.postprocess.summarize(res, hist, cfg, backend=backend)
    # If migration from v0.4.3 changed the physical model, record it as a
    # machine-readable warning (previously this only appeared in stderr and
    # configuration_migration.notes, making it invisible to automated checks).
    if loaded.migration_notes:
        _model_change_keywords = ("THIS CHANGES THE PHYSICAL MODEL",)
        for _note in loaded.migration_notes:
            if any(kw in _note for kw in _model_change_keywords):
                summ.setdefault("warnings", []).append({
                    "code": "W_MIGRATION_MODEL_CHANGE",
                    "severity": "warning",
                    "message": (
                        "Automatic configuration migration changed the physical "
                        "model.  See configuration_migration.notes in the summary "
                        "for details."
                    ),
                    "context": {"migration_note": _note},
                })
                break
    _cache_dur = getattr(backend, "cache_write_durability", None)      # cache durability warning
    if _cache_dur:
        summ.setdefault("warnings", []).append({
            "code": "W_THERMO_CACHE_DURABILITY_NOT_GUARANTEED", "severity": "warning",
            "message": ("The chamber-property cache was written, but fsync of its directory is not supported on "
                        "this file system; the cache file's durability after a crash or power loss is not "
                        "guaranteed."),
            "context": dict(_cache_dur)})
    append_run_verdict_warnings(summ)
    run_passed = run_verdict_passed(summ)
    verification_exit_code = run_verification_exit_code(
        summ,
        allow_verification_failure=args.allow_verification_failure,
    )

    thermo_md = backend.metadata()
    cea_diag = _cea_theoretical_diagnostics(cfg, backend, summ)

    outdir = loaded.output_dir
    outdir.mkdir(parents=True, exist_ok=True)
    png = outdir / loaded.output_names["figure_filename"]
    csvp = outdir / loaded.output_names["history_filename"]
    txtp = outdir / loaded.output_names["summary_text_filename"]
    jsnp = outdir / loaded.output_names["summary_json_filename"]

    text = gribs.reporting.summary_text(
        summ,
        cfg,
        thermo_md,
        cea_diag,
        loaded.migration_notes,
    )

    derived = {
        "propellant_density_kg_m3": cfg.rho_p,
        "ideal_mixture_density_kg_m3": backend.ideal_mixture_density,
        "packing_fraction": backend.packing_fraction,
        "bulk_propellant_density_kg_m3": backend.bulk_mixture_density,
        "propellant_density_model": (
            "rho_p = packing_fraction / sum_i(w_i/rho_i); independent of any CEA gas density"),
    }
    summary_json = {
        "program": {"name": PROGRAM_NAME, "version": PROGRAM_VERSION,
                    "schema_version": SCHEMA_VERSION},
        "environment": environment_info(),
        "configuration_source": cfg.configuration_source,
        "configuration_migration": ({"applied": True, "notes": list(loaded.migration_notes)}
                                    if loaded.migration_notes else {"applied": False,
                                                                    "notes": []}),
        "original_configuration": loaded.document,
        "resolved_configuration": asdict(cfg),
        "derived_inputs": derived,
        "two_phase_model": summ["two_phase_model"],
        "thermochemistry": thermo_md,
        "cea_theoretical_diagnostics": cea_diag,
        "results": summ,
        "output_files": {"figure": str(png), "history_csv": str(csvp),
                         "summary_text": str(txtp), "summary_json": str(jsnp)},
    }

    # Preserve machine-readable diagnostics before any verdict-based exit
    # and before figure generation.
    gribs.reporting.write_csv(hist, csvp)
    txtp.write_text(text, encoding="utf-8")
    jsnp.write_text(
        json.dumps(summary_json, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    gribs.plotting.make_figure(res, hist, summ, cfg, png, cea_diag)

    gribs.selftest.self_tests_post_run(
        res,
        hist,
        cfg,
        report=summ["post_run_consistency"],
    )

    print()
    print(text)
    print("Output files")
    for pth in (png, csvp, txtp, jsnp):
        print(f"  {pth}")

    if not run_passed:
        failed_codes = [
            warning["code"]
            for warning in summ["warnings"]
            if warning.get("code") in {
                "W_POST_RUN_CONSISTENCY",
                "W_MASS_LEDGER_FAILED",
                "E_OUT_OF_VALIDATED_RANGE",
            }
        ]
        if verification_exit_code == EXIT_CODE_SUCCESS:
            print(
                "Run verification failure allowed by "
                "--allow-verification-failure after diagnostic output was "
                "written: " + ", ".join(failed_codes),
                file=sys.stderr,
            )
            return EXIT_CODE_SUCCESS
        print(
            "Run verification failed after diagnostic output was written: "
            + ", ".join(failed_codes),
            file=sys.stderr,
        )
        return verification_exit_code

    return EXIT_CODE_SUCCESS
