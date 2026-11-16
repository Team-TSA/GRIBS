"""run_simulation: one simulation of a configuration file, without the command line (A-6 unit 10, U-39 (a)).

run_simulation calls the functions that gribs.cli.main calls, in the same order: load_json_configuration, validate,
initialize_thermochemistry, model_applicability_error, run_model, sample, summarize, the same summary warnings,
verdict and run-state keys, backend.metadata() and the derived inputs. The blocks between the "same text as
gribs.cli.main" markers are copies of lines of main; tests/test_a6_public_api.py fails if a block is no longer
found in main. main itself is not changed and does not call run_simulation.

Differences from main:
- nothing is printed and no output file is written (no summary, history, figure); migration notes are in the
  Result. The thermochemistry backend still reads and writes its property cache, in the directory main uses
  (<output directory>/thermo_cache, or thermochemistry.cea_python.cache_dir).
- errors are raised, not turned into exit codes: ConfigurationError (main: exit 2), ThermochemistryError (exit 3),
  ModelApplicabilityError (exit 1, before the solver), and the exceptions of the solver (exit 4).
- the backend state of gribs.thermo is restored to what it was before the call, also when an error is raised.
- gribs.cli, gribs.plotting and gribs.selftest are not imported, so matplotlib is not imported.
- not computed: the environment record, the CEA theoretical rocket diagnostics and the post-run self-test display
  (main computes them in gribs.cli and gribs.selftest).

Names that tests patch are read through their owner module (ADR decision 1), e.g.
gribs.simulation.integrator.run_model, so a patch on the owner reaches run_simulation as it reaches main.
"""
from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, Optional, Union

import gribs.config
import gribs.config.loader
import gribs.simulation.integrator
import gribs.simulation.postprocess
import gribs.thermo
from gribs.constants import PROGRAM_NAME, PROGRAM_VERSION, SCHEMA_VERSION
from gribs.errors import GribsError
from gribs.simulation.postprocess import add_run_state_keys, append_run_verdict_warnings, run_verdict_passed
from gribs.simulation.result import Result


class ModelApplicabilityError(GribsError):
    """The configuration is outside the scope of the two-phase model; the solver was not started.

    ``code`` and ``context`` are those of the structured record main prints (context["simulation_started"] is
    False).
    """

    def __init__(self, record: Dict[str, Any]):
        super().__init__(f"[{record['code']}] {record['message']}")
        self.code = record["code"]
        self.context = {**record["context"], "simulation_started": False}


def run_simulation(config: Union[str, Path], *, backend_override: Optional[str] = None) -> Result:
    """Run the configuration file `config` and return its Result.

    `backend_override` has the meaning of the --backend option of the command line.
    """
    previous = gribs.thermo._thermo_state()
    try:
        return _run(Path(config), backend_override)
    finally:
        gribs.thermo._set_thermo(previous)


def _run(path: Path, backend_override: Optional[str]) -> Result:
    if backend_override:
        loaded = gribs.config.loader.load_json_configuration(path,
                                         backend_override=backend_override)
    else:
        loaded = gribs.config.loader.load_json_configuration(path)

    cfg = loaded.config
    gribs.config.validate(cfg)

    cache_dir = loaded.output_dir / "thermo_cache"
    backend = gribs.thermo.initialize_thermochemistry(cfg, cache_dir)

    applicability_error = gribs.simulation.postprocess.model_applicability_error(cfg, backend)
    if applicability_error is not None:
        raise ModelApplicabilityError(applicability_error)

    res = gribs.simulation.integrator.run_model(cfg)
    hist = gribs.simulation.postprocess.sample(res, cfg)
    # --- begin: same text as gribs.cli.main ---
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
    add_run_state_keys(summ)
    # --- end: same text as gribs.cli.main ---

    # --- begin: same text as gribs.cli.main ---
    thermo_md = backend.metadata()
    # --- end: same text as gribs.cli.main ---

    # --- begin: same text as gribs.cli.main ---
    derived = {
        "propellant_density_kg_m3": cfg.rho_p,
        "ideal_mixture_density_kg_m3": backend.ideal_mixture_density,
        "packing_fraction": backend.packing_fraction,
        "bulk_propellant_density_kg_m3": backend.bulk_mixture_density,
        "propellant_density_model": (
            "rho_p = packing_fraction / sum_i(w_i/rho_i); independent of any CEA gas density"),
    }
    # --- end: same text as gribs.cli.main ---
    return Result(
        program={"name": PROGRAM_NAME, "version": PROGRAM_VERSION, "schema_version": SCHEMA_VERSION},
        configuration_source=cfg.configuration_source,
        configuration_migration=({"applied": True, "notes": list(loaded.migration_notes)}
                                 if loaded.migration_notes else {"applied": False, "notes": []}),
        original_configuration=loaded.document,
        resolved_configuration=asdict(cfg),
        derived_inputs=derived,
        two_phase_model=summ["two_phase_model"],
        thermochemistry=thermo_md,
        summary=summ,
        history=hist,
        integration_budget=dict(res.get("integration_budget") or {}),
        run_passed=run_passed,
    )
