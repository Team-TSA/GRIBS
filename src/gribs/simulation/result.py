"""Result of one simulation (A-6 unit 10, public API U-39 (a)).

Result holds the values that the command line program writes to summary.json and time_history.csv, computed by
the same functions in the same order (gribs.api.run_simulation). The field names follow the keys of summary.json;
`summary` is the "results" section. Solver objects (the ODE solutions) are not part of a Result.

Not included in stage 1 (values that only the command line computes, in gribs.cli): "environment",
"cea_theoretical_diagnostics" and "output_files" of summary.json.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List


@dataclass(frozen=True)
class Result:
    #: {"name", "version", "schema_version"}, as summary.json "program"
    program: Dict[str, str]
    #: absolute path of the configuration file
    configuration_source: str
    #: {"applied": bool, "notes": [str, ...]}, as summary.json
    configuration_migration: Dict[str, Any]
    #: the loaded document (after any schema migration, before a backend override)
    original_configuration: Dict[str, Any]
    #: dataclasses.asdict of the runtime configuration after the run
    resolved_configuration: Dict[str, Any]
    derived_inputs: Dict[str, Any]
    two_phase_model: Dict[str, Any]
    #: provenance of the thermochemistry backend (backend.metadata())
    thermochemistry: Dict[str, Any]
    #: performance summary, ledger, verdicts and warnings; summary.json "results"
    summary: Dict[str, Any]
    #: sampled time histories (arrays keyed by column); time_history.csv is written from these
    history: Dict[str, Any]
    #: ODE right-hand-side calls and property-table restarts of the run
    integration_budget: Dict[str, Any]
    #: combined post-run, ledger and applicability verdict (the command line exits non-zero when False)
    run_passed: bool

    @property
    def warnings(self) -> List[Dict[str, Any]]:
        """The machine-readable warnings of the run (summary["warnings"])."""
        return self.summary.get("warnings", [])

    @property
    def program_version(self) -> str:
        return self.program["version"]
