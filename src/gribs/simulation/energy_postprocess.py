"""Energy-mode run, time history, summary and output files (A-6 stage 2, C12; contract 164; reports 188, 189).

gribs.cli.main and gribs.api._run call run_energy_case right after the configuration is read and validated, with the same
text (a lazy import, only for a configuration that has an energy_model section).  run_energy_case does the whole computation
and writes nothing: it returns an EnergyOutcome.  The command line writes the files with EnergyOutcome.write_and_report; the API
builds a Result from the same fields.  gribs.simulation.postprocess and gribs.simulation.integrator.run_model are not used: an
energy run does not go through the legacy thermochemistry table, the legacy sampling, the legacy summary or the legacy figure.

SCOPE: energy mode is verified on SYNTHETIC fixtures only (an analytic reference_constant provider, ideal-gas or synthetic
providers, and one smoke run with the CEA provider of a golden propellant that checks completion, gamma_s > 1 and the closure
residuals).  No real propellant is accepted and no result is compared with measured data or, for a condensed volume, with the
legacy model.  The summary says so (energy_model.scope).

Output (D-05 (a): a dedicated minimal output): energy_time_history.csv (columns with units, the closure residuals R_M, R_B,
R_E, R_V in absolute and normalised form), summary.json and summary.txt under the configured file names with the top-level keys
model_mode = "energy" and energy_model (the interpretation, the definition of H0_out, the held items adopted, gamma_s range,
provenance, the settings, and legacy_outputs_not_produced).  The legacy figure and diagnostics are not produced, and the command
line and the summary say so.  A run passes when it completes, gamma_s > 1 in every state (a state with a gamma_s that is not
finite or not above 1 stops the run) and the normalised closure residuals are within CLOSURE_TOLERANCE.
"""
from __future__ import annotations

import json
import math
import sys
import traceback
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np

import gribs.simulation.energy_integrator as ei
from gribs.constants import PROGRAM_NAME, PROGRAM_VERSION, SCHEMA_VERSION
from gribs.energy.ledger import closure_residuals
from gribs.errors import ConfigurationError, GribsError, ThermochemistryError

#: The pass criterion of the closure residuals of an energy run: |R| / scale of each of R_M, R_B, R_E, R_V (a choice of this
#: commit, recorded in the summary; the residuals are accumulated with the same right-hand side, so they sit near 1e-13).
CLOSURE_TOLERANCE = 1.0e-6
CLOSURE_TOLERANCE_STATUS = ("provisional (held item H-09): a pass/fail threshold, not a claim of accuracy; to be fixed in the "
                            "acceptance of real cases")
CLOSURE_WARNING_CODE = "W_ENERGY_CLOSURE_FAILED"
EXIT_CODE_SUCCESS = 0
EXIT_CODE_CONFIGURATION = 2
EXIT_CODE_THERMOCHEMISTRY = 3
EXIT_CODE_SIMULATION = 4
EXIT_CODE_VERIFICATION_FAILED = 1
UNEXPECTED_ERROR_CODE = "E_ENERGY_UNEXPECTED_ERROR"
OUTPUT_WRITE_ERROR_CODE = "E_ENERGY_OUTPUT_WRITE_FAILED"
HISTORY_FILENAME = "energy_time_history.csv"

HISTORY_COLUMNS = (
    "time_s", "mass_kg", "total_energy_J", "pressure_Pa", "temperature_K", "gamma_s", "burn_depth_m", "throat_radius_m",
    "free_volume_m3", "burn_rate_m_s", "mdot_burn_kg_s", "mdot_out_kg_s", "thrust_N", "generated_mass_kg",
    "discharged_mass_kg", "impulse_Ns", "burn_surface_energy_J", "discharge_energy_J", "wall_heat_J",
    "R_M_kg", "R_B_kg", "R_E_J", "R_V_m3", "R_M_normalised", "R_B_normalised", "R_E_normalised", "R_V_normalised")

LEGACY_OUTPUTS_NOT_PRODUCED = (
    "figure (the ballistics figure of the legacy model)",
    "legacy time history (the legacy time_history.csv columns; see energy_time_history.csv)",
    "CEA theoretical rocket diagnostics",
    "equilibrium pressure diagnostic",
    "post-run consistency self-test of the legacy model",
    "legacy mass and energy ledger",
    "legacy two-phase model applicability check",
)

BURN_SURFACE_INTERPRETATION = (
    "A: the burn-surface term is written on the solid side, form (a): mdot_burn*h_s - p*A_b*r, with h_s evaluated at the far "
    "temperature of the grain (T_grain_inf) on the absolute enthalpy basis; heat through the burn surface is included by that "
    "evaluation and not added separately (contract 164 sections 1b, 2b-4.3)")
H0_OUT_DEFINITION = (
    "H0_out = h + v_inlet^2/2 with h = u + p*V_free/M (u = E_cv/M, the contents at rest) and v_inlet = 0, so the kinetic energy at "
    "the nozzle inlet plane is NOT included: this departs from the model definition (contract 164 section 2b-4.4 (i)), which "
    "requires the kinetic energy at the inlet plane to be included; the inlet plane (position, area) is not defined (held item "
    "H-18). The stagnation enthalpy per kg that leaves is not changed by Cd (an adiabatic nozzle: held item H-14)")
HELD_ITEMS_ADOPTED = (
    "H-18 (the kinetic energy at the nozzle inlet plane is not included; departs from contract 164 section 2b-4.4 (i))",
    "H-07 (quasi-steady burn surface: no storage term for ignition and tail-off)",
    "H-14 (adiabatic nozzle: Cd scales the mass flow only)",
    "H-17 (the atomic-mass table is not audited against the atomic masses of CEA)")
SCOPE = ("energy mode is verified on synthetic fixtures only; no real propellant is accepted, and no result is compared with "
         "measured data or (for a condensed volume) with the legacy model")


@dataclass
class EnergyOutcome:
    """Everything an energy run produces; run_energy_case builds it, the command line and the API consume it."""
    error: Optional[BaseException]
    program: Dict[str, str] = field(default_factory=dict)
    derived_inputs: Dict[str, Any] = field(default_factory=dict)
    two_phase_model: Dict[str, Any] = field(default_factory=dict)
    thermochemistry: Dict[str, Any] = field(default_factory=dict)
    results: Dict[str, Any] = field(default_factory=dict)
    energy_model: Dict[str, Any] = field(default_factory=dict)
    history: Dict[str, Any] = field(default_factory=dict)
    integration_budget: Dict[str, Any] = field(default_factory=dict)
    run_passed: bool = False
    source: Dict[str, Any] = field(default_factory=dict)         # configuration_source, migration, documents (for summary.json)
    #: True when `error` is an exception that run_energy_case does not know (not GribsError / ValueError / RuntimeError):
    #: the command line prints its traceback and ends with code 4 (E_ENERGY_UNEXPECTED_ERROR); the API raises `error` itself.
    unexpected: bool = False

    @property
    def summary(self) -> Dict[str, Any]:
        """The results section plus model_mode and energy_model (a Result has no field for them)."""
        return {**self.results, "model_mode": "energy", "energy_model": self.energy_model}

    def summary_text(self) -> str:
        r, e = self.results, self.energy_model
        lines = ["GRIBS energy-mode run (dedicated minimal output)", f"  scope: {SCOPE}",
                 f"  burning phase stopped by: {r['stop_reason_burning']}; blowdown: {r['stop_reason_blowdown']}",
                 f"  end time [s]: {r['end_time_s']:.9g}; burning ends at [s]: {r['burning_end_time_s']:.9g}",
                 f"  peak pressure [Pa]: {r['peak_pressure_Pa']:.9g} at t = {r['peak_pressure_time_s']:.9g} s",
                 f"  total impulse [N s]: {r['total_impulse_Ns']:.9g}; generated mass [kg]: {r['generated_mass_kg']:.9g}; "
                 f"discharged mass [kg]: {r['discharged_mass_kg']:.9g}",
                 f"  gamma_s range (CEA-type equilibrium isentropic exponent of the provider): {e['gamma_s']['min']:.9g} .. "
                 f"{e['gamma_s']['max']:.9g}",
                 f"  closure residuals (normalised, max over the history; tolerance {CLOSURE_TOLERANCE:g}): "
                 + ", ".join(f"{k} = {v:.3g}" for k, v in r["closure"]["max_normalised"].items()),
                 *[f"  latent heat check at t = 0 ({c['solid_species']} -> {c['liquid_species']}, T_m = "
                   f"{c['melting_temperature_K']:.9g} K): L_calc = {c['L_calc_J_per_kg']:.9g} J/kg, L_input = "
                   f"{c['L_input_J_per_kg']:.9g} J/kg, relative difference {c['relative_difference']:.3g} "
                   f"(tolerance {c['tolerance_rel']:g})" for c in e.get("latent_heat_check", [])],
                 f"  run passed: {self.run_passed}",
                 "  legacy outputs NOT produced in energy mode: " + "; ".join(LEGACY_OUTPUTS_NOT_PRODUCED)]
        return "\n".join(lines) + "\n"

    def write_and_report(self, outdir: Path, names: Dict[str, str], environment: Dict[str, Any]) -> int:
        """The command line side: write the files, print the report and return the exit code (as gribs.cli.main)."""
        if self.error is not None:
            if isinstance(self.error, ConfigurationError):
                print(f"Configuration error:\n{self.error}", file=sys.stderr)
                return EXIT_CODE_CONFIGURATION
            if isinstance(self.error, ThermochemistryError):          # as the legacy command line (exit code 3)
                code = getattr(self.error, "code", None)
                if code:
                    print(f"Thermochemistry backend error [{code}]:\n{self.error}", file=sys.stderr)
                    record = {"code": code, "severity": "error", "message": str(self.error),
                              "context": getattr(self.error, "context", {})}
                    print(json.dumps(record, ensure_ascii=False, sort_keys=True), file=sys.stderr)
                else:
                    print(f"Thermochemistry backend error:\n{self.error}", file=sys.stderr)
                return EXIT_CODE_THERMOCHEMISTRY
            if self.unexpected:
                print("".join(traceback.format_exception(type(self.error), self.error, self.error.__traceback__)),
                      file=sys.stderr, end="")
                print(f"Energy mode: unexpected error [{UNEXPECTED_ERROR_CODE}]: {type(self.error).__name__}: {self.error}",
                      file=sys.stderr)
                record = {"code": UNEXPECTED_ERROR_CODE, "severity": "error", "message": str(self.error),
                          "context": {"exception_type": type(self.error).__name__}}
                print(json.dumps(record, ensure_ascii=False, sort_keys=True), file=sys.stderr)
                return EXIT_CODE_SIMULATION
            print(f"Simulation failed: {type(self.error).__name__}: {self.error}", file=sys.stderr)
            return EXIT_CODE_SIMULATION
        outdir = Path(outdir)
        csvp, txtp = outdir / HISTORY_FILENAME, outdir / names["summary_text_filename"]
        jsnp = outdir / names["summary_json_filename"]
        text = self.summary_text()
        summary_json = {
            "program": self.program, "environment": environment, "model_mode": "energy", **self.source,
            "derived_inputs": self.derived_inputs, "two_phase_model": self.two_phase_model,
            "thermochemistry": self.thermochemistry, "energy_model": self.energy_model, "results": self.results,
            "output_files": {"history_csv": str(csvp), "summary_text": str(txtp), "summary_json": str(jsnp),
                             "figure": None, "legacy_history_csv": None},
        }
        try:
            outdir.mkdir(parents=True, exist_ok=True)
            csvp.write_text(history_csv_text(self.history), encoding="utf-8")
            txtp.write_text(text, encoding="utf-8")
            jsnp.write_text(json.dumps(summary_json, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        except OSError as exc:         # not a closure failure (code 1, files written): the output could not be written
            print(f"Energy mode: the output could not be written [{OUTPUT_WRITE_ERROR_CODE}]: {type(exc).__name__}: {exc}; "
                  "some output files may have been written", file=sys.stderr)
            record = {"code": OUTPUT_WRITE_ERROR_CODE, "severity": "error", "message": str(exc),
                      "context": {"path": str(getattr(exc, "filename", None) or outdir), "errno": exc.errno}}
            print(json.dumps(record, ensure_ascii=False, sort_keys=True), file=sys.stderr)
            return EXIT_CODE_SIMULATION
        print()
        print(text)
        print("Output files")
        for pth in (csvp, txtp, jsnp):
            print(f"  {pth}")
        print("Energy mode: the legacy figure and the legacy diagnostics are not produced "
              "(see energy_model.legacy_outputs_not_produced in the summary).")
        if not self.run_passed:
            print("Run verification failed after the output was written: " + CLOSURE_WARNING_CODE, file=sys.stderr)
            return EXIT_CODE_VERIFICATION_FAILED
        return EXIT_CODE_SUCCESS


def history_csv_text(history: Dict[str, Any]) -> str:
    """The CSV text of a history: the header and one row per node (floats as repr, so a file and an array agree exactly)."""
    columns = list(history)
    n = len(history[columns[0]]) if columns else 0
    rows = [",".join(columns)]
    rows += [",".join(repr(float(history[c][i])) for c in columns) for i in range(n)]
    return "\n".join(rows) + "\n"


def _history_and_closure(run: dict) -> Dict[str, np.ndarray]:
    model = run["model"]
    rows = []
    for burning, sol in ((True, run["sol_burn"]), (False, run["sol_blow"])):
        if sol is None:
            continue
        start = 1 if (not burning and rows) else 0                # the blowdown starts at the last burning node
        for j in range(start, sol.t.size):
            y = sol.y[:, j]
            M, E = float(y[ei.IM]), float(y[ei.IE])
            k = model.terms(y, burning)
            rec, V_free = k["rec"], k["V_free"]
            r = closure_residuals(M=M, M0=run["M0"], m_gen=float(y[ei.IMG]), m_out=float(y[ei.IMO]), B_e=model.B_e(M),
                                  E_cv=E, E0=run["E0"], Q_cum=float(y[ei.IQ]), W_burn_cum=float(y[ei.IWB]),
                                  H_out_cum=float(y[ei.IHO]), V_gas=rec.V_gas_m3, V_condensed=rec.V_condensed_m3,
                                  V_free=V_free)
            rows.append((float(sol.t[j]), M, E, k["p"], k["T"], k["gamma_s"], k["x"], float(y[ei.IRT]), V_free, k["r"],
                         k["mdot_burn"], k["mdot_out"], float(k["nz"]["F"]), float(y[ei.IMG]), float(y[ei.IMO]),
                         float(y[ei.IIMP]), float(y[ei.IWB]), float(y[ei.IHO]), float(y[ei.IQ]), r.R_M, r.R_B, r.R_E, r.R_V,
                         r.R_M_normalised, r.R_B_normalised, r.R_E_normalised, r.R_V_normalised))
    table = np.array(rows, dtype=float)
    return {name: table[:, i] for i, name in enumerate(HISTORY_COLUMNS)}


def run_energy_case(loaded) -> EnergyOutcome:
    """Run the energy-mode configuration `loaded` (a Configuration of gribs.config.loader) and return its EnergyOutcome.
    Nothing is written and nothing is printed.  The errors of the run (GribsError, ValueError, RuntimeError, as the command
    line treats them) are returned in EnergyOutcome.error, not raised.  Any other Exception is returned too, with
    `unexpected` set (the command line prints its traceback and ends with code 4; the API raises it unchanged); KeyboardInterrupt
    and SystemExit are not caught."""
    try:
        return _run_energy_case(loaded)
    except Exception as exc:
        return EnergyOutcome(error=exc, unexpected=True)


def _run_energy_case(loaded) -> EnergyOutcome:
    cfg = loaded.config
    try:
        run = ei.run_model_energy(cfg)
        history = _history_and_closure(run)
    except (GribsError, ValueError, RuntimeError) as exc:
        return EnergyOutcome(error=exc)
    t_all, p_all = history["time_s"], history["pressure_Pa"]
    i_peak = int(np.argmax(p_all))
    maxima = {name: float(np.max(np.abs(history[name + "_normalised"]))) for name in ("R_M", "R_B", "R_E", "R_V")}
    maxima_abs = {name: float(np.max(np.abs(history[name + unit]))) for name, unit in
                  (("R_M", "_kg"), ("R_B", "_kg"), ("R_E", "_J"), ("R_V", "_m3"))}
    closure_ok = all(math.isfinite(v) and v <= CLOSURE_TOLERANCE for v in maxima.values())
    gamma_min, gamma_max = float(run["provenance"]["gamma_s_min"]), float(run["provenance"]["gamma_s_max"])
    sol_last = run["sol_blow"] if run["sol_blow"] is not None else run["sol_burn"]
    burn_end = float(run["sol_burn"].t[-1])
    results = {
        "stop_reason_burning": run["stop_burn"], "stop_reason_blowdown": run["stop_blow"],
        "burning_end_time_s": burn_end, "end_time_s": float(sol_last.t[-1]),
        "peak_pressure_Pa": float(p_all[i_peak]), "peak_pressure_time_s": float(t_all[i_peak]),
        "initial_mass_kg": float(run["M0"]), "final_mass_kg": float(sol_last.y[ei.IM, -1]),
        "generated_mass_kg": float(sol_last.y[ei.IMG, -1]), "discharged_mass_kg": float(sol_last.y[ei.IMO, -1]),
        "total_impulse_Ns": float(sol_last.y[ei.IIMP, -1]), "nodes": int(t_all.size),
        "closure": {"tolerance_normalised": CLOSURE_TOLERANCE, "max_normalised": maxima, "max_absolute": maxima_abs,
                    "passed": closure_ok},
        "warnings": [] if closure_ok else [{"code": CLOSURE_WARNING_CODE, "severity": "error",
                                            "message": "a normalised closure residual of the energy run exceeds the tolerance",
                                            "context": {"max_normalised": maxima, "tolerance": CLOSURE_TOLERANCE}}],
    }
    energy_model = {
        "scope": SCOPE, "contract": "reports/164_P8_A6_STAGE2_CONTRACT.md", "burn_surface_form": cfg.energy.burn_surface.form,
        "burn_surface_interpretation": BURN_SURFACE_INTERPRETATION, "H0_out_definition": H0_OUT_DEFINITION,
        "held_items_adopted": list(HELD_ITEMS_ADOPTED),
        "gamma_s": {"definition": run["provenance"]["gamma_s_definition"], "min": gamma_min, "max": gamma_max},
        "provenance": json.loads(json.dumps(run["provenance"], default=str)),
        "settings": json.loads(json.dumps(asdict(cfg.energy), default=str)),
        "latent_heat_check": run["latent_heat_check"],
        "closure_tolerance_normalised": CLOSURE_TOLERANCE,
        "closure_tolerance_status": CLOSURE_TOLERANCE_STATUS,
        "legacy_outputs_not_produced": list(LEGACY_OUTPUTS_NOT_PRODUCED),
        "legacy_comparison": "none (a run with a condensed volume is a different model from the legacy model and is not compared)",
    }
    notes = list(loaded.migration_notes)
    source = {
        "configuration_source": cfg.configuration_source,
        "configuration_migration": {"applied": True, "notes": notes} if notes else {"applied": False, "notes": []},
        "original_configuration": loaded.document, "resolved_configuration": asdict(cfg),
    }
    return EnergyOutcome(
        error=None, program={"name": PROGRAM_NAME, "version": PROGRAM_VERSION, "schema_version": SCHEMA_VERSION},
        derived_inputs={"propellant_density_kg_m3": cfg.rho_p, "propellant_density_source": "Config.rho_p (B-15)"},
        two_phase_model={"mode": cfg.two_phase_mode, "condensed_volume": cfg.energy.condensed.volume},
        thermochemistry={"backend": cfg.thermo_backend, "provider": run["provenance"]["provider"],
                         "enthalpy_convention": run["provenance"]["enthalpy_convention"]},
        results=results, energy_model=energy_model, history=history,
        integration_budget=dict(run["integration_budget"]), run_passed=closure_ok, source=source)
