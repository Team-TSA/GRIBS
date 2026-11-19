# Changelog

All notable changes to GRIBS are documented in this file.

## [Unreleased]

A-6 stage 1 (package split and Python interface): G0 was approved for revision ea3f4a7 on 2026-11-10. The later items (P5 and A-9 integration, A-6 stage 2 energy mode) have not been through G0. The program version stays `0.6.2-alpha`.

### Breaking

Details of each item are under "Changed".

- P5-1: the `cea_python` cache key no longer contains absolute paths (`cea_module_path` removed; it contains the database digests), so all `cea_python` cache keys change once. `thermo.lib` / `trans.lib` in the working directory or in `CEA_DATA_DIR` are ignored and only the database of the official CEA package is used. A run that used to pick up a different database from the working directory can give different numerical results; with the package database the results are unchanged.
- P5-2: `thermochemistry.cea_python.ions = true` is rejected before any CEA call (`E_IONS_NOT_SUPPORTED`, exit code 2), also by `--validate-config`. Before, the configuration was accepted as valid and a run could end with SIGSEGV.
- P5-6: `cstar_burn_interval_m_s` is now burning-phase p0·At integral / mass discharged up to burnout (before, a copy of `cstar_eff_m_s`). Runs with a blowdown get a larger value (golden g3–g8: about +0.31 % to +0.40 %). `cstar_eff_m_s` and `Isp_s` are unchanged. New key `Isp_discharged_mass_s`.
- A9: `temperature_efficiency` (`eta_T0`) must now lie in (0, 1]. Values above 1 (accepted up to 1.5 before) are rejected with a configuration error (exit code 2). All configurations of the repository use 1.0.

### Changed

- At the end of stage 1, the program code was moved from the single file `GRIBS_v0.6.2-alpha.py` into the package `gribs` under `src/gribs/` (22 modules at that time; stage 2, below, adds the package `gribs.energy` and two modules under `gribs.simulation`). `GRIBS_v0.6.2-alpha.py` remains the program file: it adds `src/` to the import path and imports the moved names back, so existing callers of its names keep working. At the end of stage 1, command-line options, exit codes, outputs and numerical results were unchanged by the split: after each step of the split, the golden cases g1–g8 gave byte-identical time histories, CEA caches and figures, and identical normalized summaries, compared with the v0.6.2-alpha G0 reproduction. The items of P5, A-9 and stage 2 below change outputs where they say so.
- The program file now needs `src/gribs/` next to it; copying the program file alone no longer works.
- New `pyproject.toml` (package `gribs`, placeholder version `0.7.0.dev0`, not used in any output); `pytest.ini` adds `src` to the import path.
- `tools/a5_1_g4a_cea_comparison.py` records the identity of the executed code as the entry file hash combined with a digest of the package files (`body_composite_v1`), instead of the hash of the program file alone.
- P5-1 (F-18, F-26): the `cea_python` backend always uses the database of the official CEA package. `cea` is imported without its automatic database loading (`CEA_SKIP_INIT`, set only for that import), and the package `thermo.lib` is loaded explicitly (`trans.lib` only with `transport`). `thermo.lib` / `trans.lib` in the working directory or in `CEA_DATA_DIR` are no longer read; they are listed in the new warning `W_EXTERNAL_CEA_DATABASE_IGNORED`. Before, such a file was read by `import cea` while the cache key and the provenance named the package database, and a damaged `trans.lib` there ended the process even with transport off. A process whose `cea` was already loaded with another database stops with `E_CEA_DATABASE_ALREADY_INITIALIZED` (`CEADatabaseError`, exit code 3). `thermochemistry.database_identity` records the full SHA-256 (was a 16-digit prefix) and the ignored files. The cache key contains the database digests and no paths (`cea_module_path` removed), so all `cea_python` cache keys change once; with the package database the results are unchanged (golden g2 `p_max_Pa` bit-identical).
- P5-2 (F-19): `thermochemistry.cea_python.ions = true` is rejected as a configuration error `E_IONS_NOT_SUPPORTED` (`ConfigurationError`, exit code 2), also by `--validate-config` and by `gribs.run_simulation`, before any CEA call. Before, the configuration was reported valid and the run ended with SIGSEGV in `cea.EqSolution` (the reactant mixture is built without ions, so CEA found different element lists). The option is still read and still part of the cache key; `ions = false` results, cache keys and provenance are unchanged. Other backends do not read the option and do not reject it.
- P5-3 (F-06): a run whose eroding throat area exceeds the fixed exit area (At > Ae) records the first time of it in the new event `event_times.throat_exceeds_exit_s` (an empty list when it does not occur, like the other events) and adds the warning `W_THROAT_EXCEEDS_EXIT_AREA`. Before, such a run gave no signal although the model is outside its domain there (the nozzle is then computed with the exit as the smallest section). The comparison is made on the radius, Rt > Rt0·√ε, without a tolerance, so a nozzle with ε = 1 and no erosion is never reported (an area comparison would be one ulp off for golden g1). The time is located on the dense solution after the integration. Nothing else changes: the run continues, the warning does not change the run verdict or the exit code, and all computed values are unchanged (golden g1–g8: time histories bit-identical, summaries identical apart from the new empty event). The thrust and mass-flow jump by At/Ae at an unchoking with At > Ae (F-25) is not changed.
- P5-4 (F-02): the erosive-burning burn rate is solved with Newton's method from an upper bound of the root, to a relative step of 1e-13 followed by one more step. Before, a damped fixed-point iteration stopped after 60 steps and returned its last value without notice; for an erosive fraction (r − r0)/r above about 0.3 that value had not converged (relative error about 1e-10 at 0.48, 1e-6 at 0.75 and 0.8 % at 0.999). If the root is not found (at most 100 steps), `burn_rate` raises the new `ErosiveRootNotConvergedError` (`E_EROSIVE_ROOT_NOT_CONVERGED`, a `GribsError`) and returns no value; the command line reports it as a simulation failure with exit code 4. The erosive-burning equation, its coefficients and the 0.8 exponent are unchanged. Runs without erosive burning are unchanged (golden g1–g8 bit-identical); with erosive burning the burn rate changes by the former iteration error.
- P5-5 (F-16): the summary (`results` in `summary.json`, `Result.summary`) has six run-state keys at its top level: `solver_completed`, `mass_consistency_passed` (post-run mass consistency and mass ledger), `nominal_termination_reached` (burnout, then ambient pressure or no blowdown), `model_domain_passed` (no At > Ae and no `E_OUT_OF_VALIDATED_RANGE`), `warnings_present` and `results_usable` (the first four together). A run that stops with an error has no summary and gets no keys. New warnings: `W_NON_NOMINAL_TERMINATION` (any end other than the nominal one, e.g. a burning time limit or a burn ended at ambient pressure or by choking loss), `W_SUBSONIC_BURNING` (the nozzle was subsonic during part of the burning phase) and `W_PROPERTY_TABLE_EXTRAPOLATED` (the chamber pressure left the property-table range; clamped or extrapolated). Leaving the table range and subsonic burning are warnings only: they do not make `model_domain_passed` or `results_usable` false. Before, these states were visible only in `burn_stop_reason`, `blowdown_stop_reason`, `regimes_visited` and `property_extrapolation`, and `run_passed` was true for all of them. `run_passed`, the exit codes and all computed values are unchanged; the warnings and keys do not change the run verdict (golden g1–g8: time histories bit-identical, the old summary values identical; g1 and g2 now also report `W_PROPERTY_TABLE_EXTRAPOLATED` and are `results_usable = true`).
- P5-6 (F-11, F-12): `cstar_burn_interval_m_s` is now what its name says: the burning-phase integral of p0·At divided by the mass discharged up to burnout. Before, it was a copy of `cstar_eff_m_s`, whose denominator is the mass discharged over the whole run including the blowdown. Runs with a blowdown get a larger value (golden g3–g8: +0.31 % to +0.40 %); runs without a blowdown are unchanged. `cstar_eff_m_s` and `legacy_mixed_interval_cstar_m_s` keep the mixed definition, and `Isp_s` keeps the initial-propellant-mass basis (their values are unchanged). New key `Isp_discharged_mass_s`: total impulse per discharged product mass. The definitions are in README, "Effective c* and Isp definitions". Time histories, `run_passed` and the exit codes are unchanged; the tracked reference outputs `results_golden_g*/summary.json` were not regenerated.
- The post-run mass-consistency tolerance is now `max(1e-6, 100 * relative_tolerance)` instead of the fixed 1e-6 (the propellant ledger already used `max(1e-9, 100 * rtol)`). The default `relative_tolerance = 1e-10` gives the same tolerance as before; only a looser `relative_tolerance` widens it (measured: with rtol = 1e-6 golden g2-g6 exceeded the fixed 1e-6). `summary.json` `post_run_consistency` gains `tolerance_formula` and `rtol`.
- `temperature_efficiency` (`eta_T0`) must now lie in (0, 1] (listed under Breaking). Values above 1 (accepted up to 1.5 before) put T0 above the CEA adiabatic (HP equilibrium) temperature and are rejected with a configuration error (exit code 2); use `eta_T0 <= 1`, 1.0 being the standard model. All configurations in this repository use 1.0, so the golden cases are unaffected. `W_THERMAL_EFFICIENCY_PARTIAL` is unchanged for `eta_T0 < 1`.

### Added

- Python interface: `gribs.run_simulation`, `gribs.Result`, `gribs.__version__` (README, "Python API", with its limits).
- Tests: guard tests for the package split (patch targets, re-exports, module ownership, thermochemistry backend state, output version) and acceptance tests of the Python interface.
- A run that crosses `p_fit_min` / `p_fit_max` (`thermochemistry.outside_range_policy = "clamp"`) now also reports the warning `W_PROPERTY_TABLE_RANGE_CROSSED` in `summary.json` (`warnings`, with the restart records in `context`); previously only the stderr line existed. The summary schema and CSV format are unchanged.
- `W_PROPERTY_TABLE_EXTRAPOLATED` (P5-5) is the single producer of this warning; the A-9 version of it is not kept. Two items of A-9 are added to it: the `context` gains `points_total` (number of sampled points) and `basis` (`sampled_history_all_phases`: the decision covers the whole sampled history, burning and blowdown), and the message names the configuration keys that widen the table (`thermochemistry.pressure_range.minimum_Pa` / `maximum_Pa`). Measured: golden g2 spends 84 % of the burn above the 8 MPa table limit; g3-g6 stay inside. No number changes. A-9 is not released, so no external consumer of its other context key names exists.
  Context keys after the integration: `policy`, `points`, `fraction`, `property_range_min_Pa`, `property_range_max_Pa`, `pressure_min_Pa`, `pressure_max_Pa` (P5-5) and `points_total`, `basis` (A-9). The other A-9 key names (`property_policy`, `points_outside`, `fraction_outside`, `table_min_Pa`, `table_max_Pa`, `history_min_Pa`, `history_max_Pa`) are equivalent to the P5-5 keys and are not used.
- `append_migration_model_change_warning(summary, migration_notes)`: the logic that adds `W_MIGRATION_MODEL_CHANGE` is now a function in `gribs.simulation.postprocess`; `gribs.cli.main` and `gribs.run_simulation` both call it (their `summary.json` is unchanged) and a caller that runs `summarize()` itself on a migrated configuration can pass `LoadedConfiguration.migration_notes` to get the same warning.
- New tests: `tests/test_v043_anchor.py` (the current program reproduces the published v0.4.3-alpha, taken from the git tag, for a condensed AP/Al/C7H8 case: all total quantities equal to 1e-9; measured difference 0) and `tests/test_wide_table_cases.py` with the non-frozen wide-table variants in `cases/wide_table/` (golden g2 / g5 with a larger `maximum_Pa`, runs entirely inside the table).
- A-6 stage 2: an opt-in energy mode. A configuration with an `energy_model` section is run with a chamber energy balance (state: total energy and element inventories; pressure and temperature recovered from them) instead of the pressure law of the chamber-property table. Without that section nothing changes: legacy runs give the same files and values, the same exit codes, and no new summary key. The program version stays `0.6.2-alpha` and the schema version is unchanged; the section is optional, and an energy run adds keys (below) that a legacy run does not have.
- Scope: energy mode is verified on synthetic fixtures only; the number of cases qualified with a real propellant is 0. Results are not compared with measured data, nor (for a condensed volume) with the legacy model. See "Energy mode" in `MODEL_LIMITATIONS.md`.
- New configuration keys, all under `energy_model`: `initial_fill` (`temperature_K`, `composition`, `element_tolerance`), `wall` (`adiabatic`, `heat_rate_W`), `burn_surface` (`form`, `grain_far_temperature_K`, `solid_enthalpy`), `condensed` (`volume`, `phases`, `transitions`), `recovery` (`phase_fraction_tolerance`, `root_tolerance`, `plateau_edge_tolerance`, `temperature_search_K`, `start_temperatures_K`) and `propellant_formulas` (explicit formulas of the reactants; required to run). There is no key for the density of the solid: `Config.rho_p` is used.
- The CEA TP provider of the energy mode solves each state on a new `EqSolution`: a reused one failed at a solid-to-liquid change of Al2O3 and, after a failure, returned a `gamma_s` below 1 with the other values unchanged (report 198). The TP cache schema is now 3; earlier TP cache files are not read.
- Known limits of energy mode (held items H-19 to H-22 in contract 164 section 11): the `gamma_s` of states with condensed Al2O3 has not been checked against an independent value and a run through the melting plateau with the CEA provider has not been tested (H-19); the remaining history dependence of the CEA solve (H-20); the tolerance of the t = 0 latent-heat check and the status of `L_calc` as an approximation (H-22). H-21 records, for a separate U-ID, that the legacy chamber-property table solves its pressure grid on one reused `EqSolution`; it is outside A-6. `E_cv` is treated as the internal energy (held item H-05).
- `gamma_s` of the nozzle flow is the equilibrium value of the CEA TP solution (and the configured `reference_gamma` for `reference_constant`); a value that is not finite or not above 1 stops the run.
- New outputs of an energy run: `energy_time_history.csv` (columns with units, and the closure residuals), and in `summary.json` the top-level keys `model_mode = "energy"` and `energy_model` (scope, interpretation, definition of `H0_out`, held items adopted, `gamma_s` range, provenance, settings, closure threshold and `legacy_outputs_not_produced`). The legacy figure and diagnostics are not produced, and the command line and the summary say so. `gribs.run_simulation` runs an energy configuration and returns a `Result` whose `summary` also holds `model_mode` and `energy_model`.
- Rejected in energy mode (the main items; the full list is item 9 of "Energy mode" in `MODEL_LIMITATIONS.md`): an igniter (`E_ENERGY_MODE_IGNITER_REJECTED`), `equilibrium_pressure_diagnostic` other than `off`, `temperature_efficiency` other than 1, an initial fill that is not the products of the main propellant (an air fill), and `--selftest`, `--dump-thermo`, `--no-cache` or `--allow-verification-failure` with an `energy_model` section.
- The energy-mode closure threshold (normalised residuals within `1e-6`) is provisional (held item H-09); it is a pass/fail threshold, not a claim of accuracy.
- `reference_constant` in energy mode uses the test convention `h(T_ref) = h_s(T_grain,inf)`; it is for tests only.
- Exit codes of energy mode for failures outside the model: a failure to write the TP cache stops the run with code 3 (`E_THERMO_CACHE_WRITE_FAILED`; the cache directory is created first, as for the legacy cache); an unexpected exception ends with code 4 (`E_ENERGY_UNEXPECTED_ERROR`, traceback on stderr; the Python API raises the exception unchanged); a failure to write the output files ends with code 4 (`E_ENERGY_OUTPUT_WRITE_FAILED`), not code 1, which stays the closure failure; and the run is limited to 300 s of wall time over the whole run, as the legacy model (`E_INTEGRATION_BUDGET_EXCEEDED`, code 4).

### Documentation

- README: package layout, quick start, Python API; G0 approval of the v0.6.2-alpha release commit `3c0c98c`; G0 records in Release assets and in `docs/G0/`; the small margin of the total-mass-consistency check, the CI environment and the cache-key path (GitHub Issue #10, items 1–3; the issue stays open). The v0.6.2-alpha section no longer says "The README states that ...".
- MODEL_LIMITATIONS.md: G4b has not been done because no usable public data was found.
- MODEL_LIMITATIONS.md item 9 and README, "Known limitations": the At > Ae event `event_times.throat_exceeds_exit_s` and the warning `W_THROAT_EXCEEDS_EXIT_AREA` (P5-3).
- README, "Run-state keys": the six run-state keys of P5-5. `summary.txt` prints them in a new `[ run state ]` section and lists every warning code of the run in a new `[ warning codes ]` section (display only; `summary.json` and all computed values are unchanged).
- MODEL_LIMITATIONS.md item 5 describes the erosive-burning model as implemented (F-01, documentation only): one aft-end total mass flux, applied to the whole burning surface; not a port-average model; the unit of alpha.
- `docs/G0/G0_APPROVAL_3c0c98c_2026-11-07.txt`: the G0 approval record of v0.6.2-alpha with an English summary.
### Notes

- Warning-code vocabulary: `W_MIGRATION_MODEL_CHANGE` is the implemented code for a configuration migration that changes the physical model (the command line and `gribs.run_simulation` add it; a caller that runs `summarize()` itself uses `append_migration_model_change_warning`). `W_MASS_BASIS_CONFIG_MIGRATED` (improvement plan candidate) is not implemented and is covered by `W_MIGRATION_MODEL_CHANGE`; `W_DEPRECATED_MOLECULAR_WEIGHT_KEY` (review proposal) was not adopted: `gas_phase_M` is a frozen compatibility token accepted without a warning.


## [v0.6.2-alpha] - 2026-11-07

Patch release of v0.6.0-alpha "Keystone". Numerical results are unchanged; a model-applicability warning is added and the README G0 text is corrected. Status: unchanged from v0.6.1-alpha (physical validity not established; G4a not passed; G4b and G5 not complete). The G0 decision for the v0.6.2-alpha release commit is stated in the release notes.

### Release identity

- `PROGRAM_VERSION` is `0.6.2-alpha`; `--version` reports `GRIBS 0.6.2-alpha`.
- The program file is renamed from `GRIBS_v0.6.1-alpha.py` to `GRIBS_v0.6.2-alpha.py`. Scripts that call the old file name must be updated.
- The configuration schema is unchanged: `SCHEMA_VERSION` stays `0.5.0-alpha`.

### Added

- Warning `W_HIGH_CONDENSATE_G4B_UNVALIDATED` (severity `warning`): reported when `validation.model_status` is `approximate_with_limitations` and the thermochemistry table reaches `Yg_min <= 0.80` (the high-condensate band). It states that reproduction of reported or experimental data (G4b) has not been done for this band. It is reported together with `W_HIGH_CONDENSATE_MODEL_LIMITATION`, which continues to describe the unmodelled two-phase physics. `model_status`, its thresholds (0.99 / 0.65), the run verdict, the exit codes, the numerical results, the configuration schema and the thermochemistry cache are unchanged. Out-of-scope tables (`Yg_min <= 0.65`) keep stopping with `E_OUT_OF_VALIDATED_RANGE` and do not get the new warning. Golden cases g4 and g5 now report it; their values are unchanged.

### Tests

- The two model-applicability warning tests expect the new warning in their exact warning list. New boundary test at `Yg_min` = 0.65 / 0.80 / 0.8001: the warning is present at 0.80, absent at 0.8001, absent for `out_of_scope` at 0.65, and never duplicated.

### Documentation

- README "Validation status": G0 is described as approved for the v0.6.1-alpha release commit `cf91f97` (as well as `4e7896a`); the earlier text said that this approval was still to be decided.
- README "Known limitations" and MODEL_LIMITATIONS.md describe the condensed-phase scope bands and the two high-condensate warnings.

## [v0.6.1-alpha] - 2026-11-06

Patch release of v0.6.0-alpha "Keystone". The simulation code is unchanged
apart from the version identifiers and one source comment; the default
configuration, the golden set and the documentation change. Status:
unchanged from v0.6.0-alpha (physical validity not established; G4a not
passed; G4b and G5 not complete). G0 for the v0.6.1-alpha release commit is
decided by the maintainer separately; see the README section
"Validation status".

### Release identity

- `PROGRAM_VERSION` is `0.6.1-alpha`; `--version` reports `GRIBS 0.6.1-alpha`.
- The program file is renamed from `GRIBS_v0.6.0-alpha.py` to
  `GRIBS_v0.6.1-alpha.py`. Scripts that call the old file name must be
  updated.
- The configuration schema is unchanged: `SCHEMA_VERSION` stays
  `0.5.0-alpha`.
- A source comment no longer refers to `VALIDATION_REPORT.md`, which does not
  exist; it refers to the README section "Validation status".

### Changed

- The default configuration `gribs_config.json` uses graphite `C(gr)` instead of `C`. In the NASA CEA database `C` is gaseous atomic carbon (about +59.67 MJ/kg), which made the reactant enthalpy strongly positive and the default results non-physical (CEA HP at 30 bar: chamber temperature 3678 K and `c*` 1901 m/s with `C`; 2138 K and 1330 m/s with `C(gr)`). Results of the default configuration change accordingly. With `C(gr)` the minimum gas mass fraction of the default configuration is about 0.987, below the in-scope threshold 0.99, so its `validation.model_status` is `approximate_with_limitations` (previously `in_scope_not_validated`) and the run reports the warning `W_HIGH_CONDENSATE_MODEL_LIMITATION` in addition to `W_TWO_PHASE_NOZZLE_APPROXIMATION`. `examples/example_config.json` already used `C(gr)`.
- The notes of `gribs_config.json` are updated: they no longer describe the file as a v0.4.3-alpha configuration, no longer say that the backend is left empty, and no longer refer to `VALIDATION_REPORT.md`, which does not exist.

### Documentation

- README, MODEL_LIMITATIONS.md: the outside-range setting is named by its JSON key `thermochemistry.outside_range_policy` (previously the internal name `property_policy` was used, which is not accepted in a configuration file); the README "Configuration" section describes it and lists the `two_phase_model` section.
- README, MODEL_LIMITATIONS.md: GRIBS does not check whether a composition is plausible. The shape of the pressure and thrust histories is set mainly by the grain geometry and the burn law; with these fixed, normalized pressure histories of different compositions nearly overlap. The README explains the CEA species names `C` / `C(gr)` and `AL` / `AL(cr)`.
- README: golden cases g3 and g6 use gaseous `C`; they are kept unchanged as numerical regression references.
- README "Validation status": G0 is described for the v0.6.0-alpha release commit `4e7896a` (797 tests), for which it was approved after the v0.6.0-alpha entry below was written; the v0.6.1-alpha changes are stated as not covered by that approval.

### Tests

- New golden cases g7 and g8 (`cases/golden/g7_low_condensed_graphite.json`, `g8_cd_nozzle_blowdown_graphite.json`): g3 and g6 with graphite `C(gr)` instead of gaseous `C`, otherwise identical (carbon density 2200 kg/m³ as in g3 and g6). Their values are frozen in `golden_values.json`; the entries of g1-g6 are unchanged. They have no recorded output directory.
- The golden regression test also compares `two_phase_mode` and `two_phase_active` with the frozen values (previously recorded but not checked).
- `test_cd_nozzle_transition_classification` reads a dedicated fixture, `tests/fixtures/cd_nozzle_transition_eps4.json` (the v0.6.0-alpha default configuration with Ae/At = 4), instead of the default `gribs_config.json`. Its frozen jump values are unchanged, and it no longer depends on the default configuration.

## [v0.6.0-alpha] "Keystone" - 2026-11-05

**Consistent Thermochemistry Foundation.** First public release after
v0.4.3-alpha. v0.5.0-alpha (below) was an internal Phase A baseline and was
not published; this release contains it, the A-5.1-5 to A-7 changes listed in
`CHANGELOG_v0.5.0-alpha_A5.1-5_to_A7.md`, the A-8 golden regression cases, and
the changes in this section.

Status: alpha research software. Physical validity has not been established.
Gate G4a (comparison with the CEA reference `c*`) has not passed and is
documented as a known limitation with its definition and tolerance unchanged;
G4b and G5 are not complete. G0 (reproducibility) was approved for the code
basis `6ec2e6b` only, and G1 (CEA semantics) within a limited scope. See the
README sections "Validation status" and "Known limitations". Stopping at the
property-table range (M3), the safe output mode, and the Phase B work (B1,
B2-C, B2-D) are not included.

### Release identity

- `PROGRAM_VERSION` is `0.6.0-alpha`; `--version` reports `GRIBS 0.6.0-alpha`,
  and `summary.json`, `summary.txt` and the figure note show the new program
  version.
- The program file is renamed from `GRIBS_v0.5.0-alpha.py` to
  `GRIBS_v0.6.0-alpha.py`. Scripts that call the old file name must be
  updated.
- The configuration schema is unchanged: `SCHEMA_VERSION` stays
  `0.5.0-alpha`, so existing v0.5.0-alpha configurations, the golden and
  reference cases, and the chamber-property cache keys are unaffected.

### Fixed (late Phase A)

- The ODE integration no longer stalls at the property-table limits. With `property_policy = "clamp"` (table backends) dp/dt jumps at `p_fit_min` / `p_fit_max`; LSODA could shrink its step to ~1e-10 s there and, depending on its history, never finish (observed for g2 with a 1e-13 perturbation of the gamma_s table). The integration now ends a segment at each crossing of a table limit and continues with a fresh solver from the crossing state. Runs that never cross a limit, and the reference_constant backend, are integrated exactly as before. A run that restarts writes one diagnostic line `INTEGRATION_TABLE_RESTARTS {...}` to stderr; the summary and CSV formats are unchanged.
- Safety limits for one run: 2,000,000 right-hand-side evaluations, 300 s of wall time and 10,000 restarts at the table limits. Exceeding one stops the run with `E_INTEGRATION_BUDGET_EXCEEDED` (exit code 4) instead of running without end.
- Run outputs are kept when the comparison-only CEA theoretical diagnostics
  fail or are incomplete. Previously a failed diagnostic raised `KeyError`
  after the simulation and before any output file was written. The status is
  reported as `cea_theoretical_status` (`not_requested` / `ok` / `failed` /
  `incomplete` / `invalid`); no CEA value is substituted, and the simulation
  results and the run verdict are not affected.
- A gas-mass-fraction value whose `__float__` raises is reported as a
  structured phase-fraction error instead of an unhandled exception.
- Non-finite numbers are rejected before any calculation:
  - configuration files: the non-standard JSON constants `NaN`, `Infinity`
    and `-Infinity` are rejected (exit code 2);
  - any non-finite number in a configuration document, including values such
    as `1e999` that overflow to infinity and values in disabled sections, is
    rejected with its path;
  - `validate()` rejects non-finite values in `Config` fields, including
    NumPy scalars, 0-d arrays and `Decimal` values; finite values are never
    converted.
- The chamber-property cache is written atomically: a temporary file in the
  cache directory, attributes of an existing cache copied first, content
  synced, a new cache published without replacing an existing name, and the
  directory synced. A write failure stops the run with exit code 3 and the
  structured error `E_THERMO_CACHE_WRITE_FAILED`; there is no fallback to an
  in-place write. If syncing the directory is not supported, the run
  continues with `W_THERMO_CACHE_DURABILITY_NOT_GUARANTEED`. Operating
  assumptions, compatibility notes and untested conditions are in the README.

### Changed (tests)

- The golden regression test checks the three residual metrics (propellant mass balance, total mass consistency, maximum condensed mass fraction) against absolute bounds equal to the program's verification thresholds (mass-ledger propellant gate max(1e-9, 100 x rtol); post-run mass tolerance 1e-6; phase-fraction rounding tolerance 1e-12 around the golden value) instead of a relative 1e-9 comparison with the frozen values. These residuals are round-off level: a last-digit difference of the CEA table (for example on another CPU) changed them by tens of percent and failed CI on Python 3.13 while every physical metric stayed within its tolerance. Physical and timing metrics and the golden values are unchanged.

### Changed

- `validation.model_status` no longer reports `validated_in_scope`, because no
  external or experimental validation is linked yet (A-5.1 G4a/G4b open).
  The rules are checked in order; the first match wins:
  `single_phase_legacy` -> `regression_comparison_only`;
  two-phase not active -> `in_scope_not_validated`;
  `Yg_min <= 0.65` -> `out_of_scope`; `Yg_min <= 0.99` ->
  `approximate_with_limitations`; otherwise `in_scope_not_validated`.
  The thresholds and the `out_of_scope` / `approximate_with_limitations`
  results are unchanged.
- `validation.semantic_status` no longer reports `pass` from the cache schema
  number alone. A schema-4 CEA table without a linked semantic probe now
  reports `schema_ok_semantic_probe_not_linked`. Older schemas still report
  `fail`. `semantic_status_basis` and `model_status_basis` state what each
  status is based on.
- Corrected the wording of the `M` description, the `cea_py_molecular_weight`
  parameter text, and the `total_MW` rejection message. They now say M_CEA
  (`EqSolution.M`, kg-total per kmol-gas), not a gas-phase mean molecular
  weight. The token `gas_phase_M`, the rejection of `total_MW`
  (`ConfigurationError`), and the cache key are unchanged.
- `diagnostics.density_note` now describes `EqSolution.density` as total
  mixture mass per gas-phase volume (condensed volume neglected).
- `validation.semantic_status` now depends on the backend kind (cleanup P3).
  The thresholds are unchanged.
  - Backends that do not use CEA report `not_applicable_no_cea_table`. Before,
    `reference_constant` reported `fail` and the synthetic self-test backend
    reported `schema_ok_semantic_probe_not_linked`.
  - CEA table backends keep `schema_ok_semantic_probe_not_linked` for schema
    >= 4 and `fail` for schema < 4.
  - A missing or non-numeric `cache_schema` now reports `fail` with an
    "unknown cache schema" basis, instead of raising.
  - `status_previous_tokens.semantic_status` still reports what the former
    rule gave (`fail` for an unknown schema).
- `semantic_status_basis` and `model_status_basis` now describe the branch
  actually taken, instead of one fixed sentence.
- `diagnostics.sonic_velocity_status` is now classified from the stored
  column with `UNAVAILABLE` / `AVAILABLE` / `INVALID`; `PARTIAL` is replaced
  by `INVALID`. A value that cannot be converted to a float (for example an
  integer beyond the float range) is reported as `INVALID`. The
  classification never raises.
- In the text summary, the table pressure and Yg ranges are labelled
  "[not validated]".
- Thermochemistry cache files are never written or rewritten by these
  changes. The cache key, schema 4 and stored contents are unchanged.
  A reused cache keeps the diagnostic text its writer stored.
- The self-tests now build every comparison backend with an explicit
  two-phase mode, independent of the mode of the configuration they were
  started with (cleanup P4).
  - The legacy side of the homogeneous/legacy comparison uses a backend built
    in `single_phase_legacy` mode, so its legacy branch is exercised. That
    branch's recursion (fixed by P2) previously went undetected.
  - Simulation results and output fields are unchanged.

### Added

- Added machine-readable `W_HIGH_CONDENSATE_MODEL_LIMITATION` and
  `W_TWO_PHASE_NOZZLE_APPROXIMATION` diagnostics for active two-phase runs.
- Included the warnings and their numerical context in the text and JSON
  summaries. The warning-only cases do not change the physical model or run
  verdict.
- Added `W_THERMAL_EFFICIENCY_PARTIAL` when a CEA backend applies `eta_T0 != 1`:
  it explains that CEA HP properties are not re-solved after temperature scaling,
  so direct density, sonic-velocity, and CEA-performance comparisons are not
  like-for-like. This is diagnostic only; it does not change the run verdict.
- Added `E_OUT_OF_VALIDATED_RANGE` as a pre-integration gate when an active
  homogeneous-equilibrium property table has `Yg_min <= 0.65`. The error code
  and context are emitted to stderr; self-test and thermo-dump diagnostic modes
  remain available.
- Added `E_PHASE_FRACTION_VALIDATION_INVALID`, distinct from
  `E_OUT_OF_VALIDATED_RANGE`. It is raised when an active two-phase backend
  reports a gas mass fraction range that is missing, non-numeric, non-finite,
  or outside `0 < Yg_min <= Yg_max <= 1`. The pre-integration gate and direct
  `summarize()` calls both raise it; no default range is assumed. Legacy mode
  and backends without two-phase support do not use the range. The
  `reference_constant` backend now states its analytic range explicitly
  (`Yg = 1`, `Yc = 0`, matching `phase_fractions()`).
- Added `E_OLD_THERMO_CACHE` when the exact matching schema-3 cache is found
  under its legacy key, preventing a silent rebuild after the schema-4 key
  change. The error includes structured path/schema context; unrelated old
  cache entries do not block a run, and explicit rebuild/cache-disabled paths
  remain available.
- Added `E_THERMO_MASS_BASIS_INVARIANT_FAILED` checks for pressure-grid shape,
  finiteness, positivity and ordering; `R_eff*M_CEA = Ru`,
  `M_CEA*n_gas_per_kg_total = 1`, `psi = R_eff*T_used`, and phase-fraction
  sum and physical ranges. CEA tables now persist the semantics tag
  and independent `n` column; compatible pre-tag schema-4 caches remain readable
  after the numerical invariants pass.
- Added structured `E_GAMMA_S_BELOW_UNITY` when generated, cached, or interpolated
  `gamma_s <= 1`; it reports invalid pressure interval(s), gamma range, and cause
  instead of clipping or returning a generic property error.
- Added the Q18 mass-ledger CLI contract: one `W_MASS_LEDGER_FAILED` warning,
  exit code 5 for a ledger-only failure, and `--allow-verification-failure`
  limited to that case. Other verification failures remain exit code 1, and
  normal-run diagnostics are written before the exit.
- `validation.status_previous_tokens` (plus `MODEL_STATUS_PREVIOUS_TOKENS` and
  `SEMANTIC_STATUS_PREVIOUS_TOKENS`) maps each new status token to the token
  that earlier versions reported, so existing readers can compare.
- Correctly named `property_validation` keys
  `equilibrium_molecular_weight_min_kg_kmol` / `_max_kg_kmol`. The former
  `gas_phase_molecular_weight_*` keys keep the same values and are listed in
  `deprecated_keys`.
- `diagnostics.sonic_velocity_status` (`AVAILABLE` / `UNAVAILABLE` /
  `PARTIAL`) and `sonic_velocity_note`; the stored sonic column is unchanged.
- `n_gas_per_kg_total_definition` (kmol-gas per kg of total products,
  `M_CEA = 1/n`) and a unit statement in `molecular_weight_definition`.
  Nothing is renamed.
- `semantic_definitions` in the CEA backend metadata, `summary.json` and the
  thermo dump. It gives the running code's definitions of density, `n` and
  the sonic-velocity status, and observable facts about the table in use:
  `table_source`, `table_fields_missing`, and
  `table_density_note_matches_current`. A cache's creation date is not
  inferred.
- `validation.table_pressure_range_Pa`, `validation.table_Yg_range` and
  `validated_range_basis`.
  - `table_*` are filled only for a backend with a chamber-property table and
    a complete, finite Yg range; otherwise they are null.
  - `validated_range_basis` states what `validated_pressure_range_Pa` /
    `validated_Yg_range` are: the table range, or the stated analytic range
    of `reference_constant` (`[1, 1]`), or defaults. It never claims a
    validated range.
  - The values of `validated_pressure_range_Pa` and `validated_Yg_range` are
    unchanged.

### Fixed

- `--backend` now validates the configuration that will actually run:
  - The legacy/homogeneous exclusion, the sections each backend requires, and the cache
    settings apply to the overriding backend.
  - `--validate-config` reports invalid combinations (exit code 2).
  - The file is first read on its own. Errors found in that step are reported first, with
    the same message as without the override.
  - Value checks (for example, `p_a >= p_fit_min`) run afterwards on the effective
    configuration. With several problems, the order of validation decides which error
    appears first.
  - Overriding to `reference_constant` keeps the documented Phase A-5 defaults when the
    section is absent. A partial section in a file written for another backend is still
    completed with the same defaults, as before.
  - In the tested conditions, numerical results before and after the change were
    identical. See `verification/cli/CLI_BACKEND_OVERRIDE_ADR.md`.

### Validation

- Added unit coverage for warning scope, duplicate suppression, and the
  round-off threshold used to avoid a false nozzle warning at machine precision.
- Added phase-fraction validation tests:
  - missing, non-numeric, NaN/inf, and out-of-range Yg on an active two-phase
    backend, through the gate, `main()`, and direct `summarize()`
  - legacy mode and backends without two-phase support are not affected
  - genuine out-of-scope classification is unchanged
  - the reference backend's analytic phase fractions
- Added targeted coverage for matching legacy-cache rejection, structured error
  context, unrelated-cache false-positive avoidance, and rebuild/cache-disabled
  escape paths.
- Added corrupted-cache probes for the mass-basis identities, semantics-version
  mismatch, structured stderr output, and compatibility with earlier schema-4
  caches that lack the additive semantics-version and `n_gas_per_kg_total` fields.
- Added `W_THERMAL_EFFICIENCY_PARTIAL` scope, direction, idempotence, non-verdict,
  and summary-text coverage.
- Added generated/cache/interpolated probes for `gamma_s <= 1`, disjoint invalid
  pressure intervals, structured stderr context, and no-silent-clipping behavior.
- Added Q18 tests for single-warning idempotence, ledger-only exit code 5,
  explicit ledger-only downgrade, non-overridable gates, and diagnostic files
  written before exit.
- Added CLI `--backend` override tests: two-stage validation, the same first error with
  or without an override, reference defaults for an absent or partial section, the
  backend value in the file, and migrated v0.3.0/v0.4.3/v0.4.4 documents.
- Added semantic-label regression tests (`tests/test_cleanup_semantic_labels.py`)
  and G3-b falsification tests (`tests/test_a5_1_g3b_falsification.py`).
  They cover status tokens and previous-token mapping, unchanged `total_MW`
  rejection, deprecated aliases, and the density/sonic/`n` definitions.
  Evidence tool: `tools/a5_1_cea_semantic_probe.py` (independent CEA probe for
  G1). It produces evidence only and does not approve G1.
- Added `tests/test_cleanup_p3_compat.py`, which covers:
  - per-backend semantic status and its basis
  - unknown cache schemas
  - sonic-column classification, including non-convertible entries and a
    reused cache with altered fields
  - checks that the cache file's bytes and mtime are unchanged
  - range fields with incomplete Yg ranges and the reference analytic range
  - summary and thermo-dump outputs
- Updated the semantic-label tests.
- Added a direct-`summarize()` check for legacy and two-phase-inactive runs.
- Added `tests/test_selftest_backend_modes.py`.
- Added `tests/test_n6_homogeneous_equals_legacy_totals.py` (10/09 N-6).
  - On g4 and g5, the homogeneous model equals `single_phase_legacy` in total
    quantities.
  - Each mode generates its own CEA table (the cache key includes the mode).
  - The test checks that the two stored tables are identical apart from
    `cache_key`.
  - This depends on the P2 recursion fix.
- The golden regression test now writes its outputs under pytest's
  `tmp_path` (cleanup P5).
  - For the `cea_python` cases (g2-g6), the CEA cache directory is set
    explicitly to `tmp_path/results_golden_<id>/thermo_cache`. The
    `reference_constant` case (g1) has no cache.
  - The tracked `results_golden_*` directories and the frozen inputs are no
    longer written by test runs.
  - Golden values are unchanged.

## [v0.5.0-alpha] - 2026-10-02

### Changed

- Aligned the program identity, configuration schema, command examples, and
  generated result metadata with `GRIBS_v0.5.0-alpha.py`.
- Promoted `gribs_config.json` to the native `0.5.0-alpha` schema.
- Added the explicit homogeneous-equilibrium two-phase configuration to the
  native `0.5.0-alpha` configuration.
- Added an auditable `0.4.4-alpha` to `0.5.0-alpha` configuration migration
  that preserves all physical and numerical settings.
- Preserved staged migration from `0.3.0-alpha` through `0.4.3-alpha` and
  `0.4.4-alpha` to `0.5.0-alpha`.
- No calculation model, physical equation, or numerical integration setting
  was changed by the A-0 version-alignment work.

### Validation

- Python syntax compilation completed successfully under Python 3.14.4.
- `--version` reports `GRIBS 0.5.0-alpha`.
- The native `0.5.0-alpha` configuration validates without migration
  warnings.
- A legacy `0.4.3-alpha` configuration migrates successfully through
  `0.4.4-alpha` to `0.5.0-alpha`.
- The standard `Ae/At = 1` self-test configuration passes all 15 checks.
- The complete internal-ballistics simulation and post-run conservation
  checks pass.
- The pre-A-0 and post-A-0 time histories have identical dimensions of
  8,690 rows by 37 columns.
- All 295,460 compared numerical CSV values are exactly identical, with
  maximum absolute and relative differences of zero.
- All common numerical values under `summary.json` `results` agree within a
  relative tolerance of `1e-12`.

### Known limitations

- The pre-existing C-D nozzle T4 failure at `Ae/At = 4` remains:
  mass-flow continuity is maintained, but the measured thrust jump is
  `1.04e+00`.
- The T4 issue is assigned to Phase A, A-2. It was reproduced before and
  after A-0 and was not introduced or modified by the version-alignment work.
- Alpha release; independent validation remains necessary before engineering
  or safety-critical use.

## [v0.4.3-alpha] - 2026-09-20

### Changed

- Updated the public program identifier to `GRIBS v0.4.3-alpha`.
- Updated the configuration schema identifier to `0.4.3-alpha`.
- Updated command examples and user-facing version strings to use `GRIBS_v0.4.3-alpha.py`.
- Corrected comments and documentation so that they accurately match the implemented behavior.
- Corrected the solver-facing thermochemistry interface documentation to state that property evaluations made by the ODE must be side-effect free and inexpensive.
- Clarified that `cea_legacy_executable` is a transitional compatibility backend retained for regression comparison.
- No changes were made to the core calculation logic relative to the supplied v0.4.2-alpha implementation.

### Thermochemistry changes relative to the v0.3 series

- Added the recommended `cea_python` backend based on the official NASA CEA Python package.
- Renamed the former `cea2` backend to `cea_legacy_executable`.
- Removed the manual `legacy_fit` gas-property backend.
- Added explicit, mandatory backend selection with no default and no silent fallback.
- Added backend-specific cache identities and consistency checks.
- Added thermochemistry package, library, database, and cache provenance.
- Added physical validation of pressure-indexed chamber-property tables.
- Added optional theoretical CEA rocket-performance diagnostics for comparison only.

### Configuration

- Updated `schema_version` to `0.4.3-alpha`.
- Standardized bulk solid-propellant density calculation as:

  ```text
  rho_p = packing_fraction / sum_i(w_i / rho_i)
  ```

- Clarified that CEA combustion-product gas density is never used as solid propellant density.
- Added strict validation of required fields, data types, numeric ranges, backend-specific options, and unknown keys.
- Added automatic migration support for the recognized v0.3 schema while reporting every migration action.

### Command-line interface

- Added `--backend` to override the configured backend for one run.
- Added `--validate-config` to validate the JSON configuration without initializing thermochemistry.
- Added `--dump-thermo` to export the pressure-indexed chamber-property table.
- Added `--no-cache` to force rebuilding of the chamber-property table.
- Added `--version` to report the installed program version.

### Outputs and diagnostics

- Expanded `summary.json` with program, Python, platform, dependency, configuration, cache, and thermochemistry provenance.
- Added explicit separation between GRIBS effective performance values and theoretical CEA diagnostic values.
- Improved the output figure layout and automatic fitting of the key-results panel.
- Renamed the default cache output concept from the v0.3 CEA-specific cache to a backend-aware thermochemistry cache.

### Compatibility

- The v0.4 configuration schema is not a drop-in replacement for the v0.3 schema.
- Normal v0.4 configurations must select either `cea_python` or `cea_legacy_executable`.
- The old backend name `cea2` is recognized only by the explicit v0.3 migration path.
- The removed backend name `legacy_fit` is rejected and is never replaced automatically.
- The established unsteady chamber-pressure, grain-regression, burn-rate, nozzle-flow, throat-erosion, and blowdown models are retained.

### Validation

- Python syntax compilation completed successfully.
- JSON syntax validation completed successfully.
- `--version` reports `GRIBS 0.4.3-alpha`.
- `--validate-config` accepts the supplied `0.4.3-alpha` configuration.
- Full numerical results still require validation in an environment with the selected NASA CEA backend and all Python dependencies installed.

### Known limitations

- Alpha release, not certified for safety-critical use.
- Chemical equilibrium does not represent finite-rate chemistry.
- Burn-rate and loss coefficients remain formulation-specific and hardware-specific empirical inputs.
- The constituent-density mixing rule assumes additive constituent volumes before packing correction.
- External executable compatibility must be validated by each user selecting `cea_legacy_executable`.

## [v0.3.1-alpha] - 2026-09-18

### Added

- Added a unified `gribs_config.json` file for all standard user inputs.
- Added selection between `cea2` and `legacy_fit` thermochemistry backends.
- Added propellant reactant composition, constituent density, and initial-temperature input in the unified JSON configuration.
- Added ideal additive-volume mixture-density calculation:

  ```text
  rho_ideal = 1 / sum_i(w_i / rho_i)
  ```

- Added packing-fraction correction for bulk propellant density.
- Added JSON configuration for the burn law, erosive burning, grain geometry, nozzle geometry, flow separation, throat erosion, environment, igniter, thermochemistry, solver, blowdown, and outputs.
- Added configuration source, original configuration, resolved runtime configuration, and derived inputs to `summary.json`.
- Added script-relative discovery of the default JSON file and local external CEA runtime.
- Added explicit reporting of ideal mixture density, packing fraction, and bulk propellant density.

### Changed

- Moved standard input management from source-level Python edits to a complete JSON configuration.
- Moved CEA reactants from a separate reactant JSON file into `gribs_config.json`.
- Made CEA2 a selectable JSON backend rather than a source-level selection.
- Changed the runtime `Config` dataclass into an internal resolved-configuration container rather than the standard user input surface.
- Made the output directory and output filenames configurable through JSON.
- Updated documentation and source comments for the v0.3 configuration architecture.

### Fixed

- Added the missing return value from CEA batch parsing.
- Added script-relative lookup for `fcea2` to support Visual Studio Code Run and arbitrary shell working directories.
- Added explicit executable-permission diagnostics on POSIX systems.
- Reduced each legacy CEA request to eight pressure points for compatibility with the tested `fcea2` output limits while retaining 81 total grid points.
- Removed obsolete source-level command-line override examples.
- Removed obsolete `--list_params` documentation and residual CLI handling.
- Corrected duplicate source section numbering.
- Clarified that CEA product-gas density is not solid grain bulk density.

### Preserved

The following established calculation logic was retained:

- unsteady chamber-pressure equation
- cylindrical-bore grain geometry
- Saint-Robert burn law
- optional Lenoir-Robert erosive burning
- choked and unchoked nozzle branches
- optional flow separation
- optional throat erosion
- post-burnout blowdown
- pressure-dependent gas-property interpolation
- LSODA, BDF, and Radau integration options
- mass and impulse quadratures as ODE states
- numerical self-tests and verification metrics

### External dependency notice

The legacy NASA CEA executable and databases are not distributed by GRIBS. Users selecting the `cea2` backend must separately provide compatible copies of `fcea2`, `thermo.lib`, and `trans.lib` and must comply with the applicable terms.

### Known limitations

- Alpha release, not certified for safety-critical use.
- External CEA compatibility must be validated by each user.
- The constituent-density mixing rule assumes additive constituent volumes before packing correction.
- Burn-rate and loss coefficients remain formulation-specific and hardware-specific empirical inputs.
