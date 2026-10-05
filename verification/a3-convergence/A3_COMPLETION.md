# Phase A-3 Completion Record

## Status

Phase A-3 is complete.

Result: **PASS**

## Scope completed

Phase A-3 converted the historical numerical-convergence investigation into a
repeatable harness for the current GRIBS v0.5.0-alpha implementation.

The completed implementation provides:

1. a 20-case sweep executable with one command;
2. solver-tolerance sweeps;
3. thermochemistry pressure-table resolution sweeps;
4. event-policy sensitivity cases;
5. externally configured sweep values;
6. externally configured comparison metrics and tolerances;
7. CSV output containing raw values and relative differences;
8. a Markdown comparison report;
9. JSON reproducibility metadata;
10. selectable cases for smoke and diagnostic runs;
11. case-failure recording with nonzero exit status;
12. optional tolerance enforcement with nonzero exit status;
13. unit tests for case construction, configuration validation, comparison,
    failure recording, and acceptance evaluation.

## Implementation files

- `tools/convergence_sweep.py`
- `tools/convergence_sweep_config.json`
- `tests/test_convergence_sweep.py`

## Formal evidence

- `verification/a3-convergence/A3_DESIGN.md`
- `verification/a3-convergence/latest/convergence_sweep.csv`
- `verification/a3-convergence/latest/convergence_sweep.md`
- `verification/a3-convergence/latest/convergence_sweep_metadata.json`

## Formal command

The complete default sweep is reproduced with:

```bash
python tools/convergence_sweep.py \
  --sweep-config tools/convergence_sweep_config.json \
  --config gribs_config.json \
  --output-dir verification/a3-convergence/latest
```

## Case matrix

The formal run contained 20 cases:

- `rtol`: 6
- `atol_p`: 4
- `max_step_burn`: 4
- `cea_pressure_points`: 4
- `unchoked_policy`: 2

Reference cases:

- `rtol=1e-11`
- `atol_p=1e-4`
- `max_step_burn=0.005`
- `cea_pressure_points=161`
- `unchoked_policy=switch`

## Formal result

All 20 GRIBS calculations passed.

The output structure, case counts, reference selection, relative-difference
columns, diagnostic comparisons, metadata, and report sections were
independently checked after execution.

The default acceptance mode is non-enforcing. It reported:

- evaluated rows: 20
- passed rows: 16
- rows with tolerance violations: 4
- violations: 8
- unevaluated rows: 0
- case execution failures: 0

The tolerance violations are retained as measured evidence. They are not hidden
and were not removed by changing the physical model, numerical results, or
historical acceptance threshold.

## Measured convergence observations

The maximum relative differences among the four primary comparison quantities
were:

- `p_max_Pa`: `1.2771297302623438e-6`
- `burn_time_s`: `2.8181419791493772e-5`
- `F_max_N`: `1.3119576471638082e-6`
- `impulse_total_Ns`: `8.5996857112658840e-6`

The largest values occurred in the coarsest solver-tolerance case or in the
thermochemistry pressure-table resolution sweep.

The maximum absolute residuals were:

- `total_mass_consistency_error`:
  `5.8751037243915560e-5`
- `propellant_mass_balance_error`:
  `3.3669938392932863e-8`

Both remained within the configured limits of `1e-4` and `1e-7`,
respectively.

## Thermochemistry reinitialization

The historical prototype did not reliably rebuild the thermochemistry backend
for every pressure-table resolution case.

The A-3 harness explicitly initializes thermochemistry for every case.
Separate 21-point and 161-point cases produced distinct results, demonstrating
that the configured table resolution is applied.

## Event-policy coverage

The 20-case matrix includes:

- `unchoked_policy=switch`
- `unchoked_policy=stop`

Both calculations passed.

The observed blowdown stop reasons were:

- `switch`: `ambient`
- `stop`: `choking-loss`

This confirms that the event-policy cases exercise distinct calculation paths.

## Enforcement verification

A temporary configuration with `acceptance.enforce=true` was tested using the
`rtol=1e-6` and `rtol=1e-11` cases.

Both calculations passed, but the coarse case exceeded two relative
tolerances. The harness returned exit status `1` as designed.

The temporary configuration and output were removed after verification.

## Regression verification

The repository-wide validation completed successfully:

- Python compilation: PASS
- `git diff --check`: PASS
- pytest: 39 passed
- formal 20-case sweep: 20 passed
- formal output structure validation: PASS
- enforced exit-status validation: PASS

## Architectural boundary

Phase A-3 does not:

- change the GRIBS physical model;
- introduce new main-configuration event thresholds;
- remove the module-global thermochemistry backend;
- perform the Phase A-6 package reorganization;
- implement the reserved nozzle `shock` transition policy.

The harness remains sequential because GRIBS v0.5.0-alpha uses module-global
thermochemistry state. Removing that state remains Phase A-6 work.

## Completion decision

The Phase A-3 definition of done is satisfied:

1. the 20-case sweep is reproducible with one command;
2. the result includes reference-relative comparison tables;
3. CSV and Markdown outputs are generated;
4. sweep values, metrics, and tolerances are externally configured;
5. reproducibility metadata is stored;
6. optional tolerance enforcement is implemented and verified;
7. the existing GRIBS regression suite remains green.

Phase A-3 is accepted as complete.
