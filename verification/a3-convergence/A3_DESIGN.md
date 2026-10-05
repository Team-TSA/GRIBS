# Phase A-3 Numerical Convergence Harness Design

## Status

Implemented and formally verified.

The formal 20-case evidence is stored under:

- `verification/a3-convergence/latest/convergence_sweep.csv`
- `verification/a3-convergence/latest/convergence_sweep.md`
- `verification/a3-convergence/latest/convergence_sweep_metadata.json`

## Source material

The historical prototype is preserved unchanged in:

- `Specification/convergence_sweep.py`
- `Specification/convergence_sweep.json`
- `Specification/convergence_sweep.txt`

The prototype was created against an earlier GRIBS implementation. It is
treated as investigation evidence and design input, not as executable
production code for GRIBS v0.5.0-alpha.

## Scope

Phase A-3 converts the historical convergence investigation into a
repeatable harness for the current `GRIBS_v0.5.0-alpha.py`.

The harness shall:

1. run 20 cases with one command;
2. sweep solver, thermochemistry-table, and event-policy settings;
3. use the current GRIBS configuration loader and calculation functions;
4. reinitialize thermochemistry for every case;
5. write raw and relative-difference results to CSV;
6. write the same comparison as a Markdown report;
7. record reproducibility metadata;
8. return a nonzero exit status if any case fails.

Phase A-3 shall not:

- change the physical model;
- expose new event thresholds in the main configuration;
- remove the module-global thermochemistry backend;
- perform the Phase A-6 package reorganization;
- implement the reserved nozzle `shock` transition policy.

## Case matrix

### Solver convergence

- `rtol`: `1e-6`, `1e-7`, `1e-8`, `1e-9`, `1e-10`, `1e-11`
- `atol_p`: `1e-1`, `1e-2`, `1e-3`, `1e-4` Pa
- `max_step_burn`: `0.2`, `0.05`, `0.02`, `0.005` s

### Thermochemistry-table convergence

- `cea_pressure_points`: `21`, `41`, `81`, `161`

### Event-policy sensitivity

- `unchoked_policy`: `switch`, `stop`

Total: 20 cases.

## Reference cases

The convergence groups use their finest case as the reference:

- `rtol`: `1e-11`
- `atol_p`: `1e-4` Pa
- `max_step_burn`: `0.005` s
- `cea_pressure_points`: `161`

The event-policy group is categorical rather than ordered. Its reference is
the normal current policy:

- `unchoked_policy`: `switch`

Rows shall identify whether their comparison is `convergence` or
`event_sensitivity`.

## Compared quantities

Numeric quantities:

- `p_max_Pa`
- `burn_time_s`
- `F_max_N`
- `impulse_total_Ns`
- `total_mass_consistency_error`
- `propellant_mass_balance_error`

Categorical and diagnostic quantities:

- `burn_stop_reason`
- `blowdown_stop_reason`
- `property_extrapolation`
- `regimes_visited`

Relative differences are calculated only for finite numeric values. Diagnostic
quantities are compared for equality.

## Outputs

Default output directory:

`verification/a3-convergence/latest/`

Required files:

- `convergence_sweep.csv`
- `convergence_sweep.md`
- `convergence_sweep_metadata.json`

The metadata file shall include, where available:

- UTC execution time;
- GRIBS program and schema versions;
- Python version;
- CEA package version;
- Git commit;
- Git working-tree state;
- command line;
- source configuration path;
- case count.

## Failure behavior

A failed case shall be recorded with its exception type and message. Remaining
cases may continue so that all failures are visible, but the command shall
return a nonzero exit status if any case fails.

## Known architectural limitation

GRIBS v0.5.0-alpha still uses the module-global `_THERMO` backend. The A-3
harness therefore runs cases sequentially and explicitly initializes
thermochemistry before every case. Removal of this global state belongs to
Phase A-6.

## Implemented files

- `tools/convergence_sweep.py`
- `tools/convergence_sweep_config.json`
- `tests/test_convergence_sweep.py`
- `verification/a3-convergence/latest/convergence_sweep.csv`
- `verification/a3-convergence/latest/convergence_sweep.md`
- `verification/a3-convergence/latest/convergence_sweep_metadata.json`

## Acceptance interpretation

The default sweep configuration records tolerance results but does not enforce
them as a process exit condition:

- `acceptance.enforce = false`

This is intentional. The historical Phase A investigation reported differences
of approximately `1e-6` or less, but the current GRIBS v0.5.0-alpha model
produces several reproducible differences above `1e-6`.

A tolerance violation is distinct from a case execution failure:

- case execution failures always cause a nonzero exit status;
- tolerance violations are written to CSV, Markdown, and metadata;
- tolerance violations cause a nonzero exit status only when
  `acceptance.enforce` is `true`.

Residual quantities are evaluated by absolute limits rather than by relative
change against a very small reference residual.

## Formal verification result

The formal sweep ran all 20 configured cases successfully with one command.

Case distribution:

- `rtol`: 6
- `atol_p`: 4
- `max_step_burn`: 4
- `cea_pressure_points`: 4
- `unchoked_policy`: 2

All 20 GRIBS calculations passed. The non-enforced acceptance evaluation
reported:

- evaluated rows: 20
- passed rows: 16
- rows with tolerance violations: 4
- total tolerance violations: 8
- case execution failures: 0

The four rows with tolerance violations were:

1. `rtol=1e-6`
   - `p_max_Pa`: `1.277129730262344e-6`
   - `F_max_N`: `1.311957647163808e-6`
2. `cea_pressure_points=21`
   - `burn_time_s`: `8.504341403091888e-6`
   - `impulse_total_Ns`: `2.754151806793966e-6`
3. `cea_pressure_points=41`
   - `burn_time_s`: `2.818141979149377e-5`
   - `impulse_total_Ns`: `8.599685711265884e-6`
4. `cea_pressure_points=81`
   - `burn_time_s`: `9.448902077620104e-6`
   - `impulse_total_Ns`: `2.910897423664776e-6`

The maximum observed relative differences among the primary convergence
quantities were:

- `p_max_Pa`: `1.2771297302623438e-6` at `rtol=1e-6`
- `burn_time_s`: `2.8181419791493772e-5` at
  `cea_pressure_points=41`
- `F_max_N`: `1.3119576471638082e-6` at `rtol=1e-6`
- `impulse_total_Ns`: `8.5996857112658840e-6` at
  `cea_pressure_points=41`

The maximum observed absolute residuals remained within their configured
limits:

- `total_mass_consistency_error`:
  `5.8751037243915560e-5` at `atol_p=1e-1`
  (configured maximum `1e-4`)
- `propellant_mass_balance_error`:
  `3.3669938392932863e-8` at `rtol=1e-7`
  (configured maximum `1e-7`)

These results are retained as measured numerical-convergence evidence. No
physical model or solver result was changed to force agreement with the
historical prototype.

## Exit-status verification

The enforcement path was tested with a temporary copy of the sweep
configuration using:

- `acceptance.enforce = true`
- cases `rtol=1e-6` and `rtol=1e-11`

Both GRIBS calculations passed. The coarse case exceeded two configured
relative tolerances, and the harness returned exit status `1` as designed.

The temporary enforced configuration and its output were removed after
verification.

