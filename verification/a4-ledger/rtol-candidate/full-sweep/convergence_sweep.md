# Phase A-3 Numerical Convergence Sweep

## Reproducibility

- UTC execution time: `2026-10-06T02:27:21.554752+00:00`
- GRIBS version: `0.5.0-alpha`
- Schema version: `0.5.0-alpha`
- Python: `3.14.4`
- CEA package: `3.3.4`
- Git commit: `f12be541bbec4956c97e1033467886d02b079467`
- Git working tree clean before outputs: `False`
- Configuration: `/home/hikaru-kurokawa/Rocket/GRIBS-dev/verification/a4-ledger/rtol-candidate/gribs_config_rtol_1e-10.json`
- Case count: `20`

## Acceptance summary

- Enforcement enabled: `False`
- Evaluated rows: `20`
- Passed rows: `16`
- Failed rows: `4`
- Unevaluated rows: `0`
- Total violations: `8`
- All evaluated rows pass: `False`

Acceptance failures are reported even when enforcement is disabled. With enforcement disabled, tolerance violations do not change the process exit status.

## Results

| setting | value | kind | reference | status | p_max [Pa] | rel. diff | burn time [s] | rel. diff | F_max [N] | rel. diff | impulse [N s] | rel. diff | mass consistency | rel. diff | propellant balance | rel. diff | wall [s] |
|---|---:|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| rtol | 1e-06 | convergence | 1e-11 | PASS | 3.16809e+06 | 1.28e-06 | 11.236 | 3.39e-08 | 26.9011 | 1.31e-06 | 131.511 | 3.46e-08 | 3.24e-05 | 444 | 2.64e-08 | 266 | 2.51 |
| rtol | 1e-07 | convergence | 1e-11 | PASS | 3.16809e+06 | 5.8e-08 | 11.236 | 1.59e-08 | 26.9011 | 5.96e-08 | 131.511 | 3.77e-08 | 4.21e-06 | 56.9 | 3.37e-08 | 339 | 2.45 |
| rtol | 1e-08 | convergence | 1e-11 | PASS | 3.16809e+06 | 1.64e-10 | 11.236 | 2.8e-09 | 26.9011 | 1.68e-10 | 131.511 | 8.51e-09 | 9.56e-06 | 130 | 8.21e-09 | 81.9 | 2.49 |
| rtol | 1e-09 | convergence | 1e-11 | PASS | 3.16809e+06 | 7.76e-10 | 11.236 | 5.66e-11 | 26.9011 | 7.97e-10 | 131.511 | 2.36e-10 | 4.54e-07 | 5.25 | 3.05e-10 | 2.08 | 2.67 |
| rtol | 1e-10 | convergence | 1e-11 | PASS | 3.16809e+06 | 3.55e-10 | 11.236 | 2.53e-12 | 26.9011 | 3.65e-10 | 131.511 | 4.55e-11 | 9.43e-08 | 0.297 | 1.47e-10 | 0.487 | 2.58 |
| rtol | 1e-11 | convergence | 1e-11 | PASS | 3.16809e+06 | 0 | 11.236 | 0 | 26.9011 | 0 | 131.511 | 0 | 7.27e-08 | 0 | 9.9e-11 | 0 | 2.58 |
| atol_p | 1e-01 | convergence | 1e-04 | PASS | 3.16809e+06 | 7.2e-10 | 11.236 | 3.88e-11 | 26.9011 | 7.4e-10 | 131.511 | 2.4e-10 | 2.71e-07 | 5.15 | 9.95e-12 | 0.964 | 2.45 |
| atol_p | 1e-02 | convergence | 1e-04 | PASS | 3.16809e+06 | 4.82e-10 | 11.236 | 1.23e-10 | 26.9011 | 4.95e-10 | 131.511 | 6.58e-10 | 8.82e-06 | 199 | 7.62e-11 | 0.724 | 2.5 |
| atol_p | 1e-03 | convergence | 1e-04 | PASS | 3.16809e+06 | 3.58e-10 | 11.236 | 5.82e-11 | 26.9011 | 3.68e-10 | 131.511 | 1.47e-10 | 9.43e-08 | 1.14 | 1.47e-10 | 0.466 | 2.59 |
| atol_p | 1e-04 | convergence | 1e-04 | PASS | 3.16809e+06 | 0 | 11.236 | 0 | 26.9011 | 0 | 131.511 | 0 | 4.4e-08 | 0 | 2.76e-10 | 0 | 2.73 |
| max_step_burn | 0.2 | convergence | 0.005 | PASS | 3.16809e+06 | 3.41e-10 | 11.236 | 4.22e-10 | 26.9011 | 3.51e-10 | 131.511 | 3.42e-10 | 1.34e-07 | 1.28 | 1.55e-10 | 0.194 | 2.53 |
| max_step_burn | 0.05 | convergence | 0.005 | PASS | 3.16809e+06 | 1.4e-10 | 11.236 | 1.83e-10 | 26.9011 | 1.43e-10 | 131.511 | 3.2e-10 | 2.8e-08 | 0.523 | 1.53e-10 | 0.177 | 2.53 |
| max_step_burn | 0.02 | convergence | 0.005 | PASS | 3.16809e+06 | 1.92e-10 | 11.236 | 4.76e-11 | 26.9011 | 1.97e-10 | 131.511 | 3.09e-12 | 9.43e-08 | 0.607 | 1.47e-10 | 0.136 | 2.6 |
| max_step_burn | 0.005 | convergence | 0.005 | PASS | 3.16809e+06 | 0 | 11.236 | 0 | 26.9011 | 0 | 131.511 | 0 | 5.87e-08 | 0 | 1.3e-10 | 0 | 3.04 |
| cea_pressure_points | 21 | convergence | 161 | PASS | 3.16809e+06 | 4.7e-08 | 11.2358 | 8.5e-06 | 26.9011 | 5.14e-08 | 131.512 | 2.75e-06 | 1.46e-07 | 0.339 | 7.74e-11 | 0.628 | 2.54 |
| cea_pressure_points | 41 | convergence | 161 | PASS | 3.16809e+06 | 2.05e-09 | 11.2363 | 2.82e-05 | 26.9011 | 2.21e-09 | 131.51 | 8.6e-06 | 1.77e-07 | 0.622 | 2.64e-10 | 0.269 | 2.64 |
| cea_pressure_points | 81 | convergence | 161 | PASS | 3.16809e+06 | 6.48e-10 | 11.236 | 9.45e-06 | 26.9011 | 6.92e-10 | 131.511 | 2.91e-06 | 9.43e-08 | 0.138 | 1.47e-10 | 0.293 | 2.59 |
| cea_pressure_points | 161 | convergence | 161 | PASS | 3.16809e+06 | 0 | 11.2359 | 0 | 26.9011 | 0 | 131.512 | 0 | 1.09e-07 | 0 | 2.08e-10 | 0 | 2.61 |
| unchoked_policy | switch | event_sensitivity | switch | PASS | 3.16809e+06 | 0 | 11.236 | 0 | 26.9011 | 0 | 131.511 | 0 | 9.43e-08 | 0 | 1.47e-10 | 0 | 2.62 |
| unchoked_policy | stop | event_sensitivity | switch | PASS | 3.16809e+06 | 0 | 11.236 | 0 | 26.9011 | 0 | 131.509 | 1.38e-05 | 5.89e-08 | 0.375 | 1.47e-10 | 0 | 2.59 |

## Diagnostic comparison

| setting | value | burn stop | matches reference | blowdown stop | matches reference | extrapolation | matches reference | regimes | matches reference |
|---|---:|---|---|---|---|---|---|---|---|
| rtol | 1e-06 | burnout | True | ambient | True | False | True | ["choked", "subsonic"] | True |
| rtol | 1e-07 | burnout | True | ambient | True | False | True | ["choked", "subsonic"] | True |
| rtol | 1e-08 | burnout | True | ambient | True | False | True | ["choked", "subsonic"] | True |
| rtol | 1e-09 | burnout | True | ambient | True | False | True | ["choked", "subsonic"] | True |
| rtol | 1e-10 | burnout | True | ambient | True | False | True | ["choked", "subsonic"] | True |
| rtol | 1e-11 | burnout | True | ambient | True | False | True | ["choked", "subsonic"] | True |
| atol_p | 1e-01 | burnout | True | ambient | True | False | True | ["choked", "subsonic"] | True |
| atol_p | 1e-02 | burnout | True | ambient | True | False | True | ["choked", "subsonic"] | True |
| atol_p | 1e-03 | burnout | True | ambient | True | False | True | ["choked", "subsonic"] | True |
| atol_p | 1e-04 | burnout | True | ambient | True | False | True | ["choked", "subsonic"] | True |
| max_step_burn | 0.2 | burnout | True | ambient | True | False | True | ["choked", "subsonic"] | True |
| max_step_burn | 0.05 | burnout | True | ambient | True | False | True | ["choked", "subsonic"] | True |
| max_step_burn | 0.02 | burnout | True | ambient | True | False | True | ["choked", "subsonic"] | True |
| max_step_burn | 0.005 | burnout | True | ambient | True | False | True | ["choked", "subsonic"] | True |
| cea_pressure_points | 21 | burnout | True | ambient | True | False | True | ["choked", "subsonic"] | True |
| cea_pressure_points | 41 | burnout | True | ambient | True | False | True | ["choked", "subsonic"] | True |
| cea_pressure_points | 81 | burnout | True | ambient | True | False | True | ["choked", "subsonic"] | True |
| cea_pressure_points | 161 | burnout | True | ambient | True | False | True | ["choked", "subsonic"] | True |
| unchoked_policy | switch | burnout | True | ambient | True | False | True | ["choked", "subsonic"] | True |
| unchoked_policy | stop | burnout | True | choking-loss | False | False | True | ["choked", "subsonic"] | True |

## Acceptance violations

- `rtol=1e-06`: p_max_Pa: relative difference 1.277129730262344e-06 exceeds 1e-06; F_max_N: relative difference 1.311957647163808e-06 exceeds 1e-06
- `cea_pressure_points=21`: burn_time_s: relative difference 8.504344005968411e-06 exceeds 1e-06; impulse_total_Ns: relative difference 2.75437198895187e-06 exceeds 1e-06
- `cea_pressure_points=41`: burn_time_s: relative difference 2.818138995539964e-05 exceeds 1e-06; impulse_total_Ns: relative difference 8.599723375317607e-06 exceeds 1e-06
- `cea_pressure_points=81`: burn_time_s: relative difference 9.448954226523935e-06 exceeds 1e-06; impulse_total_Ns: relative difference 2.910630681742755e-06 exceeds 1e-06

## Case failures

No case failures.
