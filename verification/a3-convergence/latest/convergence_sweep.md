# Phase A-3 Numerical Convergence Sweep

## Reproducibility

- UTC execution time: `2026-10-05T05:33:05.521543+00:00`
- GRIBS version: `0.5.0-alpha`
- Schema version: `0.5.0-alpha`
- Python: `3.14.4`
- CEA package: `3.3.4`
- Git commit: `955d0d362792ffdc5f5fcbda3cf0583c03c181cb`
- Git working tree clean before outputs: `False`
- Configuration: `/home/hikaru-kurokawa/Rocket/GRIBS-dev/gribs_config.json`
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
| rtol | 1e-07 | convergence | 1e-11 | PASS | 3.16809e+06 | 5.8e-08 | 11.236 | 1.59e-08 | 26.9011 | 5.96e-08 | 131.511 | 3.77e-08 | 4.21e-06 | 56.9 | 3.37e-08 | 339 | 2.51 |
| rtol | 1e-08 | convergence | 1e-11 | PASS | 3.16809e+06 | 1.64e-10 | 11.236 | 2.8e-09 | 26.9011 | 1.68e-10 | 131.511 | 8.51e-09 | 9.56e-06 | 130 | 8.21e-09 | 81.9 | 2.67 |
| rtol | 1e-09 | convergence | 1e-11 | PASS | 3.16809e+06 | 7.76e-10 | 11.236 | 5.66e-11 | 26.9011 | 7.97e-10 | 131.511 | 2.36e-10 | 4.54e-07 | 5.25 | 3.05e-10 | 2.08 | 2.68 |
| rtol | 1e-10 | convergence | 1e-11 | PASS | 3.16809e+06 | 3.55e-10 | 11.236 | 2.53e-12 | 26.9011 | 3.65e-10 | 131.511 | 4.55e-11 | 9.43e-08 | 0.297 | 1.47e-10 | 0.487 | 2.66 |
| rtol | 1e-11 | convergence | 1e-11 | PASS | 3.16809e+06 | 0 | 11.236 | 0 | 26.9011 | 0 | 131.511 | 0 | 7.27e-08 | 0 | 9.9e-11 | 0 | 2.65 |
| atol_p | 1e-01 | convergence | 1e-04 | PASS | 3.16809e+06 | 5.99e-08 | 11.236 | 2.27e-10 | 26.9011 | 6.16e-08 | 131.511 | 3.01e-09 | 5.88e-05 | 2.78e+03 | 5.65e-11 | 0.966 | 2.49 |
| atol_p | 1e-02 | convergence | 1e-04 | PASS | 3.16809e+06 | 7.29e-10 | 11.236 | 7.6e-10 | 26.9011 | 7.49e-10 | 131.511 | 1.77e-09 | 1.89e-06 | 88.5 | 1.95e-10 | 0.882 | 2.57 |
| atol_p | 1e-03 | convergence | 1e-04 | PASS | 3.16809e+06 | 8.91e-11 | 11.236 | 5.04e-10 | 26.9011 | 9.15e-11 | 131.511 | 1.51e-09 | 4.54e-07 | 20.5 | 3.05e-10 | 0.815 | 2.67 |
| atol_p | 1e-04 | convergence | 1e-04 | PASS | 3.16809e+06 | 0 | 11.236 | 0 | 26.9011 | 0 | 131.511 | 0 | 2.12e-08 | 0 | 1.65e-09 | 0 | 2.69 |
| max_step_burn | 0.2 | convergence | 0.005 | PASS | 3.16809e+06 | 1.26e-10 | 11.236 | 1.22e-09 | 26.9011 | 1.3e-10 | 131.511 | 6.09e-10 | 1.03e-06 | 1.51 | 8.52e-10 | 0.37 | 2.6 |
| max_step_burn | 0.05 | convergence | 0.005 | PASS | 3.16809e+06 | 6.6e-10 | 11.236 | 2.83e-10 | 26.9011 | 6.78e-10 | 131.511 | 4.3e-10 | 2.99e-07 | 0.271 | 9.45e-10 | 0.3 | 2.61 |
| max_step_burn | 0.02 | convergence | 0.005 | PASS | 3.16809e+06 | 3.7e-10 | 11.236 | 5.84e-10 | 26.9011 | 3.8e-10 | 131.511 | 1.16e-09 | 4.54e-07 | 0.109 | 3.05e-10 | 0.774 | 2.69 |
| max_step_burn | 0.005 | convergence | 0.005 | PASS | 3.16809e+06 | 0 | 11.236 | 0 | 26.9011 | 0 | 131.511 | 0 | 4.1e-07 | 0 | 1.35e-09 | 0 | 3.04 |
| cea_pressure_points | 21 | convergence | 161 | PASS | 3.16809e+06 | 4.65e-08 | 11.2358 | 8.5e-06 | 26.9011 | 5.08e-08 | 131.512 | 2.75e-06 | 2.52e-07 | 0.14 | 2.77e-10 | 0.564 | 2.55 |
| cea_pressure_points | 41 | convergence | 161 | PASS | 3.16809e+06 | 2.08e-09 | 11.2363 | 2.82e-05 | 26.9011 | 2.25e-09 | 131.51 | 8.6e-06 | 9.59e-08 | 0.673 | 7.35e-10 | 0.158 | 2.63 |
| cea_pressure_points | 81 | convergence | 161 | PASS | 3.16809e+06 | 1.15e-11 | 11.236 | 9.45e-06 | 26.9011 | 3.78e-11 | 131.511 | 2.91e-06 | 4.54e-07 | 0.549 | 3.05e-10 | 0.519 | 2.68 |
| cea_pressure_points | 161 | convergence | 161 | PASS | 3.16809e+06 | 0 | 11.2359 | 0 | 26.9011 | 0 | 131.512 | 0 | 2.93e-07 | 0 | 6.34e-10 | 0 | 2.63 |
| unchoked_policy | switch | event_sensitivity | switch | PASS | 3.16809e+06 | 0 | 11.236 | 0 | 26.9011 | 0 | 131.511 | 0 | 4.54e-07 | 0 | 3.05e-10 | 0 | 2.7 |
| unchoked_policy | stop | event_sensitivity | switch | PASS | 3.16809e+06 | 0 | 11.236 | 0 | 26.9011 | 0 | 131.509 | 1.38e-05 | 2.72e-07 | 0.401 | 3.05e-10 | 0 | 2.68 |

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
- `cea_pressure_points=21`: burn_time_s: relative difference 8.504341403091888e-06 exceeds 1e-06; impulse_total_Ns: relative difference 2.754151806793966e-06 exceeds 1e-06
- `cea_pressure_points=41`: burn_time_s: relative difference 2.818141979149377e-05 exceeds 1e-06; impulse_total_Ns: relative difference 8.599685711265884e-06 exceeds 1e-06
- `cea_pressure_points=81`: burn_time_s: relative difference 9.448902077620104e-06 exceeds 1e-06; impulse_total_Ns: relative difference 2.910897423664776e-06 exceeds 1e-06

## Case failures

No case failures.
