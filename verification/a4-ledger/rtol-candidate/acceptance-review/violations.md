# A-4 rtol Candidate Acceptance Review

- Sweep rows: 20
- Failed rows: 4
- Violations: 8

## Failed rows

- `cea_pressure_points=21` (convergence)
- `cea_pressure_points=41` (convergence)
- `cea_pressure_points=81` (convergence)
- `rtol=1e-06` (convergence)

## Violations

### rtol=1e-06: p_max_Pa

- Category: `relative`
- Measured: `1.2771297302623438e-06`
- Limit: `9.9999999999999995e-07`
- Ratio to limit: `1.27712973`
- Reference value: `1e-11`

### rtol=1e-06: F_max_N

- Category: `relative`
- Measured: `1.3119576471638082e-06`
- Limit: `9.9999999999999995e-07`
- Ratio to limit: `1.31195765`
- Reference value: `1e-11`

### cea_pressure_points=21: burn_time_s

- Category: `relative`
- Measured: `8.5043440059684106e-06`
- Limit: `9.9999999999999995e-07`
- Ratio to limit: `8.50434401`
- Reference value: `161`

### cea_pressure_points=21: impulse_total_Ns

- Category: `relative`
- Measured: `2.7543719889518697e-06`
- Limit: `9.9999999999999995e-07`
- Ratio to limit: `2.75437199`
- Reference value: `161`

### cea_pressure_points=41: burn_time_s

- Category: `relative`
- Measured: `2.8181389955399643e-05`
- Limit: `9.9999999999999995e-07`
- Ratio to limit: `28.18139`
- Reference value: `161`

### cea_pressure_points=41: impulse_total_Ns

- Category: `relative`
- Measured: `8.5997233753176066e-06`
- Limit: `9.9999999999999995e-07`
- Ratio to limit: `8.59972338`
- Reference value: `161`

### cea_pressure_points=81: burn_time_s

- Category: `relative`
- Measured: `9.4489542265239350e-06`
- Limit: `9.9999999999999995e-07`
- Ratio to limit: `9.44895423`
- Reference value: `161`

### cea_pressure_points=81: impulse_total_Ns

- Category: `relative`
- Measured: `2.9106306817427547e-06`
- Limit: `9.9999999999999995e-07`
- Ratio to limit: `2.91063068`
- Reference value: `161`
