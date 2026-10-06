# A-4 Igniter Ledger Convergence Investigation

| case | rtol | max step [s] | propellant signed residual [kg] | propellant relative residual | limit | pass | total relative residual |
|---|---:|---:|---:|---:|---:|:---:|---:|
| baseline | 1.000e-09 | 2.000e-02 | 1.1180755626893912e-10 | 1.8624912520418528e-09 | 1.000e-09 | FAIL | 2.4309253994933600e-11 |
| rtol_1e-10 | 1.000e-10 | 2.000e-02 | 2.0919981025269152e-11 | 3.4848522723029594e-10 | 1.000e-09 | PASS | 8.5825276658719562e-13 |
| max_step_0p005 | 1.000e-09 | 5.000e-03 | 7.2357286828861334e-11 | 1.2053283180259971e-09 | 1.000e-09 | FAIL | 4.4239271460461145e-11 |
| tight_combined | 1.000e-10 | 5.000e-03 | 1.2993946196804274e-11 | 2.1645326960583694e-10 | 1.000e-09 | PASS | 3.2701682686059207e-12 |

## Automated interpretation

The propellant residual decreases under tighter integration settings by a factor of 8.60459. This indicates a numerical-integration contribution.
