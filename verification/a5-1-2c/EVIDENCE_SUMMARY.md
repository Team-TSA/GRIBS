# A-5.1-2c pre-fix intermediate-condensation evidence

- Baseline commit: 0c80c37ae061915e9b20cba3a6c230122f84aba5
- Product-code change in this evidence unit: none
- Thermochemistry backend: cea_python
- CEA package/library: 3.3.4
- Solver comparison setting: rtol = 1e-10
- Thermochemistry temperature efficiency: 1.0
- Property-table range: 0.1 MPa to 20 MPa
- Property-table points: 81

## Condensation range

    Yg_min = 0.7238397417377189
    Yg_max = 0.7265536177401214
    Yc_min = 0.2734463822598786
    Yc_max = 0.2761602582622811
    max |Yg + Yc - 1| = 0
    Yg_min < 0.9 = true

The case satisfies the A-5.1-2c requirement for a mechanically
distinguishable condensed pre-fix case.

## Included evidence

- derived wide configuration;
- pre-fix thermochemistry dump;
- condensation summary;
- production time history;
- text and JSON summaries;
- self-test log and exit status;
- pytest log and exit status;
- dump log, timing, and exit status;
- runtime and Git environment;
- evidence references;
- SHA-256 manifest.

## Storage normalization

The generated production history originally used CRLF line endings.
Before committing the evidence, it was normalized to LF so that the
repository `git diff --check` contract remains valid. CSV field values,
row order, and column order were not changed.

## Scope

This is pre-fix evidence. It records the behavior and thermochemistry
table before the accepted mass-basis correction. It is not corrected
golden data and is not an engineering-validation result.
