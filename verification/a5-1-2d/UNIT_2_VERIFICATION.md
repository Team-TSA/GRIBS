# A-5.1-2d unit 2 verification

## Scope

Unit 2 changes only the external inspection and failure-reporting
contract:

- failed post-run consistency is exposed as
  `W_POST_RUN_CONSISTENCY`;
- failed mass-ledger acceptance is exposed as `W_MASS_LEDGER`;
- either failed verdict produces process exit status 1;
- diagnostic CSV, summary text, and summary JSON are written before the
  verdict-based exit;
- machine-readable diagnostics are written before figure generation;
- igniter commanded and injected mass remain distinct from the
  propellant residual.

Unit 2 does not change the simulation equations, ledger equations,
ledger tolerances, thermochemistry, cache semantics, nozzle flow, exit
velocity, or thrust.

## rtol = 1e-10

- process exit: 0
- output files: 4
- post-run verdict: PASS
- mass-ledger verdict: PASS
- warning codes: none
- normalized production history: byte-identical to the pre-change
  baseline
- protected trajectory and performance results: unchanged

Normalized history SHA-256:

    9f738447e76d8453e68ac3f61a1de94ec46c52da4288f70c515fece73993fdbb

## rtol = 1e-7

- process exit: 1
- output files preserved: 4
- post-run verdict: PASS
- mass-ledger verdict: FAIL
- warning code: `W_MASS_LEDGER`
- total-system verdict: PASS
- propellant verdict: FAIL

Measured ledger values:

    total-system relative residual = 2.9416794564934683e-09
    total-system tolerance = 1e-06
    propellant relative residual = 4.621689212596392e-08
    propellant tolerance = 1e-09
    igniter commanded = 0.002 kg
    igniter injected = 0.001999999975910654 kg

The CSV, text summary, and JSON summary were written before the figure.
All four output artifacts existed and were non-empty before the process
returned status 1.

## Human-readable diagnostics

The failing text summary contains:

- `[ verification warnings ]`;
- `W_MASS_LEDGER`;
- `mass-ledger result : FAIL`;
- an explicit note that the legacy propellant metric includes igniter
  mass in its numerator;
- the invariant EOS-history metric;
- the deprecated shrinking-mass metric.
