# A-5.1-2d unit 1 verification

## Scope

Unit 1 adds one structured post-run consistency report and replaces the
EOS-history gate denominator with a time-invariant run-wide balance-mass
reference.

It does not change chamber Psi, nozzle flow, exit velocity, thrust,
thermochemistry tables, cache semantics, or ledger tolerances.

## Tests

- Python compilation: PASS
- targeted pytest: 3 passed
- full pytest: 77 passed

## rtol = 1e-10

- process exit: 0
- output files: 4
- legacy EOS-history metric:
  `3.601669931212943e-07`
- invariant EOS-history metric:
  `3.8407655452640866e-09`
- invariant reference mass:
  `0.0061107926089148945 kg`
- post-run verdict: PASS
- ledger verdict: PASS

The normalized production-history CSV is byte-identical to the
pre-change baseline:

    9f738447e76d8453e68ac3f61a1de94ec46c52da4288f70c515fece73993fdbb

Protected performance and trajectory summary values are identical.

## rtol = 1e-7

- process exit after unit 1: 0
- output files: 4
- legacy EOS-history metric:
  `3.0474506827047328e-05`
- invariant EOS-history metric:
  `3.187461057208937e-07`
- invariant reference mass:
  `0.006110793878620346 kg`
- post-run verdict: PASS
- ledger verdict: FAIL

Ledger details:

- total-system verdict: PASS
- total-system relative residual:
  `2.9416794564934683e-09`
- total-system tolerance: `1e-06`
- propellant verdict: FAIL
- propellant relative residual:
  `4.621689212596392e-08`
- propellant tolerance: `1e-09`
- ledger warning codes emitted: none
- ledger failure reflected in process exit: no

This dynamically reproduces the S-2 verdict-visibility and exit-status
defect.

## Summary comparison

The recursive comparison permits only:

- configuration-source path differences;
- output-file path differences;
- original-configuration output-directory differences;
- the four explicitly approved A-5.1-2d result fields.

All other pre-existing summary fields must remain identical.
See `unit-1-summary-comparison.txt`.
