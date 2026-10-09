# A-5.1-2d completion record

## Scope completed

A-5.1-2d changes inspection and failure-reporting contracts without
changing the ordinary accepted numerical trajectory.

Implemented:

- structured post-run mass-consistency reporting;
- one time-invariant EOS-history normalization;
- retention of the former shrinking-denominator metric as deprecated;
- failed post-run and mass-ledger warning contracts;
- nonzero CLI exit for either failed verdict;
- diagnostic output preservation before verdict-based exit;
- explicit igniter context in mass-ledger warnings;
- exact self-test name and count contract;
- reference_constant single-phase fraction protection.

Not changed:

- chamber Psi or dPsi/dp;
- chamber ODE;
- nozzle total flow;
- exit velocity;
- thrust;
- thermochemistry table construction;
- cache schema or meaning;
- mass-ledger equations;
- mass-ledger tolerances;
- empirical coefficients.

## Numerical invariance

For the accepted rtol=1e-10 wide case, the normalized history CSV is
byte-identical before and after A-5.1-2d.

Both files have SHA-256:

    9f738447e76d8453e68ac3f61a1de94ec46c52da4288f70c515fece73993fdbb

Existing summary fields are recursively identical except for:

- configuration and output path differences caused by isolated output
  directories;
- the four explicitly approved A-5.1-2d diagnostic fields.

The recursive comparison recorded:

    difference_count = 10
    unexpected_difference_count = 0

## EOS-history metric

The new gate uses a run-wide, time-invariant balance-mass reference.

For rtol=1e-10:

    legacy metric = 3.601669931212943e-07
    invariant metric = 3.8407655452640866e-09
    verdict = PASS

For rtol=1e-7:

    legacy metric = 3.0474506827047328e-05
    invariant metric = 3.187461057208937e-07
    verdict = PASS

The legacy shrinking-mass metric remains available as a deprecated
diagnostic and is not used as the gate.

## Failure-output behavior

For rtol=1e-7:

- post-run consistency passes;
- the mass ledger fails;
- the total-system ledger passes;
- the propellant ledger fails;
- `W_MASS_LEDGER` is emitted;
- CSV, text summary, JSON summary, and figure are preserved;
- all four output artifacts are non-empty;
- the process returns status 1.

Measured values:

    total-system relative residual = 2.9416794564934683e-09
    total-system tolerance = 1e-06
    propellant relative residual = 4.621689212596392e-08
    propellant tolerance = 1e-09
    igniter commanded = 0.002 kg
    igniter injected = 0.001999999975910654 kg

Immediately after the runtime verification, file timestamps showed the
CSV, text summary, and JSON summary written before the figure.

The later `failure-output-contract-after.json` timestamp values are not
used as the primary write-order proof because the CSV was subsequently
normalized from CRLF to LF for repository storage, which updated its
mtime. That JSON remains evidence of existence, nonzero size, verdict,
warning code, and exit status.

## Self-test contract

The expected self-test contract consists of:

- 14 common named results;
- one nozzle-geometry-dependent named result;
- 15 unique results in exact order.

Default converging-nozzle configuration:

    process exit = 0
    reported results = 15
    contract errors = 0

C-D nozzle configuration:

    process exit = 0
    reported results = 15
    contract errors = 0

A missing, unexpected, duplicated, or reordered result causes the
self-test contract to fail.

## reference_constant contract

The existing product implementation is protected without adding new
configuration freedom:

    phase_fractions() = (1.0, 0.0)

Metadata:

    gas_mass_fraction = 1.0
    condensed_mass_fraction = 0.0

Thermochemistry dump:

    gas_mass_fraction = [1.0, 1.0]
    condensed_mass_fraction = [0.0, 0.0]

The existing reference_constant test suite passed:

    22 passed

## Final tests

The final full test suite passed:

    83 passed

The normalized production path and protected summary values remain
unchanged for the accepted rtol=1e-10 comparison case.

## Change boundary

Product and test changes are limited to:

- `GRIBS_v0.5.0-alpha.py`;
- `tests/test_gribs.py`.

The comparison policy and A-5.1-2d verification artifacts are contained
under `verification/a5-1-2d/`.

A-5.1-4 physical mass-basis corrections are not included in this change.
