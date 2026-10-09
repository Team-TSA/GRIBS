# A-5.1-2d unit 3 verification

## Self-test name contract

The existing physical self-tests remain unchanged.

The expected contract is:

- 14 common named results;
- one nozzle-geometry-dependent result;
- 15 unique results in exact order.

For a converging nozzle, the geometry-dependent result is:

    converging-nozzle continuity (mdot, F)

For a C-D nozzle, it is:

    C-D nozzle transition mass-flow continuity

A missing, unexpected, duplicated, or reordered result causes the
self-test contract to fail. It cannot silently disappear while the
self-test command reports success.

## Runtime verification

Default configuration:

    process exit = 0
    reported results = 15
    contract errors = 0

C-D nozzle configuration:

    process exit = 0
    reported results = 15
    contract errors = 0

## reference_constant phase-fraction contract

The existing product implementation is protected without introducing
new configuration freedom:

    phase_fractions() = (1.0, 0.0)

Metadata:

    gas_mass_fraction = 1.0
    condensed_mass_fraction = 0.0

Thermochemistry dump:

    gas_mass_fraction = [1.0, 1.0]
    condensed_mass_fraction = [0.0, 0.0]

The existing reference_constant test suite passed:

    22 passed
