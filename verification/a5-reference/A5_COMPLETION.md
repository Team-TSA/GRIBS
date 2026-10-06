# Phase A-5 Reference Problems Completion

## 1. Completion status

Phase A-5 is complete at the implementation and local-verification stage.

The implementation permanently adds:

- a selectable constant-property reference backend;
- a Level 0 closed-form pressure reference;
- Level 1 quasi-steady pressure references;
- three committed reference configurations;
- deterministic thermochemistry metadata and dump output;
- CLI support for the reference backend.

Formal branch commit, Pull Request, merge, and branch cleanup remain
repository-integration actions after this completion record is committed.

## 2. Product implementation

### 2.1 Selectable backend

The following backend is now a supported explicit configuration choice:

    reference_constant

The corresponding configuration section is:

    thermochemistry.reference_constant

It requires:

    temperature_K
    gamma
    molecular_weight_kg_kmol

Unknown keys and nonphysical values are rejected.

The validated internal defaults are:

    temperature_K = 3000.0
    gamma = 1.2
    molecular_weight_kg_kmol = 25.0

Committed A-5 cases state these values explicitly.

### 2.2 Constant-property model

`ReferenceConstantBackend` provides:

    Y_gas = 1
    Y_condensed = 0
    R = R_UNIVERSAL / molecular_weight_kg_kmol
    T0 = temperature_K * eta_T0
    Theta = R * T0
    Psi = R * T0
    dTheta/dp = 0
    dPsi/dp = 0

The backend does not invoke NASA CEA or the legacy executable.

The backend does not create or reuse a thermochemistry cache.

It is an analytic verification backend and is not presented as a
combustion-chemistry prediction model.

### 2.3 Initialization and density

The backend is created by the normal `build_backend()` factory and is
registered through `initialize_thermochemistry()`.

The existing propellant-density contract remains unchanged:

    bulk mixture density =
        ideal additive-volume mixture density * packing fraction

No independent reference-only propellant-density input was added.

### 2.4 Thermochemistry export

The backend interface now provides:

    thermo_dump_table()

Existing pressure-table backends export their existing table columns.

The constant-property backend exports a deterministic two-point
representation at the configured pressure-range limits. The two points
contain identical constant properties.

This export is not a cache and is not used by the ODE solver.

### 2.5 CLI behavior

The CLI accepts:

    --backend reference_constant

The backend banner reports:

- constant analytic properties;
- effective temperature;
- gamma;
- molecular weight;
- absence of CEA use;
- disabled property caching;
- configured validity range.

`--dump-thermo` works with the constant-property backend.

When the CLI overrides a configuration that lacks an explicit reference
section, the documented internal defaults are used.

## 3. Level 0 reference

The Level 0 test uses:

- constant properties;
- a converging nozzle;
- no separation;
- no erosive burning;
- no throat erosion;
- no igniter;
- no post-burnout blowdown;
- `n_end = 0`;
- `n_burn = 0`;
- `a_burn = 1e-12 m/s`.

The existing product requirement that `a_burn` be positive was retained.

The small constant generation term is included in the exact solution:

    p(t) =
        p_inf
        + (p_initial - p_inf) * exp(-t/tau)

where:

    p_inf = rho_p * Ab * a_burn / C
    tau = V_g0 / (Psi * C)

The formal acceptance requirement is:

    maximum relative pressure error <= 1e-8

The permanent Level 0 test passed without relaxing this limit.

The test also confirmed that the ODE state retains exactly seven
components.

## 4. Level 1 quasi-steady references

The Level 1 cases use:

    n_burn = 0.01
    n_end = 0, 1, 2
    Ae/At = 1
    method = LSODA
    rtol = 1e-10

Flow separation, erosive burning, throat erosion, igniter injection, and
post-burnout blowdown are disabled.

For selected states in the central 50 percent of each numerical
trajectory, the test evaluates `equilibrium_pressure()` directly and
compares the numerical pressure with the total-mass-basis equilibrium
pressure.

The diagnostic maxima measured before fixing the acceptance criterion
were:

    n_end = 0:
        maximum relative error = 1.194670e-3

    n_end = 1:
        maximum relative error = 1.029688e-3

    n_end = 2:
        maximum relative error = 8.622144e-4

The fixed acceptance limit is:

    maximum relative error <= 1.5e-3

All three cases passed.

Every case also satisfied:

- successful integration;
- seven ODE states;
- `burnout` termination;
- at least ten valid comparison values;
- equilibrium pressures inside the configured property range.

## 5. Permanent reference configurations

The following files were added:

- `cases/reference/level1_end_faces_0.json`
- `cases/reference/level1_end_faces_1.json`
- `cases/reference/level1_end_faces_2.json`

Their final verified SHA-256 hashes are:

    b9427fa4d4d0834097d9bba80fa45b3a1c2ea06da201f89003a62177fea2518a
        cases/reference/level1_end_faces_0.json

    dd943bb66d64a0f02d0e379f86a3b82348c0ee03f5c02af95a60fcddef5d62fb
        cases/reference/level1_end_faces_1.json

    aafbf137a885b419a3f997b88b6ee8969b8813c2baad68a7e06368355a0ed43a
        cases/reference/level1_end_faces_2.json

All three files passed JSON parsing, configuration validation, backend
initialization, numerical execution, and quasi-steady comparison.

## 6. Compatibility

A-5 preserves the A-4 seven-component ODE state:

    [p0, x, Rt, m_out_total, Impulse, m_gen_total, m_igniter_total]

A-5 does not change:

- the meaning or order of existing states;
- total generated-product mass;
- integrated igniter mass;
- authoritative total discharge;
- sampled phase-discharge diagnostic status;
- energy-closure unavailability metadata;
- default scalar `rtol = 1e-10`;
- the total-system ledger limit of `1e-6`;
- the propellant ledger limit of `1e-9`;
- existing A-3 convergence values or reference values;
- existing A-3 result keys.

The protected directories remained unchanged:

- `verification/a3-convergence`
- `verification/a4-ledger`
- `Specification`

## 7. Automated verification

The final verification environment was:

    date:
        2026-10-06T22:48:32+09:00

    Python:
        3.14.4

    pytest:
        9.1.1

    local baseline commit:
        37748de3841f9f7ab0c297c0890689848ec05ba0

    origin/phase-a-foundation:
        37748de3841f9f7ab0c297c0890689848ec05ba0

The following checks completed successfully.

### 7.1 Python syntax

The checked files included:

- `GRIBS_v0.5.0-alpha.py`
- `tests/test_gribs.py`
- `tests/test_convergence_sweep.py`
- `tests/reference/test_reference_constant_backend.py`
- `tools/convergence_sweep.py`

Result:

    py_compile exit status = 0

### 7.2 Reference JSON validation

Results:

    cases/reference/level1_end_faces_0.json exit=0
    cases/reference/level1_end_faces_1.json exit=0
    cases/reference/level1_end_faces_2.json exit=0

### 7.3 A-5 tests

Result:

    22 passed in 3.87s
    exit status = 0

### 7.4 GRIBS self-tests

The existing CEA-backed self-test suite reported:

    all checks passed
    exit status = 0

The successful checks included the existing geometry, nozzle, property,
phase-fraction, derivative, pressure-factor, single-phase compatibility,
and mass-flow decomposition tests.

### 7.5 Complete regression suite

Result:

    75 passed in 18.72s
    exit status = 0

This consists of the preceding 53 tests plus 22 A-5 tests.

### 7.6 Whitespace validation

Result:

    git diff --check exit status = 0

## 8. Formal evidence

The retained A-5 evidence is stored under:

    verification/a5-reference/

The design record is:

    verification/a5-reference/A5_DESIGN.md

The completion record is:

    verification/a5-reference/A5_COMPLETION.md

The compact final evidence under `verification/a5-reference/final/`
contains:

- environment and Git state;
- syntax-check output and exit status;
- JSON validation results;
- A-5 pytest output and exit status;
- complete pytest output and exit status;
- GRIBS self-test output and exit status;
- whitespace-check output and exit status;
- reference-case hashes;
- A-5 collection count;
- protected-area status.

No large simulation histories, figures, or duplicate thermochemistry
caches are included in the A-5 evidence.

## 9. Definition of Done assessment

The A-5 Definition of Done is satisfied locally:

- `reference_constant` is selectable;
- invalid reference constants fail explicitly;
- reference tests require no CEA installation;
- the Level 0 analytic requirement passes at `1e-8`;
- all three Level 1 cases pass at `1.5e-3`;
- all permanent case files load and run;
- the seven-component A-4 state contract is preserved;
- deterministic provenance and thermochemistry export are available;
- 22 A-5 tests pass;
- all 75 tests pass;
- existing product self-tests pass;
- protected A-3 and A-4 evidence remains unchanged;
- design and final verification evidence are present.

Repository integration remains to be completed through commit, push,
Pull Request, merge into `phase-a-foundation`, synchronization, and
work-branch cleanup.
