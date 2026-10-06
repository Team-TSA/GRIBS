# Phase A-5 Reference Problems Design

## 1. Status

- Work package: A-5
- Title: Permanent Level 0/1 reference problems
- Base branch: `phase-a-foundation`
- Work branch: `phase-a-a5-reference-problems`
- Baseline commit: `37748de`
- A-4 merge commit: `37748de`
- Design status: implemented and under final verification

Phase A-0 through A-4 are complete. A-5 does not repeat or rewrite
those completed work packages.

## 2. Purpose

A-5 makes deterministic analytic and quasi-steady reference problems
part of the permanent GRIBS regression suite.

The work package provides:

1. a selectable constant-property thermochemistry backend;
2. a Level 0 analytic chamber-pressure reference;
3. a Level 1 quasi-steady equilibrium-pressure reference;
4. three committed Level 1 configurations for burning-end-face counts
   zero, one, and two;
5. deterministic thermochemistry provenance and dump output that require
   neither NASA CEA nor the legacy `fcea2` executable.

## 3. Authoritative inputs

The implementation is based on:

- `GRIBS_v0.5.0-alpha.py`
- `gribs_config.json`
- `tests/test_gribs.py`
- `tests/test_convergence_sweep.py`
- `Specification/level0_blowdown_check.py`
- `Specification/level0_blowdown_check.txt`
- `verification/a3-convergence/A3_DESIGN.md`
- `verification/a3-convergence/A3_COMPLETION.md`
- `verification/a4-ledger/A4_DESIGN.md`
- `verification/a4-ledger/A4_COMPLETION.md`
- `GRIBS改良計画_評価と詳細手順.md`

The historical Level 0 script is treated as a probe, not as production
code. Its analytic equation and demonstrated tolerance are promoted into
the current product test suite through the normal configuration and
backend-initialization paths.

## 4. Scope

A-5 includes:

- `thermochemistry.backend = "reference_constant"`;
- explicit constant temperature, gamma, and molecular weight;
- configuration validation and unknown-key rejection;
- normal backend-factory and initialization support;
- deterministic metadata;
- deterministic `--dump-thermo` output;
- CLI backend override support;
- Level 0 analytic verification;
- Level 1 quasi-steady verification;
- committed reference configurations for `n_end = 0, 1, 2`.

A-5 does not include:

- A-6 package restructuring;
- a new public `run_simulation()` API;
- A-7 JSON Schema freezing;
- a new combustion-chemistry model;
- changes to CEA equilibrium physics;
- the reserved nozzle `shock` transition;
- changes to A-3 sweep values or evidence;
- changes to A-4 ledger definitions or tolerances.

## 5. Reference backend contract

The selectable backend name is:

    reference_constant

Its required configuration block is:

    thermochemistry.reference_constant

with the keys:

    temperature_K
    gamma
    molecular_weight_kg_kmol

Validation requires:

    temperature_K > 0
    gamma > 1
    molecular_weight_kg_kmol > 0

Unknown keys are rejected.

The default values available through the internal `Config` contract and
the CLI backend override are:

    temperature_K = 3000.0
    gamma = 1.2
    molecular_weight_kg_kmol = 25.0

Committed reference configurations state these values explicitly.

## 6. Physical and numerical boundary

The backend represents a deterministic single-phase ideal gas:

    Y_gas = 1
    Y_condensed = 0

The gas constant is:

    R = R_UNIVERSAL / molecular_weight_kg_kmol

The effective constant temperature is:

    T0 = temperature_K * eta_T0

Therefore:

    Theta = R * T0
    Psi = Y_gas * R * T0 = R * T0
    dTheta/dp = 0
    dPsi/dp = 0

The backend performs no CEA calculation and creates no property cache.

It is an analytic verification backend, not a substitute combustion
chemistry model for engineering prediction.

## 7. Compatibility with A-4

A-5 preserves the seven-component ODE state:

    [p0, x, Rt, m_out_total, Impulse, m_gen_total, m_igniter_total]

No state is added, removed, reordered, or redefined.

A-5 also preserves:

- total-product mass accounting;
- the integrated igniter quadrature;
- `m_out_total` as the authoritative total discharge;
- sampled phase-discharge quantities as grid-dependent diagnostics;
- energy-closure unavailability metadata;
- scalar default `rtol = 1e-10`;
- the A-4 total-system limit of `1e-6`;
- the A-4 propellant limit of `1e-9`;
- all existing A-3 result keys.

## 8. Level 0 reference

The Level 0 case uses:

- constant gas properties;
- a converging nozzle with `Ae/At = 1`;
- no separation;
- no erosive burning;
- no throat erosion;
- no igniter;
- no post-burnout blowdown;
- `n_end = 0`;
- `n_burn = 0`;
- an effectively zero positive burn coefficient.

The existing product validation requires `a_burn > 0`. A-5 does not
change that contract. The test therefore uses:

    a_burn = 1e-12 m/s

The resulting constant generation term is retained exactly in the
closed-form solution.

For choked flow:

    mdot_out = C * p

where:

    C = Cd * At * sqrt(gamma / (R*T0))
        * (2/(gamma+1))^((gamma+1)/(2*(gamma-1)))

The pressure solution is:

    p(t) = p_inf + (p_initial - p_inf) * exp(-t/tau)

with:

    p_inf = rho_p * Ab * a_burn / C
    tau = V_g0 / (Psi * C)

Only points satisfying:

    p > 1.05 * p_choke

are compared, avoiding contamination from the nozzle transition.

The formal Level 0 acceptance limit is:

    max(abs(p_numeric / p_analytic - 1)) <= 1e-8

## 9. Level 1 quasi-steady reference

The Level 1 cases use:

- `n_burn = 0.01`;
- `n_end = 0, 1, 2`;
- a converging nozzle;
- separation disabled;
- erosive burning disabled;
- throat erosion disabled;
- igniter disabled;
- blowdown disabled;
- LSODA with `rtol = 1e-10`.

For each selected numerical state, the diagnostic equilibrium pressure
solves:

    mdot_gen_total(p_eq) = mdot_out_total(p_eq)

using the existing `equilibrium_pressure()` function on the A-4 total
mass basis.

The acceptance comparison uses direct equilibrium-pressure evaluations,
not the interpolated `sample()["p_eq"]` history.

The comparison interval is the central 50 percent of solver output
points, with at most 60 sampled states per case.

The measured pre-test diagnostic maxima were:

    n_end = 0: 1.194670e-3
    n_end = 1: 1.029688e-3
    n_end = 2: 8.622144e-4

The fixed acceptance limit is:

    maximum relative error <= 1.5e-3

The test also requires:

- successful integration;
- exactly seven ODE states;
- `burnout` termination;
- at least ten valid comparison values;
- all evaluated equilibrium pressures inside the configured property
  range.

## 10. Permanent reference configurations

The committed Level 1 configurations are:

- `cases/reference/level1_end_faces_0.json`
- `cases/reference/level1_end_faces_1.json`
- `cases/reference/level1_end_faces_2.json`

Their SHA-256 hashes at design review were:

    b9427fa4d4d0834097d9bba80fa45b3a1c2ea06da201f89003a62177fea2518a
    dd943bb66d64a0f02d0e379f86a3b82348c0ee03f5c02af95a60fcddef5d62fb
    aafbf137a885b419a3f997b88b6ee8969b8813c2baad68a7e06368355a0ed43a

Each configuration explicitly records the reference backend constants
and uses an independent output directory.

## 11. Thermochemistry dump and CLI behavior

`ThermochemistryBackend` exposes a deterministic
`thermo_dump_table()` contract.

Pressure-table backends return their existing table columns.

`ReferenceConstantBackend` returns a two-point representation at the
configured range limits. Both points carry identical constant values.
This is a deterministic export representation, not an interpolation
cache.

The CLI supports:

    --backend reference_constant

When overriding a configuration that does not contain the explicit
reference block, the documented `Config` defaults are used.

The CLI banner identifies:

- the constant analytic property model;
- the effective temperature;
- gamma;
- molecular weight;
- absence of CEA;
- disabled property caching;
- the configured validity range.

## 12. Tests

A-5 tests cover:

1. explicit backend selection;
2. required reference-property section;
3. unknown-key rejection;
4. temperature validation;
5. gamma validation;
6. molecular-weight validation;
7. backend-factory construction;
8. normal thermochemistry initialization;
9. mixture-density compatibility;
10. cache non-generation;
11. pressure-range diagnostics;
12. Level 0 analytic pressure;
13. seven-state ODE preservation;
14. Level 1 quasi-steady tracking for `n_end = 0, 1, 2`;
15. committed reference-config loading;
16. committed reference-config execution;
17. deterministic thermochemistry dumping;
18. accurate reference-backend banner output;
19. parameter documentation;
20. CLI backend override behavior.

Parameterization produces 22 collected A-5 tests.

## 13. Definition of Done

A-5 is complete when:

- `reference_constant` is selectable from configuration;
- invalid reference constants fail explicitly;
- no CEA installation is needed for reference tests;
- the Level 0 analytic test passes within `1e-8`;
- the Level 1 quasi-steady cases pass within `1.5e-3`;
- all three permanent case files load and run;
- the seven-component A-4 ODE contract remains intact;
- deterministic metadata and thermochemistry dumps are available;
- the A-5 test suite passes;
- the complete regression suite passes;
- A-3 and A-4 protected evidence remains unchanged;
- design, verification evidence, and completion records are committed.

## 14. Explicit non-goals

A-5 does not claim experimental validation.

A-5 does not make constant properties a production combustion-chemistry
recommendation.

A-5 does not alter:

- nozzle transition physics;
- two-phase particle dynamics;
- energy conservation modeling;
- CEA calculations;
- ledger tolerances;
- solver defaults;
- completed convergence evidence.
