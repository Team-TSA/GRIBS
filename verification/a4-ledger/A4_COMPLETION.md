# Phase A-4 Mass and Energy Ledger Completion

## 1. Completion status

Phase A-4 is complete.

The implementation adds an explicit product-mass ledger, a separate
propellant ledger, an integrated igniter-source quadrature, phase-discharge
diagnostics, and energy-availability metadata.

The final production verification passed for:

- the distributed LSODA configuration without an active igniter;
- an active-igniter LSODA configuration;
- an active-igniter BDF configuration;
- an active-igniter Radau configuration.

The complete automated regression suite passed with 53 tests.

## 2. Product changes

### 2.1 Mass ledger

`results.ledger.mass` now reports:

- initial chamber product inventory;
- total generated product mass;
- generated propellant-derived mass;
- commanded igniter charge;
- integrated igniter mass;
- total discharged mass;
- final chamber product inventory;
- initial solid-propellant mass;
- final unburned propellant mass;
- geometry-derived burned propellant mass;
- total-system signed and relative residuals;
- propellant signed and relative residuals;
- a sampled phase-discharge diagnostic.

The formal acceptance limits are:

    total-system relative residual <= 1e-6

    propellant relative residual <= 1e-9

### 2.2 Integrated igniter bookkeeping

The ODE state was extended from six to seven components by appending:

    m_igniter_total

Its derivative is the existing igniter source:

    dm_igniter_total/dt = igniter_mdot(t)

Existing state indices retain their original meanings.

The product ledger uses the integrated state for injected igniter mass.
The configured igniter charge remains separately available as the
commanded mass.

The history CSV exposes:

    m_igniter[kg]

### 2.3 Energy metadata

GRIBS does not integrate an independent energy state.

The product result therefore reports explicitly that:

- an energy closure residual is unavailable;
- the HP equilibrium constraint is not a complete energy-conservation law;
- the state relation `p0 * Vg = m_total * Psi` is not an independent
  energy closure.

No artificial energy residual is reported.

### 2.4 Phase-discharge diagnostic

Gas, condensed, and total discharged masses are estimated from sampled
history by trapezoidal quadrature.

These values are explicitly marked:

    method = sampled_history_trapezoid

    grid_dependent = true

The total discharged ODE state remains the authoritative total-discharge
quantity.

## 3. Compatibility

The following A-3 result keys remain available:

- `p_max_Pa`;
- `burn_time_s`;
- `F_max_N`;
- `impulse_total_Ns`;
- `total_mass_consistency_error`;
- `propellant_mass_balance_error`;
- burn and blowdown stop reasons.

The legacy propellant metric is retained for compatibility, but it
includes all generated products and therefore includes igniter mass.

The new igniter-separated result is authoritative for the A-4 propellant
ledger.

The preserved Phase A-3 evidence was not modified.

## 4. Default solver tolerance

The scalar product default was changed from:

    rtol = 1e-9

to:

    rtol = 1e-10

The distributed JSON configuration was changed to the same value.

The ledger acceptance limits were not changed.

This change was selected because:

1. the active-igniter LSODA case failed the `1e-9` propellant-ledger
   criterion at scalar `rtol=1e-9`;
2. LSODA passed at scalar `rtol=1e-10`;
3. BDF and Radau passed at scalar `rtol=1e-10`;
4. the existing A-3 sweep already included `rtol=1e-10`;
5. the `rtol=1e-10` A-3 convergence row passed against the `1e-11`
   reference;
6. the 20-case candidate sweep completed successfully;
7. principal-output differences remained below existing A-3 relative
   acceptance limits for the candidate row.

## 5. Rejected alternatives

The implementation does not:

- relax the `1e-9` propellant-ledger criterion;
- define the propellant residual to be zero by construction;
- hide failed residuals;
- use a solver-specific vector `rtol`;
- retain the experimental igniter-cutoff segmented integration;
- claim an energy closure that the model does not calculate.

The complete investigation, including rejected prototypes and large raw
histories, was archived outside the repository before formal evidence was
pruned.

## 6. Final production verification

All four final configurations passed.

### 6.1 Distributed default, LSODA

    rtol = 1e-10
    total relative residual =
        4.9247133432102780e-12
    propellant relative residual =
        1.4721688919452154e-10

### 6.2 Active igniter, LSODA

    rtol = 1e-10
    total relative residual =
        8.5825276658719562e-13
    propellant relative residual =
        3.4860116215021564e-10
    integrated igniter mass =
        1.9999999999930425e-3 kg

### 6.3 Active igniter, BDF

    rtol = 1e-10
    total relative residual =
        1.8181139647975681e-11
    propellant relative residual =
        5.2002105521511501e-10
    integrated igniter mass =
        2.0000000000016433e-3 kg

### 6.4 Active igniter, Radau

    rtol = 1e-10
    total relative residual =
        1.2662960502864900e-12
    propellant relative residual =
        2.3117631090670790e-14
    integrated igniter mass =
        2.0000000000122507e-3 kg

The final machine-readable result is:

    verification/a4-ledger/final-production/
    final_production_verification.json

Its overall result is:

    overall_passed = true

## 7. Automated verification

The following completed successfully:

    python -m py_compile \
      GRIBS_v0.5.0-alpha.py \
      tests/test_gribs.py \
      tests/test_convergence_sweep.py \
      tools/convergence_sweep.py

    git diff --check

    python GRIBS_v0.5.0-alpha.py --selftest

    pytest -q

Final pytest result:

    53 passed

## 8. Formal evidence

The formal A-4 evidence retained in the repository contains:

- the design and decision history;
- compact convergence reports;
- the final four production summaries;
- machine-readable final verification;
- configuration files;
- validation and execution exit codes;
- final self-test and pytest results.

Large time-history CSV files, figures, duplicated thermochemistry caches,
source backups, and intermediate failure logs are excluded from the
formal repository evidence.
