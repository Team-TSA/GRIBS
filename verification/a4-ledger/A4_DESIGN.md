# Phase A-4 Mass and Energy Ledger Design

## 1. Status

- Work package: A-4
- Title: Mass and energy ledger
- Base branch: `phase-a-foundation`
- Work branch: `phase-a-a4-mass-energy-ledger`
- Baseline commit: `f12be54`
- A-3 implementation commit: `eb52158`
- Design status: implementation boundary fixed
- Product implementation status: not started

A-0, A-1, A-2, and A-3 are complete. A-4 must not repeat or rewrite
those completed work packages.

## 2. Purpose

A-4 promotes the existing mass-ledger investigation into a product
diagnostic of GRIBS v0.5.0-alpha.

The ledger shall:

1. itemise the total-system mass balance;
2. distinguish propellant-derived mass from igniter-derived mass;
3. report a physically meaningful propellant-only residual;
4. preserve the existing A-3 convergence metrics for compatibility;
5. expose the ledger in machine-readable and human-readable outputs;
6. state honestly that a physical energy-closure residual is unavailable
   because GRIBS does not integrate an independent energy state.

## 3. Authoritative inputs

The implementation shall be based on:

- `GRIBS_v0.5.0-alpha.py`
- `gribs_config.json`
- `tests/test_gribs.py`
- `tools/convergence_sweep.py`
- `tools/convergence_sweep_config.json`
- `tests/test_convergence_sweep.py`
- `Specification/mass_ledger.py`
- `Specification/mass_ledger.txt`
- `GRIBS改良計画_評価と詳細手順.md`
- `verification/a3-convergence/A3_DESIGN.md`
- `verification/a3-convergence/A3_COMPLETION.md`

The investigation script is a prototype and shall not be copied
mechanically. Its equations and demonstrated limitations shall be
incorporated into the current v0.5.0-alpha engine.

## 4. Existing model facts

The ODE state is:

    y = [p0, x, Rt, m_out_total, impulse, m_gen_total]

The total generated-product rate is:

    mdot_gen_total = rho_p * r * Ab + igniter_mdot(t)

The total discharged-product rate is:

    mdot_out_total

Both integrated mass states are on the total-product basis, including
gas and condensed products in homogeneous-equilibrium mode.

The initial chamber inventory is:

    m_chamber_initial = p0_init * V_g0 / Psi(p0_init)

where:

    Psi = Y_gas * R_gas * T0

The final chamber inventory is evaluated from the EOS:

    m_chamber_final = p0_final * V_g_final / Psi(p0_final)

## 5. Existing inconsistency to correct

The current result named `propellant_mass_balance_error` is:

    abs(m_gen_total / m_propellant_initial - 1)

The numerator includes igniter-derived mass. Therefore the value is not
a propellant-only balance whenever the igniter is active.

This existing key shall remain temporarily for A-3 compatibility, but
it shall be marked as a legacy diagnostic in the A-4 output and shall
not be used as the new propellant-ledger acceptance metric.

## 6. Mass-ledger definitions

### 6.1 Primary masses

Let:

- `A = initial_chamber_inventory_kg`
- `B = generated_mass_total_kg`
- `B1 = generated_mass_propellant_kg`
- `B2 = injected_mass_igniter_kg`
- `C = discharged_mass_total_kg`
- `D = final_chamber_inventory_kg`
- `E = final_unburned_propellant_kg`
- `F = burned_propellant_geometry_kg`
- `P = initial_propellant_mass_kg`

Definitions:

    B = final ODE state m_gen_total
    C = final ODE state m_out_total
    D = final p0 * Vg / Psi
    P = propellant_mass(c)
    E = remaining solid propellant from final burn depth and geometry
    F = P - E
    B2 = actual igniter mass injected over the simulation interval
    B1 = B - B2

### 6.2 Total-system closure

The signed total-system residual is:

    r_total_kg = A + B - C - D

The normalized absolute residual is:

    epsilon_total = abs(r_total_kg) / max(abs(A + B), tiny)

The A-4 default acceptance limit is:

    epsilon_total <= 1e-6

### 6.3 Propellant-only closure

The signed propellant residual is:

    r_propellant_kg = B1 - F

The normalized absolute residual is:

    epsilon_propellant = abs(r_propellant_kg) / max(abs(P), tiny)

The A-4 default acceptance limit is:

    epsilon_propellant <= 1e-9

Both the signed residual and normalized absolute residual shall be
reported. Acceptance uses the normalized absolute value.

## 7. Igniter accounting

The ledger shall distinguish:

- commanded igniter charge:

      igniter_mass_commanded_kg = ign_mdot * ign_time

- actual injected igniter mass:

      igniter_mass_injected_kg = mass admitted by the model over the
      integrated simulation interval

The actual value, not the commanded value, shall be used to separate
propellant-derived and igniter-derived generated mass.

The implementation shall also make the sampled `mdot_gen_total` history
consistent with the ODE definition during every phase. It shall not report
zero during blowdown if the existing ODE is still adding igniter flow.

This is an output-consistency correction, not a new physical model.

## 8. Phase-resolved discharge diagnostics

The ODE contains only one integrated discharge state:

    m_out_total

Gas and condensed discharge rates are available in the sampled history.
Their accumulated masses may be reported using history-grid quadrature,
but shall be labelled:

- diagnostic only;
- sampled-history quadrature;
- grid dependent;
- not ODE states.

The ledger shall report the quadrature mismatch against the ODE total
discharge state when phase-resolved discharge masses are included.

No new gas-discharge or condensed-discharge ODE state is added in A-4.

## 9. Existing metric compatibility

The following existing result keys shall remain unchanged in A-4:

- `total_mass_consistency_error`
- `propellant_mass_balance_error`

`total_mass_consistency_error` means the maximum history-wide relative
difference between EOS chamber mass and balance chamber mass. It is not
the same as the new terminal total-system ledger residual.

`propellant_mass_balance_error` retains its existing numerical definition
temporarily so that the completed A-3 evidence remains comparable. The
new correct propellant metric shall be stored under the structured ledger.

The existing keys shall not silently change meaning.

## 10. Proposed result structure

The `summarize()` result dictionary shall contain the following hierarchy:

    ledger:
      schema_version
      mass:
        basis
        initial_chamber_inventory_kg
        generated_mass_total_kg
        generated_mass_propellant_kg
        igniter_mass_commanded_kg
        igniter_mass_injected_kg
        discharged_mass_total_kg
        final_chamber_inventory_kg
        initial_propellant_mass_kg
        final_unburned_propellant_kg
        burned_propellant_geometry_kg
        total_system:
          signed_residual_kg
          relative_residual
          tolerance
          passed
        propellant:
          signed_residual_kg
          relative_residual
          tolerance
          passed
        phase_discharge_diagnostic:
          method
          grid_dependent
          gas_mass_kg
          condensed_mass_kg
          total_quadrature_mass_kg
          ode_total_mass_kg
          relative_quadrature_mismatch
        passed
      energy:
        available
        closure_residual_available
        reason
        model
      passed

The exact spelling shall be fixed by tests before product integration.
A-4 shall not create the final A-7 JSON Schema.

## 11. Output locations

The ledger shall appear in:

1. `summary.json` under `results.ledger`;
2. `summary.txt` in a dedicated mass-ledger section;
3. the in-memory dictionary returned by `summarize()`.

No separate production ledger text file is required by A-4.

## 12. Acceptance behavior

A-4 shall report:

- each residual value;
- each tolerance;
- each individual pass/fail result;
- the combined mass-ledger pass/fail result;
- the combined ledger pass/fail result.

The initial A-4 implementation shall not add a new user-configurable
enforcement section to `gribs_config.json`.

Reasons:

1. the current parser strictly rejects unknown configuration keys;
2. configuration-contract freezing belongs to A-7;
3. aborting before output would hide the failed ledger from `summary.json`;
4. the existing post-run conservation tests continue to protect normal
   execution.

In A-4, a tolerance violation shall therefore be represented as a
machine-readable `passed: false`. Existing post-run conservation-test
failure behavior remains unchanged.

## 13. Energy-ledger boundary

GRIBS does not integrate an independent energy state. Therefore A-4 shall
not claim or compute a physical energy-closure residual.

The energy section shall state:

- `available: false`;
- `closure_residual_available: false`;
- the reason that no independent energy state is integrated;
- the existing HP-equilibrium constraint;
- the reactant-enthalpy treatment;
- the temperature-efficiency treatment;
- the definition `Psi = Yg * Rg * T0`;
- that `p * Vg = m_total * Psi` is a state-equation consistency relation,
  not a complete energy conservation law.

A-4 shall not introduce:

- a total-energy ODE;
- wall heat transfer;
- condensed-phase internal-energy dynamics;
- particle thermal lag;
- a new CEA formulation.

## 14. A-3 relationship

A-3 remains complete and its formal evidence shall not be overwritten.

A-4 shall initially preserve the A-3 metric names and behavior. After
the product ledger is implemented and tested, the convergence harness
may be extended to record the new ledger residuals while retaining the
old columns for evidence continuity.

Any new A-4 convergence evidence shall be saved separately from:

    verification/a3-convergence/latest/

## 15. Tests required before completion

Fast tests shall cover:

1. total-system ledger equations;
2. propellant-only ledger equations;
3. zero-igniter behavior;
4. active-igniter separation;
5. commanded versus actually injected igniter mass;
6. early termination before the commanded igniter duration;
7. geometry-based unburned and burned propellant masses;
8. total-system pass/fail threshold;
9. propellant pass/fail threshold;
10. signed residual preservation;
11. structured JSON-compatible ledger output;
12. legacy result-key preservation;
13. energy-unavailable metadata;
14. phase-discharge diagnostic labelling;
15. sampled generation-rate consistency across phases.

Integration tests shall cover:

1. the default configuration;
2. an active-igniter configuration matching the preserved probe;
3. `summary.json` placement under `results.ledger`;
4. `summary.txt` ledger rendering;
5. the complete pytest suite;
6. the A-3 convergence harness.

## 16. Definition of Done

A-4 is complete when:

- the mass ledger is generated on every successful production run;
- the total-system residual is reported;
- the propellant-only residual is reported;
- igniter mass is not counted as propellant error;
- tolerances and pass/fail results are machine-readable;
- the default case meets the total-system limit of `1e-6`;
- the active-igniter reference case meets the propellant limit of `1e-9`;
- legacy A-3 metrics remain available and unchanged;
- the energy section accurately records that closure is unavailable;
- no unsupported energy model is added;
- all A-4 tests pass;
- the complete regression suite passes;
- A-3 convergence behavior remains operational;
- design, verification evidence, and completion records are committed.

## 17. Explicit non-goals

A-4 does not implement:

- package restructuring from A-6;
- final JSON Schema freezing from A-7;
- final public Result API freezing;
- new gas and condensed discharge ODE states;
- particle slip or deposition;
- wall heat loss;
- a complete energy equation;
- changes to CEA equilibrium physics;
- the reserved shock-transition nozzle model;
- rewriting or replacing A-3 evidence.

## 18. Igniter quadrature amendment

### 18.1 Reason for the amendment

Production verification with an active igniter found:

- baseline solver settings:
  propellant relative residual = `1.8624912520418528e-9`;
- `rtol = 1e-10`:
  propellant relative residual = `3.4848522723029594e-10`;
- `max_step_burn = 0.005 s`:
  propellant relative residual = `1.2053283180259971e-9`;
- combined tighter settings:
  propellant relative residual = `2.1645326960583694e-10`.

The residual decreases by a factor of approximately 8.60 under tighter
integration settings. This demonstrates that the baseline failure is
dominated by numerical integration around the discontinuous igniter
cutoff, rather than by an incorrect physical mass balance.

The investigation evidence is stored in:

- `verification/a4-ledger/igniter-convergence/igniter_convergence.csv`
- `verification/a4-ledger/igniter-convergence/igniter_convergence.md`

### 18.2 Rejected responses

A-4 shall not resolve the baseline failure by:

- relaxing the `1e-9` propellant acceptance limit;
- changing the global default solver tolerance solely for the ledger;
- hiding or rounding the failed residual;
- using the commanded igniter charge as if it were an integrated state;
- redefining the propellant residual to make it pass automatically.

### 18.3 Revised numerical accounting

The ODE shall add one bookkeeping quadrature:

    m_igniter_total

Its derivative is:

    dm_igniter_total/dt = igniter_mdot(t)

The total generated-product quadrature remains:

    dm_gen_total/dt =
        rho_p * r * Ab + igniter_mdot(t)

The product ledger shall use:

    generated_mass_total_kg = m_gen_total

    igniter_mass_injected_kg = m_igniter_total

    generated_mass_propellant_kg =
        m_gen_total - m_igniter_total

The commanded charge remains a separate configuration-derived value:

    igniter_mass_commanded_kg = ign_mdot * ign_time

The new igniter quadrature is numerical bookkeeping of an existing mass
source. It does not introduce a new physical model.

### 18.4 Compatibility

Existing ODE-state indices shall retain their current meanings. The new
igniter state shall be appended after the existing total-generated-mass
state.

The existing history columns and A-3 metrics shall remain available.

A new history column may expose the integrated igniter mass, but existing
column names shall not be renamed or redefined.

The helper `igniter_injected_mass(c, simulation_end_s)` may remain as the
analytic commanded-over-simulated-interval diagnostic, but it shall not
be the authoritative source for the product ledger after the ODE
quadrature is available.

### 18.5 Additional acceptance tests

The implementation shall demonstrate that:

1. the igniter ODE state integrates the active igniter interval;
2. zero igniter input produces zero integrated igniter mass;
3. the integrated igniter mass is used by the ledger;
4. the active-igniter production case passes the `1e-9` propellant limit
   at the default solver settings;
5. existing total-generation, discharge, impulse, and A-3 outputs remain
   compatible;
6. the complete regression suite passes.

## 19. Result of the independent igniter quadrature

### 19.1 Production result

The independent igniter quadrature was implemented and verified.

The active-igniter production case at the default solver settings
reported:

    commanded igniter mass:
        2.0000000000000000e-3 kg

    integrated igniter mass:
        1.9999999999960089e-3 kg

The integrated value was exposed as `m_igniter[kg]` in the history CSV.
The CSV final value and `results.ledger.mass.igniter_mass_injected_kg`
were consistent.

The complete regression suite passed:

    51 passed

### 19.2 Rejected hypothesis

The independent igniter quadrature did not reduce the propellant-ledger
residual.

Before the quadrature:

    propellant relative residual =
        1.8624912520418528e-9

After the quadrature:

    propellant relative residual =
        1.8625577152312386e-9

The ratio of the previous value to the new value was approximately:

    0.999964316

Therefore, the baseline residual is not primarily caused by subtracting
an analytic igniter charge from the total generated-mass quadrature.

### 19.3 Revised interpretation

The remaining propellant residual compares two independently evaluated
numerical paths:

    integrated propellant generation =
        integral(rho_p * r * Ab dt)

and:

    geometry-derived burned propellant =
        initial propellant mass
        - final unburned propellant mass at x_final

These paths are analytically equivalent through:

    dx/dt = r

    dV_burned/dx = Ab

    dm_propellant/dt = rho_p * r * Ab

The observed residual is therefore a numerical consistency error between
the burn-depth integration and generated-mass integration.

The active igniter introduces a discontinuity at:

    t = ign_time

The previous convergence investigation showed that reducing `rtol` from
`1e-9` to `1e-10` reduced the propellant residual below its `1e-9`
acceptance limit. This supports a numerical-integration origin.

### 19.4 Retention of the igniter quadrature

The `m_igniter_total` ODE state shall be retained because it provides:

1. direct accounting of actually integrated igniter mass;
2. separation of commanded and injected igniter mass;
3. consistency between the ODE state, history CSV, and product ledger;
4. correct accounting if the simulation ends before the commanded
   igniter duration;
5. a basis for future time-dependent igniter source models.

The state is bookkeeping for an existing source term and does not change
the physical model.

### 19.5 Next corrective action

The next corrective action shall be to treat the igniter cutoff as an
explicit integration boundary.

For an active igniter with:

    0 < ign_time < burning-phase end time

the burning-phase integration shall be divided into:

1. an igniter-active interval ending exactly at `ign_time`;
2. an igniter-inactive interval beginning from the state at
   `ign_time` and continuing to the normal burning-phase terminal event.

This avoids integrating across the discontinuous igniter source within
one adaptive solver step.

The implementation shall preserve:

- the existing physical model;
- the existing solver defaults;
- all existing ODE-state meanings;
- event detection for burnout, choking loss, recovery, and ambient
  pressure;
- dense-output sampling;
- existing history and A-3 result keys.

### 19.6 Actions still rejected

A-4 shall not resolve this issue by:

- relaxing the `1e-9` propellant residual limit;
- changing the default global `rtol` solely for the ledger;
- defining generated propellant mass to equal the geometry result;
- suppressing the failed acceptance result;
- rounding the residual before comparison.

## 19. Result of the independent igniter quadrature

### 19.1 Production result

The independent igniter quadrature was implemented and verified.

The active-igniter production case at the default solver settings
reported:

    commanded igniter mass:
        2.0000000000000000e-3 kg

    integrated igniter mass:
        1.9999999999960089e-3 kg

The integrated value was exposed as `m_igniter[kg]` in the history CSV.
The CSV final value and `results.ledger.mass.igniter_mass_injected_kg`
were consistent.

The complete regression suite passed:

    51 passed

### 19.2 Rejected hypothesis

The independent igniter quadrature did not reduce the propellant-ledger
residual.

Before the quadrature:

    propellant relative residual =
        1.8624912520418528e-9

After the quadrature:

    propellant relative residual =
        1.8625577152312386e-9

The ratio of the previous value to the new value was approximately:

    0.999964316

Therefore, the baseline residual is not primarily caused by subtracting
an analytic igniter charge from the total generated-mass quadrature.

### 19.3 Revised interpretation

The remaining propellant residual compares two independently evaluated
numerical paths:

    integrated propellant generation =
        integral(rho_p * r * Ab dt)

and:

    geometry-derived burned propellant =
        initial propellant mass
        - final unburned propellant mass at x_final

These paths are analytically equivalent through:

    dx/dt = r

    dV_burned/dx = Ab

    dm_propellant/dt = rho_p * r * Ab

The observed residual is therefore a numerical consistency error between
the burn-depth integration and generated-mass integration.

The active igniter introduces a discontinuity at:

    t = ign_time

The previous convergence investigation showed that reducing `rtol` from
`1e-9` to `1e-10` reduced the propellant residual below its `1e-9`
acceptance limit. This supports a numerical-integration origin.

### 19.4 Retention of the igniter quadrature

The `m_igniter_total` ODE state shall be retained because it provides:

1. direct accounting of actually integrated igniter mass;
2. separation of commanded and injected igniter mass;
3. consistency between the ODE state, history CSV, and product ledger;
4. correct accounting if the simulation ends before the commanded
   igniter duration;
5. a basis for future time-dependent igniter source models.

The state is bookkeeping for an existing source term and does not change
the physical model.

### 19.5 Next corrective action

The next corrective action shall be to treat the igniter cutoff as an
explicit integration boundary.

For an active igniter with:

    0 < ign_time < burning-phase end time

the burning-phase integration shall be divided into:

1. an igniter-active interval ending exactly at `ign_time`;
2. an igniter-inactive interval beginning from the state at
   `ign_time` and continuing to the normal burning-phase terminal event.

This avoids integrating across the discontinuous igniter source within
one adaptive solver step.

The implementation shall preserve:

- the existing physical model;
- the existing solver defaults;
- all existing ODE-state meanings;
- event detection for burnout, choking loss, recovery, and ambient
  pressure;
- dense-output sampling;
- existing history and A-3 result keys.

### 19.6 Actions still rejected

A-4 shall not resolve this issue by:

- relaxing the `1e-9` propellant residual limit;
- changing the default global `rtol` solely for the ledger;
- defining generated propellant mass to equal the geometry result;
- suppressing the failed acceptance result;
- rounding the residual before comparison.

## 20. Igniter-cutoff integration-boundary design

### 20.1 Selected representation

The burning solution shall be represented internally as an ordered list
of `OdeResult` segments:

    sol_burn_segments = [
        igniter-active segment,
        igniter-inactive segment,
    ]

When no split is required, the list shall contain exactly one segment.

The existing `sol_burn` result key shall remain available for
transitional compatibility:

- for an unsplit burn, `sol_burn` is the only segment;
- for a split burn, `sol_burn` is the final burning segment.

All production post-processing shall use `sol_burn_segments` when it is
available. Old synthetic or external result dictionaries containing only
`sol_burn` shall remain accepted by `sample()`.

### 20.2 Split eligibility

A split shall be attempted only when:

    ign_mdot > 0
    ign_time > 0
    ign_time < t_max

The first burning segment shall integrate from:

    t = 0

to:

    t = ign_time

using the normal burning RHS and all normal burning events.

If a terminal physical event occurs before `ign_time`, no second burning
segment shall be created.

If the first segment reaches `ign_time` without a terminal physical
event, the second segment shall start from the first segment's final
state and integrate from:

    t = ign_time

to:

    t = t_max

using the normal burning RHS and normal burning events.

The igniter source is already zero for `t > ign_time`.

### 20.3 Boundary ownership

The state at `t = ign_time` belongs to the end of the first segment.

When histories are combined, the duplicate first sample of every segment
after the first shall be removed. Therefore the cutoff time shall appear
exactly once in the sampled history.

### 20.4 Event aggregation

The product result shall aggregate each named burning event over all
burning segments:

- burnout;
- unchoked;
- rechoked;
- ambient.

Event times shall remain absolute simulation times.

The final burning stop reason shall be determined from the segment that
terminates the burning phase, not from the artificial igniter boundary.

The igniter boundary itself is not a physical stop reason and shall not
replace `burnout`, `ambient`, `choking-loss`, or `time-limit`.

### 20.5 Blowdown handoff

Blowdown shall start from the final state and final time of the last
burning segment.

The existing blowdown conditions remain unchanged:

- blowdown is enabled;
- the burning phase ended by burnout.

### 20.6 Compatibility

The following shall remain unchanged:

- ODE state meanings and ordering;
- solver defaults;
- existing history-column meanings;
- existing summary keys;
- A-3 convergence-harness inputs and output columns;
- physical event definitions;
- post-run conservation tests.

No composite replacement for SciPy `OdeResult` shall be introduced.

## 21. Result of the igniter-cutoff integration boundary

### 21.1 Production result

The active-igniter production case was rerun with the burning integration
split exactly at:

    ign_time = 0.1 s

The result was:

    propellant relative residual before split:
        1.8625577152312386e-9

    propellant relative residual after split:
        1.6737793709211319e-9

    acceptance limit:
        1.0e-9

The split improved the residual by approximately 10.1 percent, but the
product ledger still failed its required acceptance limit.

All static checks, self-tests, A-4 fast tests, and the complete regression
suite passed. The production ledger alone remained failed.

### 21.2 Interpretation

The igniter cutoff contributes to the numerical error, but it is not the
dominant or sufficient cause.

An explicit integration boundary does not make the default solver
configuration satisfy the `1e-9` propellant consistency requirement.

The remaining difference continues to be the numerical consistency error
between:

    integral(rho_p * r * Ab dt)

and:

    geometry-derived burned propellant from x_final

### 21.3 Boundary-sample observation

The production inspection counted nine sampled times within `1e-12 s` of
the cutoff.

This count does not mean that nine exactly identical boundary rows were
written. It includes distinct adaptive and generated-grid values close to
`0.1 s`.

The segment-combination unit test confirmed that the exact first sample
of a later segment was removed. Nevertheless, the additional segmented
sampling behavior adds complexity without satisfying the A-4 requirement.

### 21.4 Decision

The igniter-cutoff burning-phase split shall be rolled back.

Reasons:

1. the product acceptance limit remains failed;
2. the numerical improvement is insufficient;
3. the change adds event aggregation and multi-segment lifecycle logic;
4. the added complexity is not justified by the observed result;
5. A-4 shall not accumulate speculative numerical changes.

The independent `m_igniter_total` bookkeeping quadrature shall remain.

The next investigation shall evaluate whether stricter numerical control
can be applied specifically to bookkeeping mass states without changing
the global pressure and geometry solver tolerance.

The global default `rtol` shall not be changed before that investigation.

## 22. Result of the state-wise relative-tolerance probe

### 22.1 Environment capability

The installed numerical environment reported:

    SciPy version: 1.18.1

The tested solver interfaces behaved as follows:

- LSODA accepted a state-wise relative-tolerance array;
- BDF rejected a state-wise relative-tolerance array;
- Radau rejected a state-wise relative-tolerance array.

Therefore, a vector `rtol` is not a solver-independent GRIBS feature.

### 22.2 Mass-state scale analysis

At the active-igniter final state:

    m_gen = approximately 6.2031e-2 kg
    m_igniter = approximately 2.0e-3 kg
    configured scalar rtol = 1e-9
    mass-state atol = 1e-12 kg

For `m_gen`, the relative tolerance contribution was approximately:

    6.2031e-11 kg

This was approximately 62 times the absolute tolerance.

For `m_igniter`, the relative tolerance contribution was approximately:

    2.0e-12 kg

This was approximately twice the absolute tolerance.

Tightening only the mass-state absolute tolerance is therefore not
expected to resolve the generated-propellant consistency residual.

### 22.3 State-wise probe results

The active-igniter case produced:

    scalar rtol 1e-9:
        propellant relative residual = 1.8625577152312386e-9
        FAIL

    x rtol 1e-10:
        propellant relative residual = 3.0487104206201458e-10
        PASS

    m_gen rtol 1e-10:
        propellant relative residual = 3.3494673334651093e-10
        PASS

    x and m_gen rtol 1e-10:
        propellant relative residual = 3.2725214542613659e-10
        PASS

    x, m_gen, and m_igniter rtol 1e-10:
        propellant relative residual = 3.2725214542613659e-10
        PASS

    all states rtol 1e-10:
        propellant relative residual = 3.4860116215021564e-10
        PASS

Tightening either the burn-depth state or the generated-mass state was
sufficient to meet the current `1e-9` ledger criterion in the tested
LSODA case.

### 22.4 Effect on principal outputs

Relative to the scalar `1e-9` reference, the tested state-wise settings
changed the principal outputs only at small numerical levels:

- maximum pressure: no more than approximately `3.6e-12` relative;
- burn time: approximately `4e-10` to `6e-10` relative;
- total impulse: approximately `1.6e-9` relative.

These differences are numerical, not model changes.

### 22.5 Decision

State-wise `rtol` shall not yet be added to the product implementation.

Reasons:

1. it is accepted by LSODA but rejected by BDF and Radau;
2. GRIBS currently exposes all three methods through one solver contract;
3. silently applying different error-control semantics by method would
   weaken reproducibility;
4. the A-3 convergence policy was established with scalar tolerances;
5. solver-independent behavior must be investigated before changing the
   product configuration.

The next investigation shall compare scalar `rtol = 1e-9` and
`rtol = 1e-10` for LSODA, BDF, and Radau.

## 23. Solver-method convergence result

### 23.1 Tested matrix

The active-igniter production case was evaluated with the three supported
solver methods at scalar relative tolerances of `1e-9` and `1e-10`.

The results were:

    LSODA, rtol = 1e-9:
        propellant relative residual = 1.8625577152312386e-9
        FAIL

    LSODA, rtol = 1e-10:
        propellant relative residual = 3.4860116215021564e-10
        PASS

    BDF, rtol = 1e-9:
        propellant relative residual = 7.8151498341568308e-10
        PASS

    BDF, rtol = 1e-10:
        propellant relative residual = 5.2002105521511501e-10
        PASS

    Radau, rtol = 1e-9:
        propellant relative residual = 4.8200260824048598e-14
        PASS

    Radau, rtol = 1e-10:
        propellant relative residual = 2.3117631090670790e-14
        PASS

All six calculations completed successfully. All six total-system mass
residuals passed the `1e-6` criterion.

### 23.2 Interpretation

The strict propellant-ledger failure is specific to the tested combination
of:

    method = LSODA
    scalar rtol = 1e-9

BDF and Radau met the `1e-9` propellant criterion at the existing scalar
`rtol = 1e-9`.

LSODA met the criterion when its scalar tolerance was reduced to `1e-10`.

Therefore:

1. the ledger equations are not generally inconsistent;
2. the active-igniter case is not generally unsolvable at the required
   precision;
3. a solver-independent scalar `rtol = 1e-10` is a viable numerical
   candidate;
4. a state-wise tolerance is unnecessary to demonstrate compliance;
5. the existing scalar solver contract can be preserved.

### 23.3 Timing observation

The single measured wall times were approximately:

    LSODA, rtol 1e-9:  5.753 s
    LSODA, rtol 1e-10: 5.673 s
    BDF, rtol 1e-9:    6.039 s
    BDF, rtol 1e-10:   6.149 s
    Radau, rtol 1e-9:  6.639 s
    Radau, rtol 1e-10: 6.781 s

These are individual wall-clock observations and shall not be interpreted
as a performance benchmark. They do not show a prohibitive cost for the
tested `1e-10` settings.

### 23.4 Candidate decision

The next candidate is a solver-independent scalar relative tolerance of:

    rtol = 1e-10

The candidate shall not be adopted until its relationship with the
completed A-3 convergence contract has been checked.

Before changing the example or default configuration, A-4 shall verify:

1. how the A-3 harness selects its reference tolerance;
2. whether `1e-9` is an A-3 reference case, a product default, or both;
3. whether changing the product example would invalidate preserved A-3
   evidence;
4. whether the A-3 sweep already includes `1e-10`;
5. whether a new A-4 validation profile can use `1e-10` without rewriting
   A-3 evidence.

The `1e-9` ledger acceptance limit shall remain unchanged.

## 24. Review of the A-3 candidate sweep acceptance

### 24.1 Candidate result

The A-4 candidate sweep used a base configuration with:

    method = LSODA
    rtol = 1e-10

All 20 calculations completed successfully.

The `rtol = 1e-10` convergence row passed all relevant A-3 acceptance
checks against the `rtol = 1e-11` reference.

Its principal results were:

    p_max_Pa = 3168091.98049835
    burn_time_s = 11.23604237447096
    F_max_N = 26.901139823064007
    impulse_total_Ns = 131.5111248533997
    total_mass_consistency_error = 9.432903467396914e-8
    propellant_mass_balance_error = 1.472169053329253e-10

### 24.2 Acceptance violations in the complete sweep

The complete candidate sweep reported four failed acceptance rows with
eight violations.

Those rows were:

- `rtol = 1e-6`;
- `cea_pressure_points = 21`;
- `cea_pressure_points = 41`;
- `cea_pressure_points = 81`.

The violations were confined to:

- maximum pressure and maximum thrust at `rtol = 1e-6`;
- burn time and total impulse for the coarser CEA tables.

No violation belonged to:

    rtol = 1e-10

The candidate's own acceptance result was therefore not invalidated by
the intentionally broad convergence and sensitivity matrix.

### 24.3 Relationship to preserved A-3 evidence

The A-3 relative-tolerance sweep already contains:

    1e-6
    1e-7
    1e-8
    1e-9
    1e-10
    1e-11

Its reference remains:

    rtol = 1e-11

Therefore `rtol = 1e-10` is already an evaluated A-3 convergence point,
not a newly invented A-4-only solver condition.

The preserved A-3 evidence was not modified by the A-4 candidate sweep.

### 24.4 Candidate disposition

A scalar product default of:

    rtol = 1e-10

is approved for final implementation review because:

1. LSODA, BDF, and Radau all accept scalar `rtol`;
2. the active-igniter ledger passes at this tolerance;
3. the `rtol = 1e-10` A-3 row passes;
4. principal-output changes against `rtol = 1e-11` are below the existing
   A-3 relative limits;
5. the complete candidate sweep executed successfully;
6. preserved A-3 evidence remains unchanged.

The following shall not change:

- the propellant-ledger acceptance limit of `1e-9`;
- the total-system ledger limit of `1e-6`;
- the A-3 sweep values;
- the A-3 reference value of `1e-11`;
- existing A-3 evidence files.

Before editing defaults, every `1e-9` occurrence shall be classified so
that solver defaults are not confused with ledger tolerances or preserved
evidence.
