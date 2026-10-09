# A-5.1-2a: Independent-review addendum to the mass-basis ADR

- Status: Accepted as a planning and verification addendum
- Phase: A-5.1-2a
- Date: 2026-10-09
- Baseline version: GRIBS v0.5.0-alpha
- Baseline commit: 4449132623804a017113bf5f25a37727bc02d11d
- Parent ADR: `A5_1_2_MASS_BASIS_ADR.md`
- Product-code modification authorized by this addendum: no

## 1. Purpose and precedence

This addendum records the changes required after independent review of the
G1 and G2 mass-basis decisions.

The parent ADR remains authoritative for the accepted physical model:

    M_CEA = M_gas / Yg
    R_eff = Ru / M_CEA
    Psi = R_eff * T_used
    p * Vg = m_total * Psi

    mdot_total =
        Cd * At * p0
        * sqrt(gamma_s / (R_eff * T0))
        * f_gamma

    mdot_gas = Yg * mdot_total
    mdot_condensed = (1 - Yg) * mdot_total

    F =
        eta_F * mdot_total * ve
        + (pe - pa) * Ae_effective

No extra multiplication of `Psi` by `Yg` and no extra division of nozzle
flow by `Yg` are permitted.

Where this addendum differs from section 5 of the parent ADR, this
addendum governs the implementation and verification sequence.

## 2. G1 evidence grading

The G1 conclusion is accepted, but its supporting observations do not all
have the same evidentiary role.

Primary semantic evidence:

1. `EqSolution.nj` is reported in kg-mol of species j per kg of total
   mixture.
2. `EqSolution.M * EqSolution.n = 1`.
3. `EqSolution.n` agrees with the sum of gas-species `nj`.
4. `p/(rho*T) = Ru/M_CEA`.

The observed relation:

    Psi_old / Psi_correct = Yg

is an algebraic consequence of substituting the old and accepted
definitions. It is useful for estimating the effect of the old formula,
but it is not independent evidence for the CEA semantics.

RP-1311 Example 14 remains a semantic probe. It is not a GRIBS product
golden case.

## 3. Meaning of `EqSolution.n`

For the observed CEA 3.3.4 runtime API:

    EqSolution.n = sum(nj over gas species)

with the amount expressed per kg of total equilibrium mixture.

Therefore:

    M_CEA = 1 / EqSolution.n

is an effective molecular-weight quantity on the total-mixture mass basis
for the gas-pressure EOS. It shall not be described as the pure-gas
mass-basis molecular weight.

The pure-gas mean molecular weight is derived as:

    M_gas = Yg / EqSolution.n
          = Yg * M_CEA

`EqSolution.MW` is a separate diagnostic that includes condensed-species
mole accounting. It shall not be substituted into the chamber EOS or the
current nozzle model.

## 4. Numerical tolerances

The semantic probe may use tight tolerances appropriate to direct
runtime-API identities and deterministic table values.

Product-path tests shall use tolerances justified for the complete path,
including:

- CEA library behavior;
- tabulation;
- PCHIP interpolation;
- ODE integration;
- event localization;
- numerical quadrature;
- floating-point environment differences.

A probe tolerance shall not automatically become a production acceptance
tolerance.

Each product-path test shall state:

1. the quantity being compared;
2. whether the comparison is absolute or relative;
3. the numerical operations included;
4. the tested environment;
5. why the tolerance distinguishes a structural mass-basis error from
   ordinary numerical error.

The pre-fix reproduction demonstrated that the current EOS-history metric
is sensitive to solver tolerance:

    rtol = 1e-7:
        final total-mass relative error = 2.94e-09, PASS
        EOS-history maximum relative error = 3.05e-05, FAIL
        process exit = 1
        output files = 0

    rtol = 1e-10:
        final total-mass relative error = 3.48e-11, PASS
        EOS-history maximum relative error = 3.60e-07, PASS
        process exit = 0
        output files = 4

This result supports replacement of the shrinking instantaneous-mass
normalization with a time-invariant reference scale in A-5.1-2d. The old
metric may remain temporarily only as a deprecated diagnostic.

## 5. Relationship between homogeneous and legacy total results

After the accepted correction, `homogeneous_equilibrium` may coincide in
total quantities with `single_phase_legacy` when both effectively use the
same `M_CEA`, `R_eff`, total-flow formula, and velocity basis.

Any such equality is not assumed from algebra alone. It shall be fixed by
a protection test against the actual product API.

The expected distinction after correction is primarily:

- equilibrium gas/condensed phase decomposition;
- phase-fraction diagnostics;
- model metadata;
- applicability and validation status.

The claim that corrected totals also coincide with the public
v0.4.3-alpha default `gas_phase_M` path remains unconfirmed at the current
pre-fix HEAD. It requires a reproducible external-version comparison.

## 6. Model capability and explicit non-capabilities

The accepted homogeneous-equilibrium model is a minimum
complete-entrainment approximation. It is not a complete two-phase-flow
model.

It assumes:

- instantaneous equilibrium phase fractions;
- complete velocity equilibrium;
- negligible condensed-phase chamber volume;
- one effective-mixture nozzle velocity;
- total product mass as the conserved chamber state.

It does not implement:

- particle slip;
- finite-rate condensation or evaporation;
- condensed-phase chamber volume;
- separate gas and condensed momentum equations;
- phase-density lag;
- particle-size distribution or particle-size effects.

`m_gas_equilibrium = Yg*m_total` is a derived instantaneous equilibrium
diagnostic. It is not an independently integrated conserved state.

These limitations shall remain visible in output metadata and validation
scope.

## 7. `eta_T0` constraint

The current implementation applies:

    T_used = eta_T0 * T_CEA

without recomputing `M_CEA` or `gamma_s` at the modified temperature.

Therefore, when `eta_T0 != 1`, a comparison of GRIBS density or sound
speed with the original CEA state is not a same-state validation.

The mass-basis correction shall not expand into a fundamental redesign of
`eta_T0`.

Direct CEA comparison in G4a is restricted to:

    eta_T0 = 1

Any other value requires separate interpretation and validation.

## 8. Blowdown-temperature constraint

The current model continues to use the combustion HP equilibrium property
table after propellant burnout and during blowdown.

This is a known energy-model limitation. The current model does not
integrate an independent energy state that would determine a physically
separate blowdown temperature history.

The limitation shall be documented in output and validation records. Its
fundamental correction is assigned to Phase B and shall not be mixed into
the A-5.1 mass-basis change.

## 9. Legacy-result quarantine

Pre-correction results shall not be discarded, but they shall not be used
as acceptance targets for the corrected physics.

Legacy classification shall be scoped by machine-readable provenance,
including where available:

- GRIBS version and commit;
- thermochemistry backend;
- two-phase mode;
- cache schema;
- mass-basis semantics version;
- `Psi` definition;
- molecular-weight definition;
- effective deviation of `Yg` from 1;
- configuration migration history.

Low-condensation historical results may remain useful as regression
context because the numerical difference can be small. That does not make
their old semantics valid.

Intermediate- and high-condensation pre-fix results shall be clearly
quarantined from corrected golden data and calibration targets.

Empirical parameters shall not be retuned to reproduce quarantined
pre-correction results.

## 10. G4a and G4b

Validation is divided into two gates.

### G4a: reproducible validation domain

G4a shall cover:

1. single-phase or effectively single-phase behavior;
2. intermediate condensation;
3. a public or otherwise reproducible alternative high-condensation case.

Direct comparison to CEA density, sound speed, or equivalent same-state
properties is restricted to `eta_T0 = 1`.

The property-table pressure range shall cover the corrected MEOP with
adequate margin, and:

    property_extrapolation_fraction = 0

is required.

The currently recorded pre-fix wide case uses a 20 MPa table and reaches
a pre-fix MEOP of approximately 10.548 MPa without extrapolation. This
does not prove that 20 MPa is sufficient after correction. The corrected
MEOP shall determine the required G4a table ceiling.

### G4b: reported AP/Al case

The originally reported high-condensation AP/Al case remains:

    out_of_scope_pending_complete_input_and_output

until the complete source input and output are obtained and reviewed.

G4b shall not be represented as validated from a surrogate configuration.

G4b remaining open does not automatically block all later work if G4a
defines and machine-readably enforces a narrower validated domain.
A-5.1-7 shall make that decision.

## 11. Environment differences and G0 evidence gap

The frozen G0 thermo dump was a low-condensation case:

    Yg_min approximately 0.997962

It could not strongly distinguish the old and accepted mass bases.

A-5.1-2c therefore added a pre-fix intermediate-condensation evidence
case with:

    points = 81
    Yg_min = 0.7238397417377189
    Yg_max = 0.7265536177401214
    Yc_min = 0.2734463822598786
    Yc_max = 0.2761602582622811
    max |Yg + Yc - 1| = 0

The evidence was obtained on:

- Python 3.14.4;
- NumPy 2.5.3;
- SciPy 1.18.1;
- Matplotlib 3.11.2;
- pytest 9.1.1;
- official CEA Python package and library 3.3.4;
- `cea.R = 8314.51`.

Independent-review numerical values may differ because of runtime,
library, interpolation, integration, and machine-environment differences.
Local reproducible measurements govern local classification. Differences
shall be recorded rather than adjusted away.

## 12. Theta and Psi

The corrected homogeneous-equilibrium solver has one authoritative
pressure-state quantity:

    Psi = R_eff * T_used

`Theta` shall not remain as a second, physically independent EOS for the
same homogeneous-equilibrium state.

If `Theta` is retained for any of the following reasons:

- historical output comparison;
- single-phase-limit testing;
- legacy backend support;
- diagnostic decomposition;

then it shall be explicitly labeled as diagnostic or legacy-scoped and
shall not independently drive the corrected pressure ODE.

The corrected solver-facing interpolation and derivative shall be built
from one authoritative `Psi` table.

Tests shall prevent divergence between:

- stored `psi_J_kg`;
- solver-facing Psi;
- `R_eff*T_used`;
- reported `psi_definition`.

## 13. Revised implementation order

The original parent-ADR sequence beginning directly with A-5.1-3 is
superseded.

The required order is:

1. Process 0: prerequisite verification and evidence-gap closure.
2. A-5.1-2a: this ADR and guidance addendum.
3. A-5.1-2b: reproduction and classification of review claims.
4. A-5.1-2c: condensed pre-fix evidence extension.
5. A-5.1-2d: inspection-contract correction without numerical trajectory
   changes.
6. A-5.1-3: independent red tests.
7. A-5.1-4: atomic EOS, nozzle-flow, and thrust correction.
8. A-5.1-5: cache, configuration, output, and provenance migration.
9. A-5.1-6: regional reruns and applicability decision.
10. A-5.1-7: independent review and G5 decision.

Process 0, A-5.1-2a, A-5.1-2b, A-5.1-2c, and A-5.1-2d shall not be
skipped before A-5.1-3.

## 14. A-5.1-2d constraints

A-5.1-2d is the only product-code change authorized before A-5.1-3.

It is restricted to inspection and failure-reporting contracts, such as:

- a time-invariant EOS consistency normalization;
- temporary retention of the old metric as deprecated;
- visibility of the mass-ledger verdict;
- explicit warnings and exit status;
- writing diagnostic output before conservation-gate termination;
- distinguishing igniter mass from conservation error;
- guaranteeing the self-test name set and count;
- a default gas mass fraction of 1.0 for `reference_constant`.

A-5.1-2d shall not alter:

- chamber `Psi`;
- nozzle flow;
- exit velocity;
- thrust;
- cache mass-basis semantics;
- ordinary default numerical trajectories.

For the default accepted comparison case, numerical values in
`time_history.csv` and `summary.json` shall be proven bit-for-bit
unchanged except for newly added diagnostic or contract fields whose
addition is explicitly reviewed.

Warning text and process exit status may change.

A-5.1-2d and A-5.1-4 shall not be placed in the same commit.

## 15. Gate status after this addendum

- G0: achieved, with the original low-condensation limitation now
  supplemented by A-5.1-2c intermediate-condensation evidence.
- G1: achieved, subject to the evidence grading and tolerance language in
  this addendum.
- G2: achieved, subject to the revised scope and implementation order in
  this addendum.
- G3: not achieved.
- G4a: not achieved.
- G4b: not achieved.
- G5: not achieved.

A-6 modularization remains suspended until the applicable A-5.1 gates
permit reconsideration.

A-8 product-golden freezing remains suspended until G5.
