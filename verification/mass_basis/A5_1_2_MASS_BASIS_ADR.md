# ADR: GRIBS CEA mass basis, chamber EOS, nozzle flow, and thrust

- Status: Accepted for implementation planning
- Phase: A-5.1-2
- Decision gate: G2
- Date: 2026-10-09
- Baseline version: GRIBS v0.5.0-alpha
- Baseline commit: 95e25731586229627e8ccd8f39a6215c7a69d60c
- Scope: model decision only
- Product-code modification authorized by this document: no

## 1. Context

GRIBS uses NASA CEA equilibrium properties to construct the chamber
equation of state and the nozzle model.

The previous implementation interpreted `EqSolution.M` as the pure-gas
mass-basis mean molecular weight and used:

    Rg = Ru / M_CEA
    Psi_old = Yg * Rg * T_used

It also calculated an isentropic nozzle reference flow with `Rg` and then
used:

    mdot_total_old = mdot_gas_reference / Yg

Phase flows were subsequently reconstructed as:

    mdot_gas = Yg * mdot_total_old
    mdot_condensed = (1 - Yg) * mdot_total_old

Phase A-5.1-1 established from the NASA CEA 3.3.4 runtime API and a
reproducible condensed-equilibrium semantic probe that:

    EqSolution.nj:
        kg-mol of species j per kg of total mixture

    n_gas = sum(nj over gas species)

    EqSolution.M = 1 / n_gas

If `Yg` is the gas mass fraction and `M_gas` is the pure-gas mass-basis
mean molecular weight, then:

    M_gas = Yg / n_gas

and therefore:

    M_CEA = M_gas / Yg

    Ru / M_CEA = Yg * Ru / M_gas

The semantic probe also established:

    p / (rho * T) = Ru / M_CEA

Consequently, `Ru / EqSolution.M` is already an effective gas constant on
the total-mixture mass basis.

## 2. Decision

### 2.1 CEA molecular-weight meaning

The value returned by `EqSolution.M` shall be treated as:

    M_CEA = 1 / n_gas

where `n_gas` is expressed per unit total-mixture mass.

It shall not be described or stored as the pure-gas mass-basis molecular
weight.

Preferred semantic name:

    cea_mixture_basis_molecular_weight_kg_kmol

The pure-gas mean molecular weight, if exposed as a diagnostic, shall be
derived explicitly as:

    M_gas = Yg * M_CEA

and named:

    gas_phase_mean_molecular_weight_kg_kmol

`EqSolution.MW` remains a separate diagnostic representing molecular
weight including condensed species. It shall not enter the chamber EOS or
the current nozzle model.

### 2.2 Effective gas constant

Define:

    R_eff = Ru / M_CEA

Its units are:

    J / (kg total mixture K)

Equivalently:

    R_eff = Yg * R_gas

where:

    R_gas = Ru / M_gas

has units:

    J / (kg gas K)

The solver-facing value currently named `gas_constant_J_kgK` shall be
renamed or explicitly versioned so that its total-mixture mass basis is
unambiguous.

Preferred name:

    effective_gas_constant_J_kg_total_K

If `R_gas` is retained as a diagnostic, it must be stored separately and
must not be confused with `R_eff`.

### 2.3 Chamber equation of state

For the homogeneous-equilibrium model with negligible condensed volume,
the chamber state relation shall be:

    Psi = R_eff * T_used

    p * Vg = m_total * Psi

where:

    T_used = eta_T0 * T_CEA

and `m_total` includes gas and condensed product mass.

No additional multiplication by `Yg` is permitted in `Psi`.

The previous expression:

    Psi_old = Yg * (Ru / M_CEA) * T_used

is rejected because it applies the gas mass fraction twice.

### 2.4 Initial chamber mass

The initial total chamber inventory shall remain:

    m_total0 = p0_init * Vg0 / Psi(p0_init)

with the corrected:

    Psi = R_eff * T_used

A diagnostic gas inventory may be reported as:

    m_gas_equilibrium = Yg * m_total

but it shall not replace the total mass state.

The former diagnostic expression based on a separately interpreted
`Theta = Rg*T` shall not be represented as a physically independent
chamber EOS after the migration.

### 2.5 Pressure derivative

The pressure ODE structure remains:

    dp/dt =
        [Psi * (mdot_gen_total - mdot_out_total) - p * Ab * r]
        / [Vg * (1 - p * Psi_prime / Psi)]

The definition changes to:

    Psi(p) = R_eff(p) * T_used(p)

The pressure-table backend shall construct `Psi` at each table node and
shall create a dedicated PCHIP interpolator in `ln(p)`.

The derivative shall remain analytic with respect to that interpolator:

    dPsi/dp =
        d[PCHIP(Psi)]/d[ln(p)] / p

It shall not be replaced by a finite-difference derivative or by an
independent product of separately interpolated derivatives.

The physicality condition remains:

    1 - p * Psi_prime / Psi > 0

over the configured pressure range.

### 2.6 Nozzle total mass flow

For the current homogeneous-equilibrium, complete-entrainment
approximation, the isentropic choking formula shall directly produce total
mixture mass flow:

    mdot_total =
        Cd * At * p0
        * sqrt(gamma_s / (R_eff * T0))
        * f_gamma

where:

    f_gamma =
        [2 / (gamma_s + 1)]
        ^ [(gamma_s + 1) / (2 * (gamma_s - 1))]

for choked flow.

The unchoked form shall use the same total-mixture mass basis.

There shall be no additional division by `Yg`.

The previous operation:

    mdot_total = mdot_gas_reference / Yg

is rejected.

The internal and output name `mdot_gas_reference` is semantically wrong
once the corrected `R_eff` is used. It shall be migrated to an explicit
total-flow name.

Preferred name:

    mdot_total_reference

or, where no separate adjustment follows:

    mdot_total

### 2.7 Phase-flow decomposition

After `mdot_total` has been calculated, the equilibrium phase bookkeeping
shall be:

    mdot_gas = Yg * mdot_total

    mdot_condensed = (1 - Yg) * mdot_total

and:

    mdot_total = mdot_gas + mdot_condensed

This split is diagnostic and supports the mass ledger. It does not alter
the total mass flow a second time.

### 2.8 Exit state and velocity

The existing isentropic area-Mach and pressure-ratio structure may remain
for the minimum consistent model.

The solver-facing thermodynamic coefficient used by that model shall be
`R_eff`, not a quantity described as the pure-gas specific gas constant.

For the accepted minimum model:

    Te = T0 / [1 + (gamma_s - 1) * Me^2 / 2]

    ve = Me * sqrt(gamma_s * R_eff * Te)

This is a homogeneous-equilibrium effective-mixture approximation. It is
not a claim that gas and particles independently obey a single-phase
perfect-gas velocity relation.

### 2.9 Thrust

Under the retained complete-velocity-equilibrium approximation:

    F =
        eta_F * mdot_total * ve
        + (pe - pa) * Ae_effective

The momentum term uses total mass flow exactly once.

The minimum model assumes:

    v_gas = v_condensed = ve

No particle slip, particle-size distribution, finite-rate phase change,
or separate phase momentum equation is introduced in A-5.1.

This limitation must remain visible in output metadata and model
limitations.

### 2.10 Empirical parameters

The following parameters shall not be recalibrated to reproduce
pre-correction results:

- discharge coefficient `Cd`
- thrust efficiency `eta_F` / `eta_thrust`
- temperature efficiency `eta_T0`
- Saint-Robert burn-law coefficients `a` and `n`
- erosive-burning coefficients
- throat-erosion coefficients

Old results are not calibration targets because they were produced with
an inconsistent mass basis.

### 2.11 Atomic implementation boundary

The following corrections form one atomic physical change:

1. CEA molecular-weight meaning and naming
2. effective gas constant meaning and naming
3. chamber `Psi`
4. analytic `dPsi/dp`
5. initial chamber inventory
6. nozzle total mass flow
7. phase-flow decomposition
8. exit velocity basis
9. momentum thrust basis
10. mass-ledger metadata
11. cache schema and rejection
12. output and provenance definitions
13. self-tests and independent regression tests

An intermediate state in which only the chamber EOS or only the nozzle is
corrected shall not be published or treated as usable.

Phase A-6 modularization shall not be included in the same change set.

#### A-5.1-4 / A-5.1-5 integration rule

A-5.1-4 and A-5.1-5 are separate verification checkpoints, but they are
not separate releasable physical states.

A-5.1-4 may implement and test the chamber EOS, nozzle-flow, and thrust
correction before all configuration and output names have been migrated.
However, its intermediate state shall not:

- be published or tagged as a usable GRIBS version;
- be merged independently as the completed mass-basis correction;
- silently read cache schema 3;
- produce new golden data;
- be used for G4 applicability claims;
- be treated as a calibration target.

If the complete cache-schema migration is not finished within A-5.1-4,
the intermediate implementation shall either reject all schema-3 cache
files or run with cache reuse disabled. It shall not execute corrected
physics using a schema-3 cache.

A-5.1-4 and A-5.1-5 shall be reviewed as parts of the same physical
correction and shall remain on the same audit branch or pull request until
both checkpoints are complete.

The atomic correction is considered complete only after:

1. EOS, nozzle flow, exit velocity, and thrust use the accepted mass basis;
2. old cache semantics are rejected;
3. configuration and output semantics are migrated;
4. all independent red tests have turned green; and
5. the full regression suite passes.

### 2.12 Cache migration

The current thermochemistry cache schema is version 3.

Schema 3 embeds or permits the rejected meanings:

- `gas_phase_molecular_weight_kg_kmol` for `EqSolution.M`
- `gas_constant_J_kgK` without a declared mass basis
- `psi_J_kg = Yg * (Ru/M_CEA) * T_used`
- provenance describing `EqSolution.M` as pure gas-phase molecular weight
- density provenance describing CEA density as based on pure-gas `Rg`

Schema 3 shall be rejected after the implementation change.

The new cache schema shall:

- have a new integer schema version;
- change the cache-key semantic version or equivalent key material;
- require regeneration rather than silent reuse;
- store `M_CEA` with an unambiguous total-mixture-basis name;
- store `R_eff` with an unambiguous total-mixture-basis name;
- store corrected `Psi = R_eff*T_used`;
- retain `Yg` and `Yc`;
- optionally store derived `M_gas` and `R_gas` as clearly marked
  diagnostics;
- record the mass basis in provenance;
- reject schema 3 with an actionable rebuild message.

### 2.13 Configuration migration

The existing configuration token:

    thermochemistry.cea_python.molecular_weight = "gas_phase_M"

is misleading because `EqSolution.M` is not pure-gas mass based when
condensed species are present.

A-5.1-5 shall either:

1. replace the token with a new explicit semantic token and reject the old
   token, or
2. remove the unnecessary choice if only `EqSolution.M` is valid for the
   accepted EOS.

`total_MW` shall remain unavailable as an EOS selection.

Any compatibility handling must be explicit and must not preserve the old
physical interpretation silently.

### 2.14 Output migration

Output and provenance fields shall be revised so their mass basis is
unambiguous.

At minimum, the migration shall cover:

- molecular-weight field names and definitions
- gas-constant field names and definitions
- `Psi` definition
- density note
- nozzle reference-flow name
- nozzle entrainment description
- mass-ledger `psi_definition`
- interpolation quantity descriptions
- self-test labels
- summary and text-report descriptions
- thermo dump schema and provenance

The following mass-ledger quantities remain total-product-mass based:

- `m_total_eos`
- `m_total_balance`
- `mdot_out_total`
- integrated discharged mass
- initial and final chamber inventory

The phase quantities remain derived decompositions:

- `m_gas_equilibrium = Yg * m_total_eos`
- `m_condensed = Yc * m_total_eos`
- `mdot_out_gas = Yg * mdot_out_total`
- `mdot_out_condensed = Yc * mdot_out_total`

The backward-compatible alias `mdot` may remain only if it is explicitly
documented as total mass flow.

### 2.15 Tests required before implementation

A-5.1-3 shall add independent tests that fail against the current
implementation.

At minimum, the red tests shall establish:

1. For `Yg < 1`, `Psi = R_eff*T`, not `Yg*R_eff*T`.
2. `M_CEA = M_gas/Yg`.
3. `R_eff = Yg*R_gas`.
4. `p/(rho*T) = R_eff`.
5. The nozzle formula directly returns `mdot_total`.
6. There is no additional `/Yg`.
7. `mdot_gas + mdot_condensed = mdot_total`.
8. `mdot_gas = Yg*mdot_total`.
9. Momentum thrust uses `mdot_total` once.
10. The single-phase limit `Yg = 1` remains unchanged.
11. Old cache schema 3 is rejected.
12. New metadata does not claim `EqSolution.M` is pure-gas based.

The tests shall use at least one clearly discriminating condensed state.
A low-condensation case with `Yg` approximately equal to 1 is not
sufficient by itself.

### 2.16 Validation regions

The corrected model shall be assessed in three regions:

1. single phase or effectively single phase;
2. intermediate condensation;
3. high condensation.

The official RP-1311 Example 14 case is accepted as a semantic probe, not
as a GRIBS product golden case.

The available GRIBS pre-fix dump has:

    Yg_min approximately 0.997962

and therefore belongs only to the low-condensation region.

The original reported high-condensation AP/Al input and output have not
been obtained. Until they are available and reviewed, that region remains:

    out_of_scope_pending_complete_input_and_output

G4 and G5 shall not declare the high-condensation AP/Al domain validated
without the complete source input and output.

## 3. Consequences

### Positive consequences

- Chamber pressure, stored mass, nozzle flow, and thrust use one explicit
  total-mixture mass basis.
- The extra chamber `Yg` factor is removed.
- The extra nozzle `/Yg` conversion is removed.
- The single-phase limit remains algebraically identical.
- Cache and output meanings become auditable.
- Phase-flow diagnostics remain available without altering total flow.

### Negative consequences and limitations

- Numerical results will change whenever `Yg < 1`.
- Old caches and old golden data cannot be reused.
- Existing comparisons based on the pre-correction model become historical
  evidence rather than acceptance targets.
- The retained nozzle model is still an approximation for particle-laden
  flow.
- High-condensation thrust accuracy remains unvalidated.
- Output consumers may require migration to new field names or schemas.

## 4. Rejected alternatives

### Keep `Psi = Yg*(Ru/M_CEA)*T`

Rejected because CEA already embeds the total-mixture mass basis in
`M_CEA = 1/n_gas`.

### Treat `EqSolution.M` as `M_gas`

Rejected by the CEA API semantics and the condensed-equilibrium probe.

### Use `EqSolution.MW` in the EOS

Rejected because `MW` includes condensed-species mole accounting and does
not satisfy the gas-pressure EOS identity established by the probe.

### Retain `mdot_total = mdot_reference/Yg`

Rejected because the reference formula already uses `R_eff = Ru/M_CEA`
and therefore already yields total-mixture mass flow under the accepted
homogeneous-equilibrium model.

### Correct only the chamber EOS

Rejected because it would leave the nozzle and thrust on an inconsistent
mass basis.

### Correct only the nozzle

Rejected because chamber mass, pressure evolution, and initial inventory
would remain inconsistent.

### Retune empirical coefficients

Rejected because it would conceal a physical-definition correction inside
calibration.

### Reuse schema-3 cache files

Rejected because their names, `Psi` column, and provenance encode the old
meaning.

### Combine with Phase A-6 modularization

Rejected because refactoring would enlarge the review surface and obscure
the atomic physical correction.

## 5. Implementation sequence after G2

After this ADR is reviewed and accepted:

1. A-5.1-3 adds independent red tests.
2. The red state is recorded.
3. A-5.1-4 applies the atomic EOS/nozzle/thrust correction.
4. A-5.1-5 performs cache, configuration, and output migration.
5. A-5.1-6 reruns the three validation regions.
6. A-5.1-7 performs independent review and decides whether A-6 and A-8 may
   resume.

No product-code change is authorized before the A-5.1-3 red-test state is
fixed and reviewed.

## 6. Decision summary

Accepted model:

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

No extra multiplication of `Psi` by `Yg` is permitted.

No extra division of nozzle flow by `Yg` is permitted.
