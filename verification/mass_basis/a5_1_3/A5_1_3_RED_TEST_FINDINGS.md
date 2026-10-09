# A-5.1-3 Independent Red-Test Findings

## Scope

A-5.1-3 defines corrected CEA mass-basis acceptance contracts before
the A-5.1-4 product correction. This unit changes no product code.

The tests are isolated under tests/red. Normal pytest discovery excludes
that directory. Explicit execution uses the following command:

    pytest -q tests/red/test_a5_1_3_mass_basis_red.py

## Validation result

- Base HEAD: e83825636a8db166ee0436b432309e8d86af2ca6
- Normal suite exit status: 0
- Normal suite result: 83 passed
- Explicit red-suite exit status: 1
- Explicit red-suite result: 6 failed, 5 passed, 1 xfailed
- Collection, import, fixture, setup, and unexpected-XPASS errors: none
- Normal collection includes red tests: no
- Product code modified by A-5.1-3: no

## Intended pre-fix failures

### Chamber Psi applies Yg twice

The product table path produced:

    actual Psi       = 698418.84 J/kg
    corrected Psi    = 997741.2 J/kg
    actual/expected  = 0.7
    Yg               = 0.7

The actual value equals the rejected legacy expression:

    Psi_legacy = Yg * R_eff * T_used

The acceptance expression is:

    Psi_corrected = R_eff * T_used

### Nozzle flow applies an additional division by Yg

The product nozzle path produced:

    actual mdot_total     = 0.04498477993188038 kg/s
    corrected mdot_total  = 0.03148934595231626 kg/s
    actual/expected       = 1.428571428571429
    1/Yg                  = 1.428571428571429

The nozzle base equation must directly return total-mixture mass flow.
No additional division by Yg is permitted.

### Corrected and legacy total-flow relation is violated

The measured pre-fix relation is:

    actual corrected/legacy ratio = 1.707469441906277
    expected ratio                = 1.195228609334394
    expected expression           = 1/sqrt(Yg)

The extra factor is consistent with the additional nozzle division by Yg.

### Cache schema 3 is still accepted

A complete schema-3 cache was loaded without raising
ThermochemistryCacheError. The corrected mass-basis implementation
must reject the legacy schema identity.

### Metadata retains obsolete mass-basis claims

The following obsolete phrases remain in the product source:

    EqSolution.M (gas-phase, kg/kmol)
    gas_phase_M (gas-phase M;

Corrected metadata must describe EqSolution.M on the total
equilibrium-mixture mass basis, not as pure-gas-basis molecular weight.

## Passing protection contracts

The following acceptance contracts already pass and must remain protected:

- M_CEA = M_gas / Yg
- R_eff = Yg * R_gas
- p / (rho * T) = R_eff
- mdot_gas = Yg * mdot_total
- mdot_condensed = Yc * mdot_total
- mdot_gas + mdot_condensed = mdot_total
- Momentum thrust uses mdot_total once
- The Yg = 1 single-phase limit is unchanged

## Synthetic-fixture boundary

The synthetic state uses Yg = 0.7 and Yc = 0.3 to make an extra factor
of Yg structurally unambiguous. This is a unit-level equation probe,
not a substitute for measured CEA evidence. The existing
medium-condensation evidence remains the empirical record.

## External comparison gap

The public v0.4.3-alpha default gas_phase_M comparison artifact has not
been obtained. Its reserved acceptance test is therefore an explicit
XFAIL. This XFAIL is not evidence for or against the corrected equations.

## A-5.1-4 handoff

A-5.1-4 must make the six intended failures pass without regressing
the five currently passing protection tests. Cache migration, metadata
migration, and physical Psi and nozzle corrections belong to A-5.1-4,
not this commit.

## Evidence-log normalization

The captured red pytest output originally contained five lines whose
content consisted of pytest failure-prefix text followed by trailing
spaces. Before commit, trailing horizontal whitespace was removed from
every line of red-pytest.txt so that git diff --check succeeds.

The normalization changed no test names, failure text, measured values,
result classification, line order, or line count. It removed trailing
spaces only.

Pre-normalization SHA-256: 3ce57c9ee1bd5c9b19fbe799a26671516709bcbaf4b7ad657f3ae1ac20a60408
Post-normalization SHA-256: 3f5a88a526ef37dc8c2363f158248011b31a055bc2791b65fc9f08fb44ffabc9
Line count before: 78
Line count after: 78
Byte count before: 5378
Byte count after: 5353
