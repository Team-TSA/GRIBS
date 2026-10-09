# Phase A-5.1-1: NASA CEA mass-basis findings

## Status

Technical finding: PASS

Gate G1 remains pending until this record and its reproducible probe have
been reviewed and committed.

No GRIBS product code was modified during this phase.

## Environment

- NASA CEA Python package: 3.3.4
- CEA library version: 3.3.4
- Universal gas constant exposed by CEA:
  `Ru = 8314.51 J/(kmol K)`
- Reproducible reference problem:
  NASA CEA RP-1311 Example 14
- Reference constraint: TP equilibrium
- Reference reactants: `H2(L)` and `O2(L)`
- Reference pressure: `0.05 atm = 0.0506625 bar`
- The reference problem contains gas and condensed product species.

## Official API semantics observed in CEA 3.3.4

Runtime descriptors from the installed official package state:

- `EqSolution.M`:
  "Molecular weight, gasses only (1/n)."
- `EqSolution.MW`:
  "Molecular weight, including condensed species."
- `EqSolution.n`:
  "Total moles."
- `EqSolution.nj`:
  "Species concentrations (kg-mol of species j per kg of mixture)."
- `EqSolution.ln_nj`:
  gas-species log concentrations, preserving amounts below the truncation
  applied to `nj`.
- `EqSolver.num_gas`:
  number of gas-phase species.
- `EqSolver.num_condensed`:
  number of condensed-phase species.
- `Mixture.species_names`:
  ordered species names for the mixture.

The official package's `matlab.py` reconstructs total moles using:

    total_moles = 1.0 / soln.M

The package's RP-1311 examples print `M` as `M, (1/n)`.

## Mass-basis derivation

For each CEA gas species, `nj` has units:

    kg-mol of gas species j / kg of total mixture

Therefore the gas mole amount on the total-mixture mass basis is:

    n_gas = sum_j(nj), over gas species

CEA reports:

    M_CEA = 1 / n_gas

Let `Yg` be the total gas mass fraction:

    Yg = kg gas / kg total mixture

The pure-gas mass-basis mean molecular weight is:

    M_gas
      = (kg gas / kg total mixture)
        / (kg-mol gas / kg total mixture)
      = Yg / n_gas

Thus:

    M_CEA = M_gas / Yg

and:

    Ru / M_CEA = Yg * Ru / M_gas

Consequently, `Ru / EqSolution.M` is already the effective gas constant on
the total-mixture mass basis. It is not the pure-gas specific gas constant.

## Density identity

The semantic probe independently verified:

    p / (rho * T) = Ru / M_CEA

at all eight RP-1311 Example 14 temperature points.

The maximum relative error was:

    4.38772888529282e-16

This confirms that CEA density uses `Ru / M_CEA` directly on the
total-mixture mass basis.

## Reproducible semantic probe

The probe is:

    verification/mass_basis/cea_mass_basis_semantic_probe.py

Its captured output is:

    verification/mass_basis/cea_mass_basis_semantic_probe.txt

Its captured process exit status is:

    verification/mass_basis/cea_mass_basis_semantic_probe_exit.txt

The probe reproduces the official RP-1311 Example 14 TP-equilibrium case at
eight temperatures:

    1000.0, 500.0, 350.0, 305.0,
    304.3, 304.2, 304.0, and 300.0 K

Four of those points contain condensed water.

## Numerical findings

At the strongest-condensation point, 300 K:

    Yg       = 0.346614384757677
    Yc       = 0.653385615242323
    M_CEA    = 64.179488892793 kg/kmol
    M_gas    = 22.2455340566374 kg/kmol
    Ru/M_CEA = 129.550891467658 J/(kg K)
    p/rho/T  = 129.550891467658 J/(kg K)

The reconstructed identity was:

    M_gas / Yg = M_CEA

and:

    Yg * Ru / M_gas = Ru / M_CEA

to near machine precision.

## Interpretation of the previous GRIBS EOS

Before correction, GRIBS uses:

    R_old = Ru / M_CEA
    Psi_old = Yg * R_old * T

However:

    Ru / M_CEA = Yg * Ru / M_gas

Therefore:

    Psi_old
      = Yg * (Yg * Ru / M_gas) * T
      = Yg^2 * (Ru / M_gas) * T

The consistent total-mixture mass-basis form is:

    R_eff = Ru / M_CEA
    Psi_correct = R_eff * T

Thus:

    Psi_old / Psi_correct = Yg

The probe confirmed this identity at every condensed point.

At 300 K:

    Psi_old / Psi_correct
      = 0.346614384757676
      approximately equals Yg
      = 0.346614384757677

The previous chamber EOS therefore applies the gas mass fraction twice.

## Tolerance rationale

The probe uses two acceptance levels.

For all eight points:

    tolerance = 1.0e-8

This accommodates the documented distinction between:

- `n`, which retains the internal total gas-mole amount, and
- `nj`, which is truncated for trace species.

The largest all-point discrepancy was approximately `2.33e-9` and occurred
at the non-condensed 1000 K point.

For the four condensed points:

    tolerance = 1.0e-10

The mass-basis identities at those points passed with errors from
approximately `1.6e-15` to `4.3e-14`.

The tolerance adjustment does not conceal a condensed-phase discrepancy.
The decisive condensed-point identities pass far more tightly than the
acceptance limit.

## Required model consequence

The first candidate for the internally consistent chamber EOS is:

    M_CEA = M_gas / Yg
    R_eff = Ru / M_CEA
    Psi = R_eff * T_used
    p * V = m_total * Psi

No additional multiplication by `Yg` belongs in `Psi`.

For a homogeneous-equilibrium nozzle using total mixture mass flow, the
corresponding first candidate is:

    mdot_total =
        Cd * At * p0
        * sqrt(gamma_s / (R_eff * T0))
        * f_gamma

with phase bookkeeping:

    mdot_gas = Yg * mdot_total
    mdot_condensed = (1 - Yg) * mdot_total

There must be no additional division of the reference flow by `Yg`.

These nozzle consequences are model decisions to be recorded in A-5.1-2.
They are not product-code changes authorized by this A-5.1-1 finding.

## Scope and limitations

This probe establishes the meaning of `EqSolution.M` and the associated EOS
mass basis using an official CEA 3.3.4 condensed-equilibrium reference case.

It does not by itself validate:

- the complete GRIBS chamber model,
- the GRIBS nozzle and thrust model,
- particle slip or finite-rate phase behavior,
- the reported high-condensation AP/Al case,
- calibrated performance over the full intended operating range.

The original high-condensation AP/Al input and output remain unavailable.
That domain remains:

    out_of_scope_pending_complete_input_and_output

The RP-1311 Example 14 probe is a semantic discriminator, not a GRIBS product
golden case.

## G1 technical conclusion

The following identities are established for NASA CEA 3.3.4:

    EqSolution.M = M_gas / Yg
    Ru / EqSolution.M = Yg * Ru / M_gas
    p / (rho * T) = Ru / EqSolution.M

The previous GRIBS chamber EOS multiplies by `Yg` one extra time.

Technical result:

    PASS

G1 may be marked complete after the findings, probe source, captured output,
exit status, and hashes have been reviewed and committed.
