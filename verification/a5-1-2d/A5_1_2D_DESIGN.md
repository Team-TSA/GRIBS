# A-5.1-2d inspection-contract design

- Baseline phase: after A-5.1-2c and before product-code modification
- Product physics authorized to change: no
- Numerical trajectory authorized to change: no
- Allowed external changes: diagnostics, warnings, and exit status

## EOS-history consistency metrics

The gate metric shall use one time-invariant run-wide reference mass:

    m_ref = max(max_t |m_total_balance(t)|, tiny)

    e_invariant =
        max_t |m_total_eos(t) - m_total_balance(t)| / m_ref

The existing shrinking-denominator metric shall remain temporarily as a
deprecated diagnostic:

    e_legacy =
        max_t |m_total_eos(t) - m_total_balance(t)|
        / max(|m_total_eos(t)|, tiny)

The existing output field `total_mass_consistency_error` shall retain its
old meaning during A-5.1-2d. New explicitly named fields shall carry the
invariant metric, its reference mass, and its definition version.

## Structured post-run report

One authoritative function shall calculate both metrics and the final
total-mass conservation result. Printing, summary generation, warnings,
and CLI exit behavior shall consume the same report.

## Output-before-failure contract

After a simulation has produced a history, diagnostic CSV and summary
JSON shall be written before a conservation or ledger failure produces a
nonzero process exit.

Summary text shall also be written before the verdict exit. The figure is
generated after the machine-readable diagnostic outputs.

## Mass-ledger contract

The existing ledger equations and tolerances remain unchanged in
A-5.1-2d.

A false ledger verdict shall:

- be represented in warnings;
- identify total-system and propellant sub-verdicts;
- be reflected in the CLI exit status;
- not suppress diagnostic output.

Igniter commanded and injected mass remain explicitly separated from
propellant-derived generated mass.

## Self-test contract

The physical calculations of the existing self-tests remain unchanged.
The test-name set and result count shall become inspectable and protected.
A conditionally inapplicable test shall be represented explicitly rather
than silently disappearing.

## reference_constant contract

The existing single-phase behavior is protected without changing the
numerical path:

    gas_mass_fraction = 1.0
    condensed_mass_fraction = 0.0

No new configurable phase-fraction freedom is introduced.

## Numerical invariance contract

For the rtol=1e-10 wide baseline:

- all CSV headers, row counts, phases, regimes, booleans, and numerical
  cells shall remain identical;
- existing summary JSON fields shall remain recursively identical after
  excluding only approved path fields and newly added A-5.1-2d
  diagnostic or warning fields;
- chamber Psi, nozzle flow, velocity, thrust, cache semantics, and
  physical configuration shall not change.

Approved path-dependent exclusions:

- configuration_source;
- output_files.

New A-5.1-2d fields may be excluded only by explicit name in the
comparison script. Existing fields shall not be broadly excluded.

## Prohibited changes

A-5.1-2d shall not modify:

- Psi or dPsi/dp;
- chamber state equations;
- nozzle mass flow;
- exit velocity;
- thrust;
- cache schema or mass-basis semantics;
- CEA table construction;
- ordinary default numerical trajectories;
- mass-ledger tolerances.

A-5.1-2d shall not be combined with A-5.1-4.

## Baseline evidence storage normalization

The generated `rtol=1e-7` reproduction log originally contained one
extra blank line at EOF. Before finalizing the baseline evidence commit,
only that trailing blank line was removed to satisfy the repository
`git diff --check` contract. Log messages, ordering, numerical values,
and the recorded exit status were not changed.
