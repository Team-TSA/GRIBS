# Independent-review claim classification

Classification vocabulary: `reproduced`, `not_reproduced`, `environment_dependent`.

Base commit: `4449132623804a017113bf5f25a37727bc02d11d`.

| # | Claim | Classification | Local evidence and disposition |
|---:|---|---|---|
| 1 | Post-run EOS metric uses a small instantaneous-mass denominator and can fail | reproduced | rtol=1e-7 gave EOS-history error 3.05e-05, exit 1, and zero output files, while final mass error was 2.94e-09 and passed. |
| 2 | Ledger thresholds are not coupled to solver rtol and verdict can change | environment_dependent | Fixed thresholds are confirmed in source. Local rtol=1e-10 ledger passed. The rtol=1e-7 run was stopped before summary/ledger output, so the claimed full reversal matrix was not locally reproduced. |
| 3 | ledger.passed=False is not surfaced through warning or CLI exit | reproduced | Static product-path inspection confirms main does not consume ledger.passed. A ledger-failing local case has not yet been generated, so only the implementation defect is classified as reproduced. |
| 4 | Cache reader does not explicitly compare stored schema | reproduced | Static inspection confirms no explicit stored-schema comparison. Runtime rejection behavior remains assigned to A-5.1-3/A-5.1-5. |
| 5 | Stored psi_J_kg is not authoritative and current consistency checking is self-referential | reproduced | Static inspection confirms solver-facing Psi is reconstructed and the stored column is not the interpolation authority. |
| 6 | After the physical fix, homogeneous and legacy total outputs coincide | not_reproduced | Requires the post-fix implementation, which is intentionally absent at current HEAD. Preserve as an A-5.1-4 protection-test requirement. |
| 7 | Public v0.4.3-alpha gas_phase_M totals coincide with the corrected implementation | not_reproduced | Requires both the post-fix implementation and the external v0.4.3 comparison artifact. It is not testable at current pre-fix HEAD. |
| 8 | Legacy quarantine must be limited by version, backend, mode, and effective Yg deviation | reproduced | The logical need for scoped quarantine follows from the low-condensation evidence and the review analysis. Implementation is deferred to A-5.1-5. |
| 9 | eta_T0 changes temperature without recomputing M_CEA or gamma | reproduced | Static source inspection confirms temperature-only scaling. Current G4a evidence config uses eta_T0=1. |
| 10 | Blowdown continues to use the combustion HP property table | reproduced | Static source inspection confirms the same pressure-dependent HP table remains active through blowdown. Root correction is deferred to Phase B. |
| 11 | The fix strongly changes MEOP/burn time while impulse or Isp may be weak discriminators | not_reproduced | Quantitative claim requires the post-fix run. Current pre-fix MEOP is 10.547965 MPa; no corrected local result exists yet. |
| 12 | Selftest can omit checks rather than reporting skip/fail | reproduced | Conditional test structure is confirmed in source. Exact 13-versus-15 output reproduction is deferred to the A-5.1-3 test-count contract. |

## Summary

- reproduced: 8
- not_reproduced: 3
- environment_dependent: 1

A classification of `not_reproduced` here does not reject the claim. It means the claim requires a corrected build or external comparison that is deliberately unavailable at this pre-fix HEAD.
