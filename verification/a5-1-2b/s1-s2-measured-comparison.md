# A-5.1-2b measured S-1/S-2 comparison

Environment: Python 3.14.4, SciPy 1.18.1, CEA package/library 3.3.4.

| rtol | selftest | production exit | final mass rel. error | EOS-history max rel. error | post-run verdict | output files |
|---:|:---:|---:|---:|---:|:---:|---:|
| 1e-7 | PASS | 1 | 2.94e-09 | 3.05e-05 | FAIL | 0 |
| 1e-10 | PASS | 0 | 3.48e-11 | 3.60e-07 | PASS | 4 |

Both production runs reused cache key `62977944a749076ca9a584b2`.
The first 1e-7 selftest created that cache after a cache miss.
Existing top-level results files remained SHA-256 identical.
Both derived configurations were migrated from schema 0.4.4-alpha to 0.5.0-alpha.

Evidence-storage note: the generated production history originally used
CRLF line endings. Before committing the evidence, it was normalized to
LF so that the repository-wide `git diff --check` gate remains valid.
CSV field values, row order, and column order were not changed. The
evidence SHA-256 manifest records the normalized committed artifact.

Raw-log storage note: three generated logs originally contained one
extra blank line at EOF. Before finalizing the evidence commit, only that
trailing blank line was removed so that the repository
`git diff --check` contract is satisfied. Log messages, ordering,
numerical values, and exit-status records were not changed. The evidence
SHA-256 manifest records the normalized committed logs.
