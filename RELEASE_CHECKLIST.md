# Public release checklist

Checklist for publishing a GRIBS release (current target: v0.6.0-alpha
"Keystone"). Each step that runs tests, changes history or publishes anything
needs the maintainer's explicit approval.

## 1. Decisions recorded before preparation

- [ ] Scope of the release decided (v0.6.0-alpha: G4a kept as a documented
      known limitation (G4a-0); M3 and the safe output mode deferred to the
      next release; `SCHEMA_VERSION` unchanged).
- [ ] Release name and title decided.

## 2. Release-preparation changes

- [ ] `PROGRAM_VERSION` updated; `SCHEMA_VERSION` changed only if the
      configuration schema changed (with a migration step).
- [ ] Program file renamed to `GRIBS_v<version>.py`; hard-coded references in
      tests and tools updated; historical text ("since v…") left unchanged.
- [ ] README: current release, file names, schema example, validation status,
      known limitations and supported Python versions are current.
- [ ] CHANGELOG: release section with date; earlier release sections unchanged.
- [ ] CITATION.cff: `version` and `date-released`.
- [ ] Golden inputs, golden values, recorded golden outputs and `verification/`
      are not changed by the preparation.

## 3. Verification of the release revision

- [ ] G0 reproduction on the exact release revision, outside the working
      checkout (export without `.git`): environment record, tracked tests,
      golden cases through the CLI, output comparison with the previous G0
      record (expected differences limited to the version strings and the
      figure), G1 probes. The maintainer decides the G0 approval.
- [ ] Python versions: the GitHub CI runs `--selftest` and `pytest -q` on
      Python 3.11, 3.12 and 3.13; Python 3.14 is checked locally (not in the
      CI matrix). Record skipped tests (for example ACL tests without
      `setfacl`) and the installed `cea` version for each run.
- [ ] Push the release branch and confirm that the GitHub CI passes on it before merging (golden residual metrics are checked against the verification thresholds; physical metrics against the frozen values).
- [ ] Long-running or stalled integrations: confirm that the integration safety limits (2,000,000 RHS evaluations, 300 s wall time, 10,000 table-limit restarts; exit 4 with E_INTEGRATION_BUDGET_EXCEEDED) are documented and not reached by the golden and example cases.
- [ ] Results from Python or dependency versions other than the G0 environment
      are reported as CI or reference results, not as G0 evidence.

## 4. Repository content and history

- [ ] Only tracked files are published: confirm with `git ls-files` that no
      unrelated working files are tracked.
- [ ] Review the Git history to be published for secrets and restricted
      material.
- [ ] Author information: confirm the author names and e-mail addresses that
      will be published
      (`git log --format='%an <%ae> | %cn <%ce>' main..HEAD | sort | uniq -c`).
      History is not rewritten for this.
- [ ] The internal technical manual and internal work records stay outside
      the public repository.
- [ ] Copyright ownership and the MIT License are confirmed.
- [ ] Default values are clearly identified as examples.

## 5. Backup, merge, tag and release

- [ ] Bundle of all local refs updated and verified by restoring it; external
      backup updated.
- [ ] Merge into `main` by fast-forward or an ordinary merge commit only (no
      squash, no rebase), so that the commit IDs referenced by the approval
      records stay reachable.
- [ ] CI passes on `main`.
- [ ] Annotated (or signed) tag `v<version>` on the G0-approved release
      commit.
- [ ] GitHub Release with the title and the release notes, stating the
      validation status and known limitations (physical validity not
      established; G4a not passed; G4b / G5 not complete; features not
      included).
- [ ] Recheck the public pages while signed out of GitHub.

## 6. Repository settings (once)

- [ ] Branch protection for `main` with required checks.
- [ ] Dependabot alerts, secret scanning, push protection and private
      vulnerability reporting where available.
