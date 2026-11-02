# ADR: Validation of the CLI `--backend` override

- Status: Accepted
- Date: 2026-11-02
- Scope: configuration loading and validation only. No physical, numerical, tolerance,
  default-value, output-format, or cache schema/key change.
- Related contracts:
  - Phase A-5 reference backend: `verification/a5-reference/A5_DESIGN.md` §5 and its CLI
    paragraph, and `A5_COMPLETION.md` §2.5.
  - `_backend_choice`: an explicit backend is required and is never assumed or
    substituted.

## 1. Context

`--backend NAME` overrides `thermochemistry.backend` for one run. Before this change, the
override was assigned to the runtime `Config` after the document had been validated for the
backend written in the file. Rules that depend on the backend were therefore never applied
to the backend that actually ran:

- the legacy/homogeneous exclusion: `cea_legacy_executable` cannot run
  `two_phase_model.mode = "homogeneous_equilibrium"`;
- the sections and values each backend requires;
- the check of the legacy executable and data directory;
- the source of the chamber-property cache settings, which came from the section of the
  backend written in the file.

As a result, `--validate-config` reported such combinations as valid (exit code 0).

## 2. Decision

### 2.1 Two-stage validation

When `--backend` is given, `load_json_configuration` runs `configuration_from_document`
twice:

1. **On the loaded document by itself** (after any schema migration), with the same rules
   as without the override. If `configuration_from_document` rejects the file, the override
   does not make it acceptable, and the error is the same as without the override. For
   example, the file must name a valid backend of its own and meet that backend's own
   section rules.
2. **On the effective document.** This is a copy whose `thermochemistry.backend` is
   replaced by the override.
   - Every rule that depends on the backend is applied to the backend that will run.
   - The cache settings are taken from that backend's section; `reference_constant` uses
     no cache.

The value checks in `validate(Config)` run once, after both stages, on the effective
runtime configuration. They are not repeated for the file by itself. They are the same
checks as before, for example `R_i0 < R_p` and `p_a >= p_fit_min`.

### 2.2 Order of reported errors

Validation stops at the first error, so each run reports one problem. When a configuration
has several problems, the validation order decides which one appears first:

- **An error that `configuration_from_document` finds in the file by itself is reported
  first.** It is the same error as without the override, and errors in the effective
  combination (stage 2) are not reported until it is fixed.
- **This does not hold for errors that only `validate(Config)` detects.** Suppose a file's
  only problem is such a value check, for example `p_a < p_fit_min`, and the effective
  combination is also invalid, for example a homogeneous file overridden to
  `cea_legacy_executable`. Then the stage-2 error is reported first. Without the override,
  the value-check error is reported. Both cases exit with code 2.
- Within each `configuration_from_document` run, the existing order applies: section
  structure and unknown keys, then the reactants, the backend value, the legacy/homogeneous
  exclusion, the pressure range, the backend sections, and so on.
- Fixing one problem may reveal another. Only the exit code (2) and the first message are
  defined. A complete list of problems, and their order, is not part of this contract.

### 2.3 `reference_constant` values not stated in the file

When the override is `reference_constant`, two different sources supply values that the
file does not state. Both are kept.

**Case 1: the file has no `thermochemistry.reference_constant` section.**

- **Behaviour:** the documented `Config` defaults are used: `temperature_K = 3000.0`,
  `gamma = 1.2`, and `molecular_weight_kg_kmol = 25.0`.
- **Source:** the documented Phase A-5 contract. `A5_DESIGN.md` states: "When overriding a
  configuration that does not contain the explicit reference block, the documented
  `Config` defaults are used". The test
  `tests/reference/test_reference_constant_backend.py::test_cli_backend_override_uses_documented_reference_defaults`
  pins it.

**Case 2: the file is written for another backend and has a partial or empty section.**

- **Behaviour:** the keys that are present are used, and the missing keys take the same
  defaults.
- **Source:** behaviour that existed before this change and is not documented by A-5. The
  keys were optional when the file's backend was not `reference_constant`. This change
  keeps the behaviour, this ADR records it, and `tests/test_cli_backend_override_edge.py`
  pins it.

A file written for `reference_constant` must already have a complete section, because
stage 1 checks it. Without an override, a partial section in such a file is rejected, as
before.

### 2.4 Recorded configuration

The meaning of the recorded configuration does not change:

- `original_configuration` and `Configuration.document` hold the loaded document after any
  migration and before the override. Because a migration may have been applied, this is not
  always the file as written.
- `resolved_configuration` holds the effective runtime configuration, with the overridden
  backend.
- `configuration_migration.notes` holds the migration notes, followed by
  `thermochemistry.backend overridden on the command line to '<name>'.` The text is
  unchanged.

## 3. Consequences

| Input | Before | After |
|---|---|---|
| No override, or the override equals the file's backend | — | Same outcome |
| Valid override to another backend | exit 0 | exit 0; cache settings come from the target backend's section |
| `configuration_from_document` accepts the file, but the effective combination is invalid (for example, a homogeneous file overridden to `cea_legacy_executable`, or a missing section for the target backend) | exit 0 with `--validate-config`; at run time, exit 3 or silent defaults | exit 2 at validation |
| `configuration_from_document` rejects the file by itself, with any override | exit 2 | exit 2, with the same first error line |
| The file's only problem is a `validate(Config)` value check | exit 2 | exit 2; if the effective combination is also invalid, that error may be reported first (§2.2) |
| Override to `reference_constant`, with no section, or with a partial section in a file written for another backend | exit 0 with defaults | exit 0 with defaults (unchanged) |

Exit code 2 is the existing configuration-error code. This change adds no new message text.

**Numerical results under the conditions tested.** Before and after the change, the
following values were bit-identical:

- `burn_time_s` and `impulse_total_Ns`;
- for the legacy runs, also `p_max_Pa` and the numeric columns of the property table.

The conditions tested were:

- **Golden g2 and g1 inputs:**
  - g2 without an override;
  - g2 with `--backend cea_python`;
  - g2 with a reference section, overridden to `reference_constant`;
  - g1 overridden to `cea_python`;
  - g2 without a reference section, overridden to `reference_constant`.
- **Migrated documents:** v0.3.0, v0.4.3, and v0.4.4 documents overridden to
  `reference_constant`.
- **Legacy backend:** one composition with `single_phase_legacy`, with and without an
  override, run in an isolated sandbox with a substitute `libgfortran`.

These results are evidence for these conditions only. They do not prove identity in
general, and the legacy results say nothing about the original fcea2 runtime.

## 4. Rejected alternatives

- **Validate only the effective document.** This accepts files that are invalid by
  themselves once they are overridden. One candidate of this kind also bypassed
  `_backend_choice`.
- **Reject a missing or partial `reference_constant` section under an override.** This
  breaks the documented Phase A-5 contract and its tracked test, and removes the existing
  partial-section behaviour.
- **Re-run the backend-dependent checks after the old post-hoc assignment.** This
  duplicates the rules in a second place.

## 5. Verification

- New tests:
  - `tests/test_cli_backend_override.py` (13 tests)
  - `tests/test_cli_backend_override_edge.py` (41 tests)
  - `tests/test_cli_backend_override_boundary.py` (20 tests)
- Fixtures, taken from git history, in `tests/fixtures/cli_backend_override/`:
  `gribs_config.json` at commits `1e50615` (schema 0.3.0) and `b1c9bfd` (schema 0.4.3).
- The full test suite was run on two different trees. The results are not comparable,
  because the trees contain different code and tests:
  - **The committed tree** (the parent commit plus this change only): 169 passed. The parent
    alone gives 95 passed; the difference is the 74 new tests.
  - **A copy of the developer working tree** at the time of the change. That tree also
    contains other uncommitted work and additional test files. With this change, 435
    passed, and every test module loaded the changed program file.
