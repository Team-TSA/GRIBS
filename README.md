# GRIBS

GRIBS is a JSON-configured transient internal-ballistics simulator for solid rocket motors. It models cylindrical-bore grain regression, unsteady chamber pressure, choked and unchoked nozzle flow, optional flow separation, optional erosive burning, optional throat erosion, and post-burnout blowdown.

**Current release:** v0.6.2-alpha — patch release of v0.6.0-alpha "Keystone" (Consistent Thermochemistry Foundation)  
**Release status:** Alpha research software. Physical validity has not been established; the comparison with the CEA reference `c*` (gate G4a) has not passed, and reported-data / experimental reproduction (G4b) and the release gate (G5) are not complete. Results must not be used for production, safety-critical, or design decisions without independent validation. See [Validation status](#validation-status) and [Known limitations](#known-limitations).

## Highlights

- Complete user configuration through `gribs_config.json`
- Direct integration with the official NASA CEA Python package through `import cea`
- Transitional support for an external `fcea2` executable
- Explicit thermochemistry backend selection with no silent fallback
- Pressure-dependent chamber properties tabulated and interpolated with PCHIP in `ln(p)`
- Constituent-density-based propellant mixture-density calculation
- Packing-fraction correction for bulk propellant density
- Homogeneous-equilibrium two-phase chamber model on a total-product mass basis (condensed products carried with the gas)
- Numerical self-tests on request (`--selftest`)
- Atomic writes of the chamber-property cache, with a structured error on write failure
- PNG, CSV, text, and JSON outputs
- Detailed configuration, environment, cache, and thermochemistry provenance in `summary.json`
- Optional theoretical CEA rocket-performance diagnostics for comparison only

## What changed in v0.6.2-alpha

v0.6.2-alpha is a patch release of v0.6.0-alpha "Keystone". The program version is 0.6.2-alpha and the program file is `GRIBS_v0.6.2-alpha.py`; the configuration schema is unchanged at `0.5.0-alpha`. Numerical results are unchanged. See `CHANGELOG.md` for the full list.

- **High-condensate warning.** Runs whose thermochemistry table reaches `Yg_min <= 0.80` with model status `approximate_with_limitations` also report `W_HIGH_CONDENSATE_G4B_UNVALIDATED`: reproduction of reported or experimental data (G4b) has not been done for that band. `model_status`, the run verdict and the exit codes are unchanged (see [Known limitations](#known-limitations)).
- **Validation status.** The README states that G0 was approved for the v0.6.1-alpha release commit `cf91f97`.

## What changed in v0.6.1-alpha

v0.6.1-alpha is a patch release of v0.6.0-alpha. The program version is 0.6.1-alpha and the program file is `GRIBS_v0.6.1-alpha.py`; the configuration schema is unchanged at `0.5.0-alpha`. The simulation code is unchanged apart from the version identifiers and one source comment. See `CHANGELOG.md` for the full list.

- **Default configuration.** `gribs_config.json` uses graphite `C(gr)` instead of `C`, which in the NASA CEA database is gaseous atomic carbon. Results of the default configuration change accordingly (see [Propellant composition and bulk density](#propellant-composition-and-bulk-density)).
- **Golden cases.** g7 and g8 are g3 and g6 with `C(gr)`. The values of g1–g6 are unchanged; g3 and g6 keep `C` as numerical regression references.
- **Documentation.** The outside-range setting is named by its JSON key `thermochemistry.outside_range_policy`; GRIBS does not check whether a composition is plausible.

## What changed in v0.6.0-alpha

v0.6.0-alpha is the first public release after v0.4.3-alpha. The intermediate v0.5.0-alpha was an internal Phase A baseline and was not published. The program version is 0.6.0-alpha; the configuration schema is unchanged at `0.5.0-alpha`. See `CHANGELOG.md` for the full list.

- **Two-phase model B and corrected mass basis.** The conserved chamber mass is the total product mass (gas + condensed). The equation of state acts on the gas phase only, and the CEA molecular weight is used as total mass per kmol of gas, without a second gas-fraction factor. Thermochemistry cache columns and cache identity were renamed accordingly; caches written under the old mass basis are not silently reused.
- **Validation and status reporting.** `validation.model_status` no longer reports `validated_in_scope`; statuses state what they are based on. Model-applicability checks fail closed on missing or non-finite phase fractions.
- **Input robustness.** NaN and infinite numbers are rejected in configuration files and in `Config` values (including NumPy and `Decimal` values), before any calculation.
- **Diagnostics.** A failed CEA theoretical diagnostic no longer stops the run or the summary; its status is reported.
- **Chamber-property cache.** The cache is written atomically; a write failure stops the run with exit code 3 (`E_THERMO_CACHE_WRITE_FAILED`). See [Chamber-property cache writes and concurrent runs](#chamber-property-cache-writes-and-concurrent-runs) for the operating assumptions and what is not tested.
- **Verification assets.** Golden regression cases (`cases/golden/`), reference problems, a mass and energy ledger, a convergence harness, and the A-5.1 mass-basis evidence (`verification/`).
- **CLI.** Self-tests run only with `--selftest`; normal runs no longer run them first. The `--backend` override is validated in two stages (see below).

Not included in this release (see [Known limitations](#known-limitations)): stopping at the property-table range with a dedicated exit code (M3), the safe output mode, the Phase B physical extensions (B1), and the Phase B comparison layers (B2-C / B2-D).

## What changed from the v0.3 series (v0.4 series)

The v0.4 series retains the established GRIBS unsteady internal-ballistics, grain-regression, burn-rate, nozzle-flow, throat-erosion, and blowdown models. The main changes are in thermochemistry integration, configuration validation, reproducibility, diagnostics, and provenance reporting.

- Added the recommended `cea_python` backend using the official NASA CEA Python package.
- Renamed the external executable backend from `cea2` to `cea_legacy_executable`.
- Removed the manual `legacy_fit` gas-property backend.
- Made thermochemistry backend selection explicit and mandatory.
- Standardized bulk propellant density calculation from constituent mass fractions, constituent solid densities, and packing fraction.
- Added stronger thermochemistry cache identity and validation.
- Added strict JSON validation, including rejection of unknown keys.
- Added optional theoretical CEA `c*`, `Cf`, and `Isp` diagnostics for comparison only.
- Added `--backend`, `--validate-config`, `--dump-thermo`, `--no-cache`, and `--version` command-line options.
- Expanded execution and thermochemistry provenance in `summary.json`.

The update from the supplied v0.4.2-alpha implementation to v0.4.3-alpha changes version identifiers and corrects comments and user-facing documentation. It does not change the core calculation logic.

## Repository files

```text
GRIBS/
├── GRIBS_v0.6.2-alpha.py
├── gribs_config.json
├── requirements.txt
├── requirements-dev.txt
├── README.md
├── CHANGELOG.md
├── MODEL_LIMITATIONS.md
├── LICENSE
├── CITATION.cff
├── CONTRIBUTING.md
├── SECURITY.md
├── cases/            golden and reference input cases
├── examples/         example configuration
├── tests/            pytest suite
├── tools/            verification tools
├── verification/     recorded verification evidence
└── results_golden_g1 … g6/   recorded golden outputs (regression references)
```

The repository root also keeps the patch files (`*.diff`) of earlier Phase A steps as historical records.

External CEA runtime files are not distributed by GRIBS.

## Requirements

- Python 3.11 or later. The GitHub CI runs the test suite on Python 3.11, 3.12 and 3.13. Python 3.14 is checked locally only (it is the environment of the recorded G0 evidence; see [Validation status](#validation-status)); it is not part of the CI matrix.
- NumPy
- SciPy
- Matplotlib
- Official NASA CEA Python package, API version 3.0 or later, for the recommended `cea_python` backend

Create and activate a virtual environment, then install the dependencies:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

Check the installed packages:

```bash
python -m pip check
python -c "import numpy, scipy, matplotlib, cea; print(cea.__version__)"
```

## Quick start

Place `GRIBS_v0.6.2-alpha.py` and `gribs_config.json` in the same directory.

Validate the configuration without initializing CEA:

```bash
python GRIBS_v0.6.2-alpha.py --config gribs_config.json --validate-config
```

Run the numerical self-tests:

```bash
python GRIBS_v0.6.2-alpha.py --config gribs_config.json --selftest
```

Run the simulation:

```bash
python GRIBS_v0.6.2-alpha.py --config gribs_config.json
```

Display the program version:

```bash
python GRIBS_v0.6.2-alpha.py --version
```

The default configuration path is resolved relative to the Python file. This allows the program to be started with the Visual Studio Code Run command even when the terminal working directory differs from the source directory.

## Configuration

All standard user inputs are stored in `gribs_config.json`. The main sections are:

- `propellant`
- `grain`
- `nozzle`
- `environment`
- `igniter`
- `thermochemistry`
- `solver`
- `output`
- `two_phase_model`

Field names include SI units where applicable. Edit the JSON configuration rather than the runtime `Config` dataclass in the Python source.

`thermochemistry.outside_range_policy` (required) selects how chamber properties are treated outside the tabulated pressure range `thermochemistry.pressure_range`: `"clamp"` holds them at the nearest table limit, and `"extrapolate"` extrapolates the table. The supplied configurations use `"clamp"`. In the Python source this setting is the `Config` field `property_policy`; that name is not a JSON key, and a configuration that contains it is rejected as an unknown key.

The current configuration schema is:

```json
{
  "schema_version": "0.5.0-alpha"
}
```

The configuration schema did not change in v0.6.0-alpha, v0.6.1-alpha or v0.6.2-alpha: it is still `0.5.0-alpha`, so v0.5.0-alpha configurations are used as they are. Configurations with schema `0.4.4-alpha`, `0.4.3-alpha` or `0.3.0-alpha` are migrated automatically when loaded where possible, and the migration notes are recorded in `summary.json`; a configuration that selects a removed backend is rejected. Other schema values are rejected.

## Thermochemistry backend selection

The backend must be selected explicitly. GRIBS does not choose a default backend and does not silently substitute another backend if the selected backend is unavailable.

### Recommended: official NASA CEA Python package

```json
{
  "thermochemistry": {
    "backend": "cea_python"
  }
}
```

The `cea_python` backend:

- imports the official package as `cea`
- performs HP equilibrium calculations through the Python API
- does not launch an external executable
- does not require `thermo.lib` or `trans.lib` beside the GRIBS script
- records package, library, database, and cache provenance

Install the backend with:

```bash
python -m pip install "cea>=3.0"
```

Official NASA CEA project:

- https://github.com/nasa/CEA
- https://nasa.github.io/cea/

### Transitional: external CEA executable

```json
{
  "thermochemistry": {
    "backend": "cea_legacy_executable"
  }
}
```

The `cea_legacy_executable` backend retains the v0.3 external `fcea2` workflow for regression comparison during the alpha series. It requires a compatible executable and thermodynamic database files.

When using this backend, place the following files in the same directory as the selected GRIBS Python file unless alternative paths are configured:

```text
fcea2
thermo.lib
trans.lib
```

On Linux, the executable may require permission:

```bash
chmod +x fcea2
```

NASA CEA is not developed, maintained, or distributed by the GRIBS project. Users are responsible for obtaining compatible CEA software and databases from an authorized source and complying with the applicable terms.

### Removed backend

The v0.3 `legacy_fit` backend is no longer supported. GRIBS does not use manually fitted pressure correlations as a substitute for CEA.

## Propellant composition and bulk density

Reactant names, mass percentages, constituent solid densities, and initial temperatures are read from `propellant.reactants`.

Example:

```json
{
  "name": "NH4NO3(IV)",
  "wt_percent": 88.7,
  "temperature_K": 298.15,
  "density_kg_m3": 1720.0
}
```

Requirements:

- Reactant names must exist in the thermochemistry database used by the selected backend. Use the name of the intended phase: in the NASA CEA database `C` is gaseous atomic carbon (about +59.67 MJ/kg at 298.15 K) and graphite is `C(gr)`; likewise `AL` is aluminium vapour and `AL(cr)` is solid aluminium. GRIBS does not check whether the chosen species are plausible propellant ingredients.
- `wt_percent` values must sum to 100.
- Every `density_kg_m3` value must be positive and represent the relevant solid constituent density.
- Every `temperature_K` value must be positive.
- `packing_fraction` must satisfy `0 < packing_fraction <= 1`.

GRIBS calculates the ideal additive-volume mixture density as:

```text
rho_ideal = 1 / sum_i(w_i / rho_i)
```

The bulk propellant density used by the internal-ballistics model is:

```text
rho_p = rho_ideal * packing_fraction
```

Equivalently:

```text
rho_p = packing_fraction / sum_i(w_i / rho_i)
```

CEA combustion-product gas density is never used as the solid propellant density.

This mixing rule is an approximation. Constituent interaction, processing, porosity, phase changes, and non-additive volume effects may require measured bulk-density data or a more specialized model.

## Inputs not determined by CEA

CEA equilibrium calculations do not determine all solid-propellant motor inputs. Users must provide defensible values for quantities including:

- burn-rate coefficient
- pressure exponent
- burn-rate temperature sensitivity
- erosive-burning coefficients
- constituent solid densities
- packing fraction
- grain geometry
- nozzle discharge coefficient
- thrust efficiency
- throat-erosion parameters

These inputs should be based on appropriate material data, experiments, calibration, or validated references.

## CEA theoretical rocket diagnostics

When enabled for `cea_python`, GRIBS can report theoretical CEA rocket-performance values, including theoretical `c*`, `Cf`, and `Isp`.

These values are for comparison only. They are not fed back into the GRIBS unsteady internal-ballistics or nozzle-flow calculations. They must not be confused with the GRIBS effective quantities `cstar_eff_m_s`, `CF_eff`, and `Isp_s`.

## Command-line options

```text
--config FILE          Select a complete JSON configuration file
--backend NAME         Override the configured thermochemistry backend for one run
--selftest             Run numerical self-tests only
--validate-config      Validate the JSON configuration and exit
--dump-thermo FILE     Write the chamber-property table and exit
--no-cache             Force rebuilding of the chamber-property table
--allow-verification-failure
                       Return 0 only for a mass-ledger-only verification failure
--version              Display the GRIBS version
```

### Process exit codes

| Code | Meaning |
|---:|---|
| 0 | Run passed, or a mass-ledger-only failure was explicitly allowed. The diagnostic outputs still record the failed ledger. |
| 1 | Run verification gate or standalone self-test failed; also used for model-applicability rejection. |
| 2 | Configuration is invalid. |
| 3 | Thermochemistry backend initialization failed, including a failure to write the chamber-property cache (`E_THERMO_CACHE_WRITE_FAILED`). |
| 4 | Simulation/integration failed, including an integration that exceeds a safety limit (`E_INTEGRATION_BUDGET_EXCEEDED`: more than 2,000,000 right-hand-side evaluations, more than 300 s of wall time, or more than 10,000 solver restarts at the property-table limits). |
| 5 | Mass-ledger verification failed while the other run-verification gates passed. |

`--allow-verification-failure` applies only when the mass ledger is the sole
failed run-verification gate. It does not override post-run consistency or
model-applicability failures. In a failed normal run, the CSV, text summary,
and JSON summary are written before the nonzero exit; the ledger warning and
residuals remain recorded even when the flag permits exit code 0.

### Backend override (`--backend`)

`--backend NAME` is validated in two stages:

1. The file is read on its own, with the same rules as without the override. The file must
   name a valid backend and satisfy that backend's section rules.
2. The combination with `NAME` is read. For example, `cea_legacy_executable` cannot run
   `two_phase_model.mode = "homogeneous_equilibrium"`, and `NAME`'s own section must be
   present.

Value checks such as `R_i0 < R_p` and `p_a >= p_fit_min` run once afterwards, on the
effective configuration.

Both `--validate-config` and a normal run report an invalid configuration as a
configuration error (exit code 2). Validation stops at the first error:

- An error found while reading the file on its own (stage 1) is reported first, with the
  same message as without the override.
- If the file's only problem is a value check and the combination is also invalid, the
  combination error may be reported first.
- After one problem is fixed, another may be reported.

The cache settings are taken from the section of the backend that runs.

When the override is `reference_constant`:

- If the file has no `reference_constant` section, the documented defaults are used
  (3000 K, 1.2, 25 kg/kmol). This is the Phase A-5 contract.
- If the file is written for another backend and has a partial `reference_constant`
  section, the missing values take the same defaults. This keeps the earlier behaviour.

See `verification/cli/CLI_BACKEND_OVERRIDE_ADR.md`.

Example thermochemistry-table export:

```bash
python GRIBS_v0.6.2-alpha.py \
  --config gribs_config.json \
  --dump-thermo thermo_table.json
```

## Outputs

The default output directory is `results/`. Output names are configurable in `gribs_config.json`.

```text
results/
├── ballistics.png
├── time_history.csv
├── summary.txt
├── summary.json
└── thermo_cache/
```

- `ballistics.png`: combined publication-style result figure
- `time_history.csv`: sampled time history
- `summary.txt`: human-readable result summary
- `summary.json`: configuration, resolved runtime inputs, derived values, environment information, thermochemistry provenance, diagnostics, and results
- `thermo_cache/`: reusable pressure-dependent chamber-property tables

### Chamber-property cache writes and concurrent runs

With caching enabled, a run writes the cache only when it builds the table: when no cache
exists yet, or when it rebuilds the table (`rebuild_cache = true` or `--no-cache`). An
existing cache is replaced only in the second case; otherwise it is read and left
unchanged. GRIBS writes the table to a temporary file in the cache directory
(`.<name>.<pid>.<random>.tmp`), syncs it, publishes it under the cache name and syncs the
cache directory. Under the operating assumptions below, a run that reads the cache sees
either the previous complete file or the new complete file, never a partially written one.
Reading is not changed by this: a run checks that the cache name is a regular file and then
opens it normally (following symbolic links), without the special-file protection that the
writer has.

- A failure while writing the cache stops the run with exit code 3 and the structured
  error `E_THERMO_CACHE_WRITE_FAILED`. Its `context` gives `phase`, `step`, `replaced`
  (whether this run changed the cache name) and `temp_cleanup` (what happened to this
  run's temporary file). There is no fallback to an in-place write.
- If syncing the cache directory is not supported (`ENOTSUP` / `EOPNOTSUPP`), the run
  continues and records the warning `W_THERMO_CACHE_DURABILITY_NOT_GUARANTEED`.
- Survival of the cache after a crash or power loss is not guaranteed. In particular, when
  GRIBS creates a new cache directory, the parent of that directory is not synced, so the
  new directory and its cache may not survive a power loss.

Operating assumptions:

1. Restrict access to the cache directory: only trusted users of this cache may write to
   it. Temporary files are visible in the directory while a table is written; restrict
   read and search access to the directory if their content must not be readable by
   others.
2. Use the same permission settings (umask, default ACL of the directory) for all runs
   that write to the same cache directory. A setgid directory aligns the group only, not
   the mode.
3. Runs that write the same cache (same key) must produce the same content from the same
   inputs.

Concurrent runs under these assumptions (verified for one user on ext4, with threads in
one process and with 2 to 4 synchronised processes):

- If every run observed that the cache did not exist, the first run to publish keeps its
  file (inode and permissions). A later run with identical content adopts that file,
  removes its own temporary file and syncs the directory itself. New caches are published
  without replacing an existing name: by a hard link, or by
  `renameat2(RENAME_NOREPLACE)` where hard links are not available (Linux, glibc 2.28 or
  later). If neither is available, nothing is published and the run exits with code 3.
- A run that rebuilds the table while a cache already exists replaces it. The cache name
  then refers to a new file (new inode) with the same content; group, extended attributes
  and mode are copied from the existing cache, so permissions are not widened.
- Runs that write different content for the same cache are outside these assumptions.
  Depending on timing, a later run either stops with exit code 3 without replacing the
  cache, or replaces it.

Compatibility notes (differences from earlier versions):

- Exit code 3 is returned when rebuilding an existing cache that this user cannot read,
  when creating a new cache that this user could not read back (for example with umask
  0466), when rebuilding a cache inside a read-only directory, when a FIFO or another
  non-regular file is found where the writer checks the cache name (the writer does not
  wait on it), and when creating a new cache on a file system that supports neither hard
  links nor `renameat2(RENAME_NOREPLACE)`.
- If the cache name is a symbolic link when the run starts, the file it points to is
  rewritten (in the target's directory) and the link is kept. A symbolic link that appears
  at a checked name during the write makes the run stop with exit code 3. Device files
  have not been tested.
- The inode of the cache file changes when it is rebuilt, and a hard-linked alias is not
  updated. After a forced kill (SIGKILL) a temporary file may remain; it is not removed
  automatically.

Not tested: other users (UIDs), other file systems (for example NFS), other hosts, file
systems without hard links (simulated only), `security.*` / `trusted.*` extended
attributes, SELinux, device files (tested with a FIFO only), and power loss.

### Example output

The following figure is an example of the combined diagnostic output generated
by GRIBS. It shows the chamber-pressure and thrust histories, mass-flow balance,
grain-geometry evolution, ignition transient, nozzle operating point, grain
regression, and a summary of the principal calculated results.

![docs/images/ballistics-example.png](docs/images/ballistics-example.png)

*Example combined diagnostic figure generated by an earlier version (GRIBS v0.4.3-alpha). Click the
image to view the full-resolution version. The displayed values are from an
example calculation and are not validated design recommendations.*

## Validation status

GRIBS v0.6.2-alpha is alpha research software. **Its physical validity has not been established.** The evidence below fixes the reproducibility of the software and the internal numerical consistency of its use of NASA CEA within stated limits. It is not experimental validation, not independent verification, and not a qualification for any application. The gates G0–G5 are the project's Phase A acceptance steps.

- **G0 — evidence fixation (reproducibility): approved for revision `cf91f97` (the v0.6.1-alpha release commit) and, separately, for revision `4e7896a` (the v0.6.0-alpha release commit).** Earlier code revisions, including `6ec2e6b` and `236ad1e`, were approved separately. The approval of `4e7896a` is limited to that revision, the recorded environment and the following objects:
  - environment: Ubuntu 26.04.1 LTS, Linux 7.0.0-42, glibc 2.43, Python 3.14.4, cea 3.3.4 (with recorded library and data-file hashes), NumPy 2.5.3, SciPy 1.18.1, Matplotlib 3.11.2, pytest 9.1.1;
  - tests: 797 tracked tests passed;
  - golden cases g1–g6 rerun through the CLI: all 13 result metrics of each case pass the golden regression test described below; g1 and g3–g6 are identical to the earlier reproduction (time histories, figures and CEA caches byte-identical; `summary.json` and `summary.txt` identical after replacing run-directory paths); g2 differs only through the solver restart at the property-table limit (physical and timing metrics within 3.5e-10 of the golden values);
  - G1 probes: numerical difference 0 against the earlier record.

  The v0.6.1-alpha release commit `cf91f97` was reproduced separately in the same recorded environment (801 tracked tests passed; golden cases g1–g8 rerun through the CLI; G1 probes with numerical difference 0 against the `4e7896a` record), and G0 was approved for it with the same limits. G0 does not extend automatically to any later revision, including changes made after v0.6.1-alpha. The G0 decision for the v0.6.2-alpha release commit is stated in the v0.6.2-alpha release notes, because a commit cannot record the approval of itself. G0 does not cover physical validity, other kernels, other Python or dependency versions (including other `cea` versions; `requirements.txt` accepts `cea>=3.0`), other CPUs or floating-point environments, the external `fcea2` executable, or CI.
- **G1 — CEA semantics (internal consistency): approved within a limited scope.** It covers the numerical internal consistency of cea 3.3.4 outputs (mass-fraction sums, molecular-weight and density identities, species reconstruction, nozzle stations) for golden cases g2–g5 and two mid-condensed AP/Al evidence compositions (AP70/Al5/HTPB25 and AP70/Al7/HTPB23): 60 equilibrium points and 54 nozzle-station points. Conditions: some checks at a stricter tolerance (1e-8) remain recorded mismatches; some interpretations of CEA quantities are inferred and not confirmed against a CEA specification; G1 is not independent verification or physical validation.
- **G4a — comparison with the CEA reference `c*`: not passed (known limitation).** The criterion (relative difference within 2e-3 between the GRIBS constant-gamma nozzle `c*` and the CEA infinite-area-chamber `c*`) and its tolerance are unchanged. 9 of 18 points exceed it: g3 at 10/30/70 bar (−2.0e-3 to −2.1e-3), g4 at 70 bar (+2.05e-3), AP70/Al5 at 10 and 30 bar (+2.1e-3 to +2.6e-3), and AP70/Al7 at 10/30/70 bar (about −2.1% to −2.3%). The larger differences coincide with condensed species appearing or changing phase between the chamber and the throat; this is supported by intervention tests but not proven. The equilibrium-expansion alternative studied during Phase A is not part of this release.
- **G4b — reproduction of reported or experimental data: not done.** The complete inputs of the reported AP/Al cases, experimental data and an independent two-phase model are not available. Runs whose thermochemistry table reaches the high-condensate band (Yg_min <= 0.80, model status `approximate_with_limitations`) report the warning `W_HIGH_CONDENSATE_G4B_UNVALIDATED`; see [Known limitations](#known-limitations).
- **Golden regression test:** the physical and timing metrics are compared with the frozen golden values (relative 1e-6 / 1e-9). The three residual metrics (propellant mass balance, total mass consistency, maximum condensed mass fraction) are round-off-level numerical residuals that change by tens of percent with last-digit differences of the CEA table; they are checked against the program's own verification thresholds instead: propellant mass balance <= max(1e-9, 100 x solver rtol) (the mass-ledger gate), total mass consistency <= 1e-6 (the post-run consistency tolerance), maximum condensed mass fraction within max(1e-9 x |golden|, 1e-12) of the golden value. The test also checks the two-phase model mode and activity of each case. The golden values of g1–g6 are unchanged; g7 and g8 were added in v0.6.1-alpha. Golden cases g3 and g6 use the reactant `C` (gaseous atomic carbon, see [Propellant composition and bulk density](#propellant-composition-and-bulk-density)); they are kept unchanged as numerical regression references, and their results are not representative of a real propellant. Golden cases g7 and g8 are the same cases with graphite `C(gr)` (carbon density 2200 kg/m³ as in g3 and g6; the default configuration uses 2260 kg/m³). Their values are frozen in `cases/golden/golden_values.json`; unlike g1–g6 they have no recorded output directory (`results_golden_g7`, `results_golden_g8`).
- **G5 — release gate: not done.** The golden reference values in `cases/golden/` (g1–g8) and the recorded outputs in `results_golden_*` (g1–g6) were frozen before G5; g7 and g8 were added before G5 with the maintainer's approval. They are regression references, not validated or approved values.
- **Continuous integration:** the GitHub workflow runs the self-tests and the test suite on Python 3.11, 3.12 and 3.13. Python 3.14 is checked locally only. CI results are not part of the G0 evidence.

## Known limitations

- This is alpha-stage research software. It is not certified for safety-critical use and is not intended for production use; physical validity has not been established (see [Validation status](#validation-status)).
- The nozzle `c*` uses a constant-gamma isentropic relation. Against the CEA equilibrium reference it differs by up to about 2.3% in the tested condensed AP/Al compositions (G4a not passed).
- The chamber gas phase is treated as an ideal gas; non-ideal effects at high pressure are not modelled or bounded. The condensed-phase volume is neglected.
- Condensed-phase scope, by the minimum gas mass fraction `Yg_min` of the thermochemistry table (code constants, not linked to validation evidence):
  - `Yg_min <= 0.65`: `model_status = out_of_scope`; the run stops before the simulation with `E_OUT_OF_VALIDATED_RANGE`.
  - `0.65 < Yg_min <= 0.99`: `model_status = approximate_with_limitations`; the run continues and reports `W_HIGH_CONDENSATE_MODEL_LIMITATION`. This warning describes the physics the homogeneous-equilibrium model does not include (particle slip, finite-rate condensation, condensed-phase volume, phase-density lag).
  - `0.65 < Yg_min <= 0.80` (high-condensate band, for example AP/Al/HTPB with about 15–20 % aluminium): in addition, `W_HIGH_CONDENSATE_G4B_UNVALIDATED`. This warning is about validation evidence, not an additional model error: reproduction of reported or experimental data (G4b) has not been done for this band. It does not change `model_status`, the run verdict or the exit code.
- With `thermochemistry.outside_range_policy = "clamp"` (used by the supplied configurations), chamber properties outside the tabulated pressure range are held at the nearest table limit. Results in that region depend on this approximation; the run is not stopped. The pressure derivative of the properties jumps to zero at the table limits, so dp/dt is discontinuous there; the ODE solver is restarted at each crossing of a table limit so that it does not stall at the jump.
- GRIBS does not decide whether a composition is combustible, ignites or sustains combustion, and it does not warn about implausible compositions. The burn rate is the configured law, so the shape of the pressure and thrust histories is set mainly by the grain geometry and the burn law. With the geometry, burn law and nozzle fixed, the normalized pressure histories of different compositions nearly overlap: the composition changes mainly their level (through the chamber properties and `c*`), and their shape only slightly.
- Chemical equilibrium does not model finite-rate chemistry.
- The additive-volume density rule is an approximation.
- Burn-rate and erosive-burning coefficients are empirical inputs.
- The external executable backend depends on the user's compatible executable and database files.
- Theoretical CEA rocket diagnostics do not replace GRIBS results.
- The chamber-property cache relies on the operating assumptions in [Chamber-property cache writes and concurrent runs](#chamber-property-cache-writes-and-concurrent-runs). Survival after a crash or power loss is not guaranteed (the parent of a newly created cache directory is not synced), and other users (UIDs), other file systems, other hosts, file systems without hard links, `security.*` / `trusted.*` extended attributes, SELinux and device files have not been tested. Publishing a new cache requires hard links or Linux `renameat2` (glibc 2.28 or later).
- Not included in this release and not approved for use: stopping at the property-table range with a dedicated exit code (M3), the safe output mode (output-directory locking and overwrite refusal), the Phase B physical extensions (B1: energy / equation-of-state closure), and the Phase B CEA comparison layers (B2-C, B2-D). None of them is connected to the program.
- `MODEL_LIMITATIONS.md` summarises the model assumptions. It was written for the earlier v0.3 / v0.4 model and has been revised only minimally for v0.6.0-alpha, v0.6.1-alpha and v0.6.2-alpha; where it differs, this README and `CHANGELOG.md` take precedence.
- Model outputs are not substitutes for static firing tests, material characterization, or independent review.

## License

GRIBS is licensed under the MIT License.

Copyright (c) 2026 Hikaru Kurokawa.

The GRIBS license applies only to files distributed by the GRIBS project. External NASA CEA software and databases are not included in this repository and remain subject to their respective terms.

## Citation

Citation metadata is provided in `CITATION.cff`.

## Contributing and support

See `CONTRIBUTING.md` for contribution guidance. Use GitHub Issues for reproducible bug reports and feature proposals. See `SECURITY.md` for private security-reporting guidance.
