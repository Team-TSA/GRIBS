# Changelog

All notable changes to GRIBS are documented in this file.

## [v0.4.3-alpha] - 2026-09-20

### Changed

- Updated the public program identifier to `GRIBS v0.4.3-alpha`.
- Updated the configuration schema identifier to `0.4.3-alpha`.
- Updated command examples and user-facing version strings to use `GRIBS_v0.4.3-alpha.py`.
- Corrected comments and documentation so that they accurately match the implemented behavior.
- Corrected the solver-facing thermochemistry interface documentation to state that property evaluations made by the ODE must be side-effect free and inexpensive.
- Clarified that `cea_legacy_executable` is a transitional compatibility backend retained for regression comparison.
- No changes were made to the core calculation logic relative to the supplied v0.4.2-alpha implementation.

### Thermochemistry changes relative to the v0.3 series

- Added the recommended `cea_python` backend based on the official NASA CEA Python package.
- Renamed the former `cea2` backend to `cea_legacy_executable`.
- Removed the manual `legacy_fit` gas-property backend.
- Added explicit, mandatory backend selection with no default and no silent fallback.
- Added backend-specific cache identities and consistency checks.
- Added thermochemistry package, library, database, and cache provenance.
- Added physical validation of pressure-indexed chamber-property tables.
- Added optional theoretical CEA rocket-performance diagnostics for comparison only.

### Configuration

- Updated `schema_version` to `0.4.3-alpha`.
- Standardized bulk solid-propellant density calculation as:

  ```text
  rho_p = packing_fraction / sum_i(w_i / rho_i)
  ```

- Clarified that CEA combustion-product gas density is never used as solid propellant density.
- Added strict validation of required fields, data types, numeric ranges, backend-specific options, and unknown keys.
- Added automatic migration support for the recognized v0.3 schema while reporting every migration action.

### Command-line interface

- Added `--backend` to override the configured backend for one run.
- Added `--validate-config` to validate the JSON configuration without initializing thermochemistry.
- Added `--dump-thermo` to export the pressure-indexed chamber-property table.
- Added `--no-cache` to force rebuilding of the chamber-property table.
- Added `--version` to report the installed program version.

### Outputs and diagnostics

- Expanded `summary.json` with program, Python, platform, dependency, configuration, cache, and thermochemistry provenance.
- Added explicit separation between GRIBS effective performance values and theoretical CEA diagnostic values.
- Improved the output figure layout and automatic fitting of the key-results panel.
- Renamed the default cache output concept from the v0.3 CEA-specific cache to a backend-aware thermochemistry cache.

### Compatibility

- The v0.4 configuration schema is not a drop-in replacement for the v0.3 schema.
- Normal v0.4 configurations must select either `cea_python` or `cea_legacy_executable`.
- The old backend name `cea2` is recognized only by the explicit v0.3 migration path.
- The removed backend name `legacy_fit` is rejected and is never replaced automatically.
- The established unsteady chamber-pressure, grain-regression, burn-rate, nozzle-flow, throat-erosion, and blowdown models are retained.

### Validation

- Python syntax compilation completed successfully.
- JSON syntax validation completed successfully.
- `--version` reports `GRIBS 0.4.3-alpha`.
- `--validate-config` accepts the supplied `0.4.3-alpha` configuration.
- Full numerical results still require validation in an environment with the selected NASA CEA backend and all Python dependencies installed.

### Known limitations

- Alpha release, not certified for safety-critical use.
- Chemical equilibrium does not represent finite-rate chemistry.
- Burn-rate and loss coefficients remain formulation-specific and hardware-specific empirical inputs.
- The constituent-density mixing rule assumes additive constituent volumes before packing correction.
- External executable compatibility must be validated by each user selecting `cea_legacy_executable`.

## [v0.3.1-alpha] - 2026-09-18

### Added

- Added a unified `gribs_config.json` file for all standard user inputs.
- Added selection between `cea2` and `legacy_fit` thermochemistry backends.
- Added propellant reactant composition, constituent density, and initial-temperature input in the unified JSON configuration.
- Added ideal additive-volume mixture-density calculation:

  ```text
  rho_ideal = 1 / sum_i(w_i / rho_i)
  ```

- Added packing-fraction correction for bulk propellant density.
- Added JSON configuration for the burn law, erosive burning, grain geometry, nozzle geometry, flow separation, throat erosion, environment, igniter, thermochemistry, solver, blowdown, and outputs.
- Added configuration source, original configuration, resolved runtime configuration, and derived inputs to `summary.json`.
- Added script-relative discovery of the default JSON file and local external CEA runtime.
- Added explicit reporting of ideal mixture density, packing fraction, and bulk propellant density.

### Changed

- Moved standard input management from source-level Python edits to a complete JSON configuration.
- Moved CEA reactants from a separate reactant JSON file into `gribs_config.json`.
- Made CEA2 a selectable JSON backend rather than a source-level selection.
- Changed the runtime `Config` dataclass into an internal resolved-configuration container rather than the standard user input surface.
- Made the output directory and output filenames configurable through JSON.
- Updated documentation and source comments for the v0.3 configuration architecture.

### Fixed

- Added the missing return value from CEA batch parsing.
- Added script-relative lookup for `fcea2` to support Visual Studio Code Run and arbitrary shell working directories.
- Added explicit executable-permission diagnostics on POSIX systems.
- Reduced each legacy CEA request to eight pressure points for compatibility with the tested `fcea2` output limits while retaining 81 total grid points.
- Removed obsolete source-level command-line override examples.
- Removed obsolete `--list_params` documentation and residual CLI handling.
- Corrected duplicate source section numbering.
- Clarified that CEA product-gas density is not solid grain bulk density.

### Preserved

The following established calculation logic was retained:

- unsteady chamber-pressure equation
- cylindrical-bore grain geometry
- Saint-Robert burn law
- optional Lenoir-Robert erosive burning
- choked and unchoked nozzle branches
- optional flow separation
- optional throat erosion
- post-burnout blowdown
- pressure-dependent gas-property interpolation
- LSODA, BDF, and Radau integration options
- mass and impulse quadratures as ODE states
- numerical self-tests and verification metrics

### External dependency notice

The legacy NASA CEA executable and databases are not distributed by GRIBS. Users selecting the `cea2` backend must separately provide compatible copies of `fcea2`, `thermo.lib`, and `trans.lib` and must comply with the applicable terms.

### Known limitations

- Alpha release, not certified for safety-critical use.
- External CEA compatibility must be validated by each user.
- The constituent-density mixing rule assumes additive constituent volumes before packing correction.
- Burn-rate and loss coefficients remain formulation-specific and hardware-specific empirical inputs.
