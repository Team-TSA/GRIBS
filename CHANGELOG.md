# Changelog

All notable changes to GRIBS are documented in this file.

## [v0.3.1-alpha] - 2026-09-18

### Added

- A unified `gribs_config.json` file for all standard user inputs.
- Selection between `cea2` and `legacy_fit` thermochemistry backends.
- Propellant reactant composition, constituent density, and initial-temperature input in the unified JSON configuration.
- Ideal additive-volume mixture-density calculation:
  `rho_ideal = 1 / sum(w_i / rho_i)`.
- Packing-fraction correction for bulk propellant density.
- JSON configuration for the burn law, erosive burning, grain geometry, nozzle geometry, flow separation, throat erosion, environment, igniter, thermochemistry, solver, blowdown, and outputs.
- Configuration source, original configuration, resolved runtime configuration, and derived inputs in `summary.json`.
- Script-relative discovery of the default JSON file and local external CEA runtime.
- Explicit reporting of ideal mixture density, packing fraction, and bulk propellant density.

### Changed

- Standard input management moved from source-level Python edits to a complete JSON configuration.
- CEA reactants moved from a separate reactant JSON file into `gribs_config.json`.
- CEA2 became a selectable JSON backend rather than a source-level selection.
- The runtime `Config` dataclass is now an internal resolved-configuration container rather than the standard user input surface.
- Output directory and output filenames became configurable through JSON.
- Documentation and source comments were updated for the v0.3 configuration architecture.

### Fixed

- Added the missing return value from CEA batch parsing.
- Added script-relative lookup for `fcea2` to support Visual Studio Code Run and arbitrary shell working directories.
- Added explicit executable-permission diagnostics on POSIX systems.
- Reduced each legacy CEA request to eight pressure points for compatibility with the tested `fcea2` output limits while retaining 81 total grid points.
- Removed obsolete source-level command-line override examples.
- Removed the obsolete `--list_params` documentation and residual CLI handling.
- Corrected duplicate source section numbering.
- Clarified that CEA product-gas density is not solid grain bulk density.

### Preserved

The following established calculation logic was retained:

- unsteady chamber pressure equation
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
- Burn-rate and loss coefficients remain formulation- and hardware-specific empirical inputs.
