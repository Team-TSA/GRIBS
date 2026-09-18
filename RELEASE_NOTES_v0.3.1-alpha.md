# GRIBS v0.3.1-alpha

This alpha release introduces a unified JSON-based configuration system for GRIBS.

## Major changes

- Moved all standard user inputs from the Python source to `gribs_config.json`.
- Added selectable thermochemistry backends:
  - `cea2`
  - `legacy_fit`
- Integrated propellant composition and constituent properties into the unified JSON configuration.
- Added propellant-density calculation from constituent mass fractions and densities.
- Added packing-fraction correction for bulk propellant density.
- Added JSON configuration for:
  - propellant burn law
  - erosive burning
  - grain geometry
  - nozzle geometry and efficiency
  - flow separation
  - throat erosion
  - ambient and initial conditions
  - igniter operation
  - CEA2 execution and cache settings
  - manual thermodynamic-property correlations
  - numerical solver and blowdown
  - output directory and filenames
- Added configuration provenance, resolved inputs, and derived values to `summary.json`.
- Preserved compatibility with Visual Studio Code Run.
- Preserved script-relative discovery of external CEA files.
- Retained the existing internal-ballistics and numerical calculation logic.

## External NASA CEA dependency

> [!IMPORTANT]
> The external NASA CEA executable and database files are not included in this repository or release.

To use the `cea2` thermochemistry backend, users must separately obtain compatible versions of the following files:

- `fcea2`
- `thermo.lib`
- `trans.lib`

Place these files in the same directory as the GRIBS Python file.

A typical local directory structure is:

```text
GRIBS/
├── GRIBS_v0.3.1-alpha.py
├── gribs_config.json
├── fcea2
├── thermo.lib
├── trans.lib
└── results/
```

On Linux, the external CEA executable may require execute permission:

```bash
chmod +x fcea2
```

NASA CEA is not developed, maintained, or distributed by the GRIBS project. Users are responsible for obtaining NASA CEA from an authorized source and complying with all applicable license, redistribution, usage, and export-control requirements.

Official NASA CEA project:

https://github.com/nasa/CEA

NASA CEA installation documentation:

https://nasa.github.io/cea/installation.html

## Running without NASA CEA

GRIBS can run without the external CEA executable and database files.

Set the following option in `gribs_config.json`:

```json
{
  "thermochemistry": {
    "backend": "legacy_fit"
  }
}
```

In this mode, GRIBS uses the manually specified propellant density and pressure-dependent thermodynamic-property correlations in the `legacy_fit` section.

The supplied example coefficients are formulation-specific and must not be assumed valid for other propellants.

## Running GRIBS

Place the following files in the same directory:

```text
GRIBS_v0.3.1-alpha.py
gribs_config.json
```

Run:

```bash
python3 GRIBS_v0.3.1-alpha.py
```

An alternative complete configuration can be specified with:

```bash
python3 GRIBS_v0.3.1-alpha.py --config another_config.json
```

## Validation status

The unified JSON configuration and `legacy_fit` execution paths were tested through complete simulation and output generation.

The following output files were successfully generated:

- `ballistics.png`
- `time_history.csv`
- `summary.txt`
- `summary.json`

The CEA2 execution path requires a compatible external `fcea2` executable and CEA database files.

## Alpha-release notice

This version is an alpha release. Results must be independently validated before engineering or experimental use.

Material properties, burn-rate coefficients, constituent densities, packing fraction, nozzle-loss coefficients, geometry, and initial conditions must correspond to the actual propellant and motor being analyzed.
