# GRIBS

GRIBS is a JSON-configured internal-ballistics simulator for solid rocket motors. It models cylindrical-bore grain regression, unsteady chamber pressure, choked and unchoked nozzle flow, optional flow separation, optional erosive burning, optional throat erosion, and post-burnout blowdown.

> **Release status:** `v0.3.1-alpha` is a pre-release. Results must be independently validated before engineering or experimental use.

## Highlights

- Complete user configuration through `gribs_config.json`
- Optional external NASA CEA thermochemistry backend
- Manual pressure-dependent property correlations when CEA is not used
- Constituent-density-based propellant mixture-density calculation
- Packing-fraction correction for bulk propellant density
- Script-relative file discovery for Visual Studio Code Run compatibility
- Numerical self-tests before production calculations
- PNG, CSV, text, and JSON outputs
- Configuration provenance and derived inputs recorded in `summary.json`

## Repository files

```text
GRIBS/
├── GRIBS_v0.3.1-alpha.py
├── gribs_config.json
├── requirements.txt
├── README.md
├── CHANGELOG.md
├── LICENSE
├── CITATION.cff
├── CONTRIBUTING.md
├── SECURITY.md
└── .gitignore
```

External NASA CEA files are not included.

## Requirements

- Python 3.10 or later is recommended
- NumPy
- SciPy
- Matplotlib

Install Python dependencies with:

```bash
python3 -m pip install -r requirements.txt
```

## Quick start

Place `GRIBS_v0.3.1-alpha.py` and `gribs_config.json` in the same directory, then run:

```bash
python3 GRIBS_v0.3.1-alpha.py
```

An alternative complete configuration can be selected with:

```bash
python3 GRIBS_v0.3.1-alpha.py --config another_config.json
```

Run numerical self-tests only with:

```bash
python3 GRIBS_v0.3.1-alpha.py --selftest
```

The default `gribs_config.json` path is resolved relative to the Python file, so the program can be started with Visual Studio Code Run even when the terminal working directory differs from the source directory.

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

Field names include SI units where applicable. Users should edit the JSON file rather than the runtime `Config` dataclass in the Python source.

### Thermochemistry backend selection

Use external NASA CEA:

```json
{
  "thermochemistry": {
    "backend": "cea2"
  }
}
```

Use manually supplied pressure-dependent correlations:

```json
{
  "thermochemistry": {
    "backend": "legacy_fit"
  }
}
```

### CEA2 mode

When `thermochemistry.backend` is `cea2`, GRIBS reads the reactant names, mass percentages, constituent densities, and initial temperatures from the `propellant.reactants` array.

Each reactant entry has the following form:

```json
{
  "name": "NH4NO3(IV)",
  "wt_percent": 88.7,
  "temperature_K": 298.15,
  "density_kg_m3": 1720.0
}
```

Requirements:

- Reactant names must match names available in the user's compatible CEA thermodynamic database.
- `wt_percent` values must sum to 100.
- Every `density_kg_m3` value must be positive and must represent the relevant solid constituent density.
- Every `temperature_K` value must be positive.
- `packing_fraction` must satisfy `0 < packing_fraction <= 1`.

GRIBS calculates the ideal additive-volume mixture density as:

```text
rho_ideal = 1 / sum(w_i / rho_i)
```

The bulk propellant density used by the internal-ballistics model is:

```text
rho_p = rho_ideal * packing_fraction
```

This mixing rule is an approximation. Constituent interaction, processing, porosity, phase changes, and non-additive volume effects may require measured bulk-density data or a more specialized model.

### External NASA CEA dependency

> [!IMPORTANT]
> The external NASA CEA executable and database files are not included in this repository or release.

The legacy CEA runtime currently expected by the `cea2` backend uses these compatible local files:

- `fcea2`
- `thermo.lib`
- `trans.lib`

Place the files beside the GRIBS Python file:

```text
GRIBS/
├── GRIBS_v0.3.1-alpha.py
├── gribs_config.json
├── fcea2
├── thermo.lib
├── trans.lib
└── results/
```

On Linux, the executable may require permission:

```bash
chmod +x fcea2
```

NASA CEA is not developed, maintained, or distributed by the GRIBS project. Users are responsible for obtaining compatible CEA software and databases from an authorized source and for complying with all applicable license, redistribution, usage, and export-control requirements.

Official NASA CEA project:

- https://github.com/nasa/CEA
- https://nasa.github.io/cea/installation.html

GRIBS uses CEA HP equilibrium calculations to construct pressure-dependent chamber-gas properties. GRIBS retains its own unsteady chamber, burn-rate, grain-regression, nozzle-flow, thrust, and blowdown models.

### Legacy-fit mode

The `legacy_fit` backend does not require `fcea2`, `thermo.lib`, or `trans.lib`.

The following quantities are supplied manually in `thermochemistry.legacy_fit`:

- propellant density
- gas-constant correlation coefficients
- chamber-temperature correlation coefficients
- specific-heat-ratio correlation coefficients

The implemented correlations are:

```text
R(p)     = R_a + R_b * ln(p / 1e5)
T0(p)    = eta_T0 * [T_a + T_b * ln(p / 1e5)]
gamma(p) = g_a + g_b * ln(p / 1e5)
```

The example coefficients are formulation-specific and must not be assumed valid for other propellants.

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

## Outputs

The default output directory is `results/`. Output names are configurable in `gribs_config.json`.

```text
results/
├── ballistics.png
├── time_history.csv
├── summary.txt
├── summary.json
└── cea_cache/
```

- `ballistics.png`: combined publication-style result figure
- `time_history.csv`: sampled time history
- `summary.txt`: human-readable result summary
- `summary.json`: original configuration, resolved runtime inputs, derived values, and results
- `cea_cache/`: reusable pressure-dependent CEA property tables

## Validation status

For the `v0.3.1-alpha` release:

- Python syntax compilation was checked.
- The unified JSON configuration loader was exercised.
- The `legacy_fit` path completed the numerical self-tests, simulation, post-processing, and all four primary output files.
- CEA2 operation was developed against the external legacy `fcea2` workflow with eight assigned pressure points per batch and 81 total pressure-grid points.
- External CEA runtime files are not included, so users must validate CEA compatibility in the target environment.

## Known limitations

- This is alpha software and is not certified for safety-critical use.
- Chemical equilibrium does not model finite-rate chemistry.
- The additive-volume density rule is an approximation.
- Burn-rate and erosive-burning coefficients are empirical inputs.
- CEA compatibility depends on the user's external executable and databases.
- The current CEA integration expects the legacy interactive executable and plot-file workflow.
- Model outputs are not substitutes for static firing tests, material characterization, or independent review.

## License

GRIBS is licensed under the MIT License.

Copyright (c) 2026 Hikaru Kurokawa.

The GRIBS license applies only to files distributed by the GRIBS project.
External NASA CEA software and databases are not included in this repository
and remain subject to their respective terms.

## Citation

Citation metadata is provided in `CITATION.cff`.

## Contributing and support

See `CONTRIBUTING.md` for contribution guidance. Use GitHub Issues for reproducible bug reports and feature proposals. See `SECURITY.md` for private security-reporting guidance.
