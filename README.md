# GRIBS

**General Rocket Internal Ballistics Solver**

GRIBS is a Python-based, lumped-parameter internal-ballistics simulator for educational and research use. It integrates grain regression, chamber pressure, pressure-dependent gas properties, nozzle flow, thrust, optional throat erosion, and post-burnout blowdown as an initial-value problem.

> [!CAUTION]
> GRIBS is research software. It is not certified for flight hardware, pressure-vessel design, motor qualification, ignition safety, or operational safety assessment. Numerical results depend on model assumptions and input validity. Independently verify all results and comply with applicable laws, regulations, institutional rules, and safety procedures.

## Status

- Development status: **Research preview / alpha**
- Experimental validation: **Not established by this repository**
- Units: **SI unless explicitly stated otherwise**
- Interface stability: Command-line options and output schemas may change before v1.0

## Features

- Cylindrical-bore grain with 0, 1, or 2 burning end faces
- Saint-Robert pressure-law regression with temperature sensitivity
- Optional Lenoir-Robert-type erosive-burning correction
- Pressure-dependent fitted gas properties
- Choked, subsonic, and simplified separated nozzle-flow regimes
- Converging and converging-diverging nozzle geometries
- Optional throat erosion and igniter gas injection
- LSODA, BDF, and Radau integration through SciPy `solve_ivp`
- Event detection for burnout, choking transitions, and ambient-pressure termination
- Built-in numerical self-tests and mass-consistency diagnostics
- PNG, CSV, text, and JSON outputs

## Model scope

GRIBS uses a spatially uniform chamber model. It does not resolve pressure waves, combustion instability, multidimensional flow, local thermal gradients, structural response, material failure, or asymmetric nozzle loads. Gas properties are fitted as pressure-only functions, and the blowdown model does not solve an independent energy equation.

The default values are examples, not validated design recommendations.

## Requirements

- Python 3.10 or newer
- NumPy
- SciPy
- Matplotlib

## Installation

Clone the repository and create an isolated virtual environment:

```bash
git clone https://github.com/Team-TSA/GRIBS.git
cd GRIBS
python -m venv .venv
```

Activate it, then install the dependencies:

```bash
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

## Quick start

Run the built-in checks first:

```bash
python GRIBS.py --selftest
```

List all adjustable parameters:

```bash
python GRIBS.py --list_params
```

Run the default example:

```bash
python GRIBS.py
```

Run with the included JSON configuration:

```bash
python GRIBS.py --config examples/example_config.json --outdir results/example
```

Configuration precedence is: built-in defaults, JSON settings, then command-line arguments.

## Outputs

A normal run writes:

- `ballistics.png`: combined diagnostic figure
- `time_history.csv`: sampled time history with units in the header
- `summary.txt`: human-readable configuration and performance summary
- `summary.json`: machine-readable inputs and results

Generated results are ignored by Git by default. Preserve the input JSON and `summary.json` together when reproducibility matters.

## Verification before interpreting results

1. Confirm that `--selftest` passes.
2. Check the actual inputs recorded in `summary.json`.
3. Inspect the termination reason and nozzle regimes.
4. Check property-range warnings and mass-consistency errors.
5. Repeat the case with tighter tolerances and, where practical, another solver.
6. Compare against independent analytical, numerical, or experimental evidence.

A passed self-test demonstrates selected implementation consistencies. It does not validate the physical model or supplied properties.

## Documentation

The detailed technical manual is maintained internally and is not included in the public repository. Public model assumptions and limitations are summarized here and in [`MODEL_LIMITATIONS.md`](MODEL_LIMITATIONS.md).

## Contributing and reporting problems

See [`CONTRIBUTING.md`](CONTRIBUTING.md). For defects that could materially affect safety-related interpretation, do not open a public issue. Follow [`SECURITY.md`](SECURITY.md).

## Citation

Citation metadata is provided in [`CITATION.cff`](CITATION.cff). Replace the placeholder repository URL and release metadata before the first public release.

## License

Copyright (c) 2026 Hikaru Kurokawa.

This project is released under the MIT License. See [`LICENSE`](LICENSE). The license disclaimer does not replace the safety notice above.
