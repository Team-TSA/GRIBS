# Contributing to GRIBS

Thank you for considering a contribution to GRIBS.

## Before opening a change

- Search existing Issues and Pull Requests.
- Use a focused branch and keep each change limited in scope.
- Do not commit external NASA CEA executables or databases.
- Do not commit generated `results/` files or CEA caches unless a maintainer explicitly requests a specific test fixture.

## Development setup

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
```

For `legacy_fit` development, edit a copy of `gribs_config.json` and set:

```json
{
  "thermochemistry": {
    "backend": "legacy_fit"
  }
}
```

For `cea2` development, separately provide compatible `fcea2`, `thermo.lib`, and `trans.lib` files beside the Python script.

## Required checks

Before submitting a Pull Request:

```bash
python3 -m py_compile GRIBS_v0.3.1-alpha.py
python3 GRIBS_v0.3.1-alpha.py --selftest
python3 GRIBS_v0.3.1-alpha.py
```

Confirm that the expected output files are generated and review `summary.json` for the resolved configuration and verification metrics.

## Pull Requests

A Pull Request should include:

- a concise description of the problem
- the reason for the proposed change
- a list of affected models or configuration fields
- validation steps and results
- any change to assumptions, units, numerical behavior, or compatibility
- documentation updates when behavior changes

Changes to physical equations, numerical algorithms, thermochemistry interpretation, or output definitions must be clearly identified. Avoid combining model changes with unrelated formatting changes.

## Bug reports

Include:

- GRIBS version and Git commit
- operating system
- Python version
- NumPy, SciPy, and Matplotlib versions
- selected thermochemistry backend
- sanitized configuration needed to reproduce the issue
- complete traceback or terminal output
- whether Visual Studio Code Run or a shell command was used

Do not post private credentials, proprietary formulations, export-controlled information, or third-party files that cannot be redistributed.

## Coding style

- Preserve SI units in internal calculations.
- Use explicit unit suffixes in new JSON field names.
- Validate new inputs before beginning production calculations.
- Keep user-facing output in English.
- Preserve deterministic cache keys when adding thermochemistry inputs.
- Add or update self-tests for numerical changes.
