# Contributing

Thank you for helping improve GRIBS.

## Before opening an issue

- Run `python GRIBS.py --selftest`.
- Reproduce the problem with the latest revision.
- Check the effective inputs in `summary.json`.
- Remove personal, confidential, proprietary, and hazardous operational information.
- Use the private process in `SECURITY.md` for potentially safety-significant defects.

## Development setup

```bash
python -m venv .venv
python -m pip install -r requirements-dev.txt
python GRIBS.py --selftest
pytest -q
ruff check GRIBS.py tests
```

## Pull requests

1. Create a focused branch from `main`.
2. Keep physical-model changes separate from formatting-only changes.
3. Add or update tests for modified behavior.
4. Document assumptions, units, references, compatibility effects, and validation evidence.
5. Do not silently alter defaults, equations, safety notices, or output definitions.
6. Confirm that self-tests and automated tests pass.

By contributing, you agree that your contribution is licensed under the repository's MIT License.
