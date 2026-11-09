"""A-6 unit 3b (report 120, condition 3 of U-48 / U-49): CEALegacyExecutableBackend._resolve_executable looks for a
relative fcea2 path first in the program directory (the directory of the entry, gribs.thermo._PROGRAM_DIR), then in
the working directory, then on PATH. Before the move the first place was Path(__file__).resolve().parent of the entry;
these tests check that it is still the entry's directory, without patching anything."""
import importlib
import os
import shutil
import subprocess
import sys
import types
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
NAME = "a6_fcea2_resolve_probe"  # no such file next to the entry, in the working directory or on PATH


def _entry_dir():
    bodies = sorted(ROOT.glob("GRIBS_v*.py"))
    if len(bodies) != 1:
        raise RuntimeError(f"expected exactly one GRIBS_v*.py, found {[b.name for b in bodies]}")
    return bodies[0].resolve().parent


def _probe(thermo):
    backend = thermo.CEALegacyExecutableBackend.__new__(thermo.CEALegacyExecutableBackend)
    backend.c = types.SimpleNamespace(cea_legacy_executable_path=NAME)
    return backend


@pytest.mark.fast
def test_relative_executable_is_searched_next_to_the_entry_first(tmp_path, monkeypatch):
    thermo = importlib.import_module("gribs.thermo")
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("PATH", str(tmp_path / "empty"))
    assert not (_entry_dir() / NAME).exists()
    with pytest.raises(thermo.BackendUnavailableError) as err:
        _probe(thermo)._resolve_executable()
    searched = str(err.value).split("Searched: ", 1)[1].split(", and the system PATH.", 1)[0]
    assert searched == f"{_entry_dir() / NAME}, {tmp_path / NAME}", searched


@pytest.mark.integration
def test_relative_executable_found_in_the_program_directory_of_a_copied_package(tmp_path):
    """A copy <tmp>/prog/src/gribs: an executable <tmp>/prog/NAME is found, before one with the same name in the working
    directory (fresh process, so this copy of the package is the one imported)."""
    prog, cwd = tmp_path / "prog", tmp_path / "cwd"
    shutil.copytree(ROOT / "src" / "gribs", prog / "src" / "gribs", ignore=shutil.ignore_patterns("__pycache__"))
    cwd.mkdir()
    for d in (prog, cwd):
        (d / NAME).write_text("#!/bin/sh\nexit 0\n")
        os.chmod(d / NAME, 0o755)
    child = ("import sys, types\n"
             "sys.path.insert(0, sys.argv[1])\n"
             "import gribs.thermo as t\n"
             "b = t.CEALegacyExecutableBackend.__new__(t.CEALegacyExecutableBackend)\n"
             "b.c = types.SimpleNamespace(cea_legacy_executable_path=sys.argv[2])\n"
             "print(t.__file__)\n"
             "print(b._resolve_executable())\n")
    env = {k: v for k, v in os.environ.items() if k not in ("PYTHONPATH", "PYTEST_ADDOPTS")}
    r = subprocess.run([sys.executable, "-c", child, str(prog / "src"), NAME], cwd=cwd, capture_output=True,
                       text=True, timeout=120, env=env)
    assert r.returncode == 0, r.stdout + r.stderr
    module_file, found = r.stdout.splitlines()
    assert Path(module_file).resolve() == (prog / "src" / "gribs" / "thermo.py").resolve()
    assert Path(found) == (prog / NAME).resolve()
