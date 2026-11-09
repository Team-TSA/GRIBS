"""A-6 (report 118 Medium-1, report 119): the conftest state guard restores the thermochemistry backend between tests.
Each helper runs in a child pytest (fresh process, fixed order): with the conftest both of its tests pass; with
--noconftest (control) the second one fails, so the pass is due to the guard.
- tests/_a6_state_guard_cases.py: gribs.thermo is loaded at collection; the guard restores the saved backend.
- tests/_a6_state_guard_late_import_cases.py: gribs.thermo is first imported inside a test; the guard resets it to None."""
import os, re, subprocess, sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
HELPERS = {
    "_a6_state_guard_cases.py": ("test_1_installs_a_backend_and_leaves_it", "test_2_sees_the_state_from_before_test_1"),
    "_a6_state_guard_late_import_cases.py": ("test_1_loads_the_entry_and_installs_a_backend",
                                             "test_2_sees_the_state_of_a_fresh_import"),
}


def _child_pytest(helper, *extra):
    env = dict(os.environ)
    env.pop("PYTEST_ADDOPTS", None)
    r = subprocess.run([sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider", "-rA", *extra,
                        str(Path(__file__).with_name(helper))],
                       cwd=ROOT, capture_output=True, text=True, timeout=300, env=env)
    outcome = dict((name, status) for status, name in
                   re.findall(r"^(PASSED|FAILED) \S*::(test_\w+)", r.stdout, flags=re.M))
    return r, outcome


@pytest.mark.integration
@pytest.mark.parametrize("helper", sorted(HELPERS))
def test_conftest_guard_restores_the_backend_between_tests(helper):
    first, second = HELPERS[helper]
    r, outcome = _child_pytest(helper)
    assert r.returncode == 0, r.stdout + r.stderr
    assert outcome == {first: "PASSED", second: "PASSED"}, r.stdout
    r, outcome = _child_pytest(helper, "--noconftest")  # control: without the guard the backend of test 1 is still installed
    assert r.returncode == 1, r.stdout + r.stderr
    assert outcome == {first: "PASSED", second: "FAILED"}, r.stdout
