"""A-6 unit 3a (report 118, Medium-1): the conftest state guard restores the thermochemistry backend between tests.
tests/_a6_state_guard_cases.py runs in a child pytest (fresh process, fixed order): with the conftest both tests pass;
with --noconftest (control) the second one fails, so the pass is due to the guard."""
import os, re, subprocess, sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
CASES = Path(__file__).with_name("_a6_state_guard_cases.py")


def _child_pytest(*extra):
    env = dict(os.environ)
    env.pop("PYTEST_ADDOPTS", None)
    r = subprocess.run([sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider", "-rA", *extra, str(CASES)],
                       cwd=ROOT, capture_output=True, text=True, timeout=300, env=env)
    outcome = dict((name, status) for status, name in
                   re.findall(r"^(PASSED|FAILED) \S*::(test_\w+)", r.stdout, flags=re.M))
    return r, outcome


@pytest.mark.integration
def test_conftest_guard_restores_the_backend_between_tests():
    r, outcome = _child_pytest()
    assert r.returncode == 0, r.stdout + r.stderr
    assert outcome == {"test_1_installs_a_backend_and_leaves_it": "PASSED",
                       "test_2_sees_the_state_from_before_test_1": "PASSED"}, r.stdout
    r, outcome = _child_pytest("--noconftest")  # control: without the guard the backend of test 1 is still installed
    assert r.returncode == 1, r.stdout + r.stderr
    assert outcome == {"test_1_installs_a_backend_and_leaves_it": "PASSED",
                       "test_2_sees_the_state_from_before_test_1": "FAILED"}, r.stdout
