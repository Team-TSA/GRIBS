"""P2 (pending approval): CEA backend + single_phase_legacy must not recurse.

Regression introduced by A-6 commit 3c24fbc (theta/psi unification).  The
fixed path must return exactly the pre-A-6 value PCHIP(z, R*T).
"""
import copy, importlib.util, json, math, sys
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = sorted(ROOT.glob("GRIBS_v*.py"))[0]
if "gribs" in sys.modules:
    gribs = sys.modules["gribs"]
else:
    spec = importlib.util.spec_from_file_location("gribs", MODULE_PATH)
    gribs = importlib.util.module_from_spec(spec); sys.modules["gribs"] = gribs
    spec.loader.exec_module(gribs)


def _doc(tmp_path):
    doc = json.loads((ROOT / "cases" / "golden" / "g2_single_phase_cea.json").read_text())
    doc = copy.deepcopy(doc)
    doc["output"]["directory"] = str(tmp_path / "out")
    doc["two_phase_model"]["mode"] = "single_phase_legacy"
    doc["thermochemistry"]["cea_python"]["cache_directory"] = str(tmp_path / "cache")
    return doc


@pytest.mark.integration
def test_cea_backend_single_phase_legacy_psi_has_no_recursion(tmp_path):
    pytest.importorskip("cea")
    conf = gribs.configuration_from_document(_doc(tmp_path)).config
    b = gribs.initialize_thermochemistry(conf, tmp_path / "cache")
    p = 3.0e6
    psi, dpsi = b.psi_and_derivative(p)
    theta, dtheta = b.theta_and_derivative(p)
    z = math.log(p)
    assert psi == theta == float(b.psii(z))      # pre-A-6 thi was PCHIP(z, R*T)
    assert dpsi == dtheta


@pytest.mark.integration
def test_cea_backend_single_phase_legacy_run_completes(tmp_path):
    pytest.importorskip("cea")
    cfg = tmp_path / "c.json"; cfg.write_text(json.dumps(_doc(tmp_path)))
    assert gribs.main(["--config", str(cfg)]) == 0



class _NoTwoPhaseSynthetic(gribs._SyntheticTestBackend):
    """Synthetic table with supports_two_phase=False: the same branch that
    CEALegacyExecutableBackend (fcea2) takes in every mode, without fcea2."""
    supports_two_phase = False


def _ref_config(mode):
    doc = json.loads((ROOT / "cases" / "reference" / "level1_end_faces_0.json").read_text())
    doc["two_phase_model"] = dict(doc.get("two_phase_model", {}), mode=mode)
    return gribs.configuration_from_document(doc).config


@pytest.mark.fast
@pytest.mark.parametrize("mode", [gribs.TWO_PHASE_SINGLE_LEGACY])
def test_table_backend_without_two_phase_support_has_no_recursion(mode):
    # homogeneous + supports_two_phase=False is rejected by configuration
    # validation (ConfigurationError), so only legacy reaches this branch.
    b = _NoTwoPhaseSynthetic(_ref_config(mode))
    assert b.c.two_phase_mode == mode
    p = 1.0e6
    psi, dpsi = b.psi_and_derivative(p)
    theta, dtheta = b.theta_and_derivative(p)
    assert psi == theta == float(b.psii(math.log(p)))
    assert dpsi == dtheta
    assert b.phase_fractions(p) == (1.0, 0.0)


@pytest.mark.fast
def test_synthetic_backend_built_in_legacy_mode_has_no_recursion():
    """Self-test T15 keeps the backend's own mode homogeneous, so it never
    reaches this branch; here the backend itself is built in legacy mode."""
    b = gribs._SyntheticTestBackend(_ref_config(gribs.TWO_PHASE_SINGLE_LEGACY))
    assert b.c.two_phase_mode == gribs.TWO_PHASE_SINGLE_LEGACY
    psi, _ = b.psi_and_derivative(1.0e6)
    assert psi == float(b.psii(math.log(1.0e6)))
