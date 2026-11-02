"""P4 v2 (proposal): self-test comparisons use explicitly-moded backends,
independent of the start configuration's two_phase_mode; _ThermoSwap
restores the global backend on normal exit and on exceptions."""
import importlib.util, json, sys, tempfile
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


@pytest.mark.fast
def test_thermoswap_restores_on_normal_exit_and_exception():
    sentinel_outer, inner = object(), object()
    saved = gribs._THERMO
    try:
        gribs._THERMO = sentinel_outer
        with gribs._ThermoSwap(inner):
            assert gribs._THERMO is inner
        assert gribs._THERMO is sentinel_outer
        with pytest.raises(RuntimeError):
            with gribs._ThermoSwap(inner):
                assert gribs._THERMO is inner
                raise RuntimeError("boom")
        assert gribs._THERMO is sentinel_outer
        with gribs._ThermoSwap(inner):             # nested
            with gribs._ThermoSwap(sentinel_outer):
                assert gribs._THERMO is sentinel_outer
            assert gribs._THERMO is inner
        assert gribs._THERMO is sentinel_outer
    finally:
        gribs._THERMO = saved


@pytest.mark.integration
@pytest.mark.parametrize("start_mode", [gribs.TWO_PHASE_HOMOGENEOUS, gribs.TWO_PHASE_SINGLE_LEGACY])
def test_selftest_backend_modes_independent_of_start_mode(start_mode, monkeypatch, tmp_path):
    pytest.importorskip("cea")
    doc = json.loads((ROOT / "gribs_config.json").read_text())
    doc["two_phase_model"]["mode"] = start_mode
    doc["thermochemistry"]["cea_python"]["cache_directory"] = str(tmp_path / "cache")
    c = gribs.configuration_from_document(doc).config
    gribs.initialize_thermochemistry(c, tmp_path / "cache")
    outer = gribs._THERMO
    swaps, splits = set(), set()
    enter = gribs._ThermoSwap.__enter__
    def rec_enter(self):
        b = self._replacement
        swaps.add((type(b).__name__, b.c.two_phase_mode, getattr(b, "_ycond", None)))
        return enter(self)
    split = gribs.NozzleTwoPhaseModel.split
    def rec_split(self, ref, yg):
        splits.add((self.mode, round(yg, 9)))
        return split(self, ref, yg)
    monkeypatch.setattr(gribs._ThermoSwap, "__enter__", rec_enter)
    monkeypatch.setattr(gribs.NozzleTwoPhaseModel, "split", rec_split)
    assert gribs.self_tests(c, verbose=False)          # includes the name contract
    assert gribs._THERMO is outer
    assert swaps == {
        ("_SyntheticTestBackend", gribs.TWO_PHASE_HOMOGENEOUS, 0.0),
        ("_SyntheticTestBackend", gribs.TWO_PHASE_HOMOGENEOUS, 0.3),
        ("_SyntheticTestBackend", gribs.TWO_PHASE_SINGLE_LEGACY, 0.0),
    }
    # T14 really takes the condensed homogeneous path (Yg = 1 - 0.3)
    assert (gribs.TWO_PHASE_HOMOGENEOUS, 0.7) in splits
