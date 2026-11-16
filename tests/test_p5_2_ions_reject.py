"""P5-2 (F-19): thermochemistry.cea_python.ions = true is rejected as a configuration error before any CEA call.

The real ``cea`` package is never reached: ``sys.modules["cea"]`` is replaced by a tripwire module whose every
attribute access fails the test, so a regression shows as a test failure and not as a crash of the test process.
"""
import dataclasses
import json
import sys
import types
from pathlib import Path

import pytest

import gribs.api
import gribs.cli
import gribs.config
import gribs.config.loader
import gribs.thermo
from gribs.errors import ConfigurationError

ROOT = Path(__file__).resolve().parents[1]
G2 = ROOT / "cases" / "golden" / "g2_single_phase_cea.json"


class _Tripwire(types.ModuleType):
    def __getattr__(self, name):
        raise AssertionError(f"cea.{name} was reached; ions = true must be rejected before any CEA call")


@pytest.fixture
def no_cea(monkeypatch):
    monkeypatch.setitem(sys.modules, "cea", _Tripwire("cea"))


def _case(tmp_path, *, ions):
    doc = json.loads(G2.read_text(encoding="utf-8"))
    doc["thermochemistry"]["cea_python"]["ions"] = ions
    doc["output"]["directory"] = str(tmp_path / "out")
    doc["thermochemistry"]["cea_python"]["cache_directory"] = str(tmp_path / "cache")
    path = tmp_path / "case.json"
    path.write_text(json.dumps(doc), encoding="utf-8")
    return path


@pytest.mark.fast
def test_validate_rejects_ions_for_cea_python(tmp_path):
    cfg = gribs.config.loader.load_json_configuration(_case(tmp_path, ions=True)).config
    with pytest.raises(ConfigurationError, match=r"\[E_IONS_NOT_SUPPORTED\] thermochemistry\.cea_python\.ions"):
        gribs.config.validate(cfg)


@pytest.mark.fast
@pytest.mark.parametrize("ions, backend", [(False, None), (True, "reference_constant")])
def test_validate_accepts_ions_false_and_other_backends(tmp_path, ions, backend):
    cfg = gribs.config.loader.load_json_configuration(_case(tmp_path, ions=ions)).config
    if backend is not None:
        cfg = dataclasses.replace(cfg, thermo_backend=backend)    # other backends do not read the option
    gribs.config.validate(cfg)


@pytest.mark.fast
@pytest.mark.parametrize("extra", [["--validate-config"], []])
def test_cli_exits_2_before_cea(tmp_path, capsys, no_cea, extra):
    rc = gribs.cli.main(["--config", str(_case(tmp_path, ions=True))] + extra)
    out, err = capsys.readouterr()
    assert rc == 2
    assert "E_IONS_NOT_SUPPORTED" in err and "is valid" not in out
    assert not (tmp_path / "cache").exists()


@pytest.mark.fast
def test_cli_validate_config_still_accepts_ions_false(tmp_path, capsys, no_cea):
    assert gribs.cli.main(["--config", str(_case(tmp_path, ions=False)), "--validate-config"]) == 0
    assert "is valid" in capsys.readouterr().out


@pytest.mark.fast
def test_python_api_raises_before_cea(tmp_path, no_cea):
    with pytest.raises(ConfigurationError, match="E_IONS_NOT_SUPPORTED"):
        gribs.api.run_simulation(_case(tmp_path, ions=True))


@pytest.mark.fast
def test_backend_guard_rejects_a_config_that_skipped_validate(tmp_path):
    cfg = gribs.config.loader.load_json_configuration(_case(tmp_path, ions=True)).config
    with pytest.raises(ConfigurationError, match="E_IONS_NOT_SUPPORTED"):
        gribs.thermo.CEAPythonBackend._check_option_consistency(types.SimpleNamespace(c=cfg))
