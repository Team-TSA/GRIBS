"""P5-1: the cea_python backend pins the database of the official package.

Fast tests with a stand-in ``cea`` module: the real package is not imported and no database is loaded.
"""
import hashlib
import json
from types import SimpleNamespace

import pytest

import gribs.simulation.postprocess as postprocess
import gribs.thermo as thermo
from gribs.errors import CEADatabaseError, ThermochemistryError


class _FakeCEA:
    """Records init calls; raises ValueError like cea 3.3.4 when a different file is already loaded."""

    def __init__(self, root, loaded_thermo=None, loaded_trans=None):
        self.__file__ = str(root / "__init__.py")
        self._thermo, self._trans = loaded_thermo, loaded_trans
        self.calls = []

    def is_initialized(self):
        return self._thermo is not None

    def init_thermo(self, path):
        self.calls.append(("init_thermo", path))
        if self._thermo not in (None, path):
            raise ValueError("Invalid thermofile name.")
        self._thermo = path

    def init_trans(self, path):
        self.calls.append(("init_trans", path))
        if self._trans not in (None, path):
            raise ValueError("Invalid transfile name.")
        self._trans = path


@pytest.fixture
def package(tmp_path):
    root = tmp_path / "site-packages" / "cea"
    (root / "data").mkdir(parents=True)
    (root / "data" / "thermo.lib").write_bytes(b"bundled thermo")
    (root / "data" / "trans.lib").write_bytes(b"bundled trans")
    return root


@pytest.fixture
def clean_cwd(tmp_path, monkeypatch):
    cwd = tmp_path / "cwd"
    cwd.mkdir()
    monkeypatch.chdir(cwd)
    monkeypatch.delenv("CEA_DATA_DIR", raising=False)
    return cwd


@pytest.mark.fast
def test_pins_bundled_thermo_only_and_records_the_full_digest(package, clean_cwd):
    cea = _FakeCEA(package)
    identity = thermo.pin_bundled_cea_database(cea, transport=False)
    assert cea.calls == [("init_thermo", str(package / "data" / "thermo.lib"))]
    assert identity["thermo_lib_sha256"] == hashlib.sha256(b"bundled thermo").hexdigest()
    assert len(identity["thermo_lib_sha256"]) == 64
    assert identity["trans_lib_loaded_by_gribs"] is False and identity["trans_lib_sha256"] is None
    assert identity["process_already_initialized"] is False
    assert identity["ignored_external_files"] == []


@pytest.mark.fast
def test_transport_loads_the_bundled_trans_lib(package, clean_cwd):
    cea = _FakeCEA(package)
    identity = thermo.pin_bundled_cea_database(cea, transport=True)
    assert ("init_trans", str(package / "data" / "trans.lib")) in cea.calls
    assert identity["trans_lib_sha256"] == hashlib.sha256(b"bundled trans").hexdigest()


@pytest.mark.fast
def test_external_files_are_listed_not_used(package, clean_cwd, tmp_path, monkeypatch, capsys):
    (clean_cwd / "thermo.lib").write_bytes(b"other thermo")
    env_dir = tmp_path / "env"
    env_dir.mkdir()
    (env_dir / "trans.lib").write_bytes(b"damaged")
    monkeypatch.setenv("CEA_DATA_DIR", str(env_dir))
    cea = _FakeCEA(package)
    identity = thermo.pin_bundled_cea_database(cea, transport=False)
    assert cea.calls == [("init_thermo", str(package / "data" / "thermo.lib"))]
    listed = {(f["location"], f["name"]) for f in identity["ignored_external_files"]}
    assert listed == {("cwd", "thermo.lib"), ("CEA_DATA_DIR", "trans.lib")}
    assert thermo.W_EXTERNAL_CEA_DATABASE_IGNORED in capsys.readouterr().err


@pytest.mark.fast
def test_cea_data_dir_pointing_at_the_package_data_is_not_external(package, clean_cwd, monkeypatch):
    monkeypatch.setenv("CEA_DATA_DIR", str(package / "data"))
    identity = thermo.pin_bundled_cea_database(_FakeCEA(package), transport=False)
    assert identity["ignored_external_files"] == []


@pytest.mark.fast
def test_process_initialized_with_the_package_file_is_accepted(package, clean_cwd):
    cea = _FakeCEA(package, loaded_thermo=str(package / "data" / "thermo.lib"))
    identity = thermo.pin_bundled_cea_database(cea, transport=False)
    assert identity["process_already_initialized"] is True


@pytest.mark.fast
def test_process_initialized_with_another_database_is_refused(package, clean_cwd):
    cea = _FakeCEA(package, loaded_thermo="/elsewhere/thermo.lib")
    with pytest.raises(CEADatabaseError) as info:
        thermo.pin_bundled_cea_database(cea, transport=False)
    assert isinstance(info.value, ThermochemistryError)
    assert info.value.code == thermo.E_CEA_DATABASE_ALREADY_INITIALIZED
    assert info.value.context["database_file"] == "thermo.lib"
    assert cea._thermo == "/elsewhere/thermo.lib"


@pytest.mark.fast
def test_cache_key_options_hold_digests_but_no_paths(package, clean_cwd):
    backend = thermo.CEAPythonBackend.__new__(thermo.CEAPythonBackend)
    backend.c = SimpleNamespace(
        cea_py_ions=False, cea_py_transport=False, cea_py_trace=0.0, cea_py_products_from_reactants=True,
        cea_py_product_species=[], cea_py_omit_species=[], cea_py_insert_species=[],
        cea_py_molecular_weight="gas_phase_M", cea_py_smooth_truncation=False, cea_py_truncation_width=0.0)
    backend.cea = SimpleNamespace(__version__="3.3.4", __file__=str(package / "__init__.py"), lib_version=lambda: "3.3.4")
    backend.cea_database = thermo.pin_bundled_cea_database(_FakeCEA(package), transport=False)
    blob = json.dumps(backend.option_signature(), sort_keys=True)
    assert str(package) not in blob and str(clean_cwd) not in blob and '"/' not in blob
    assert backend.cea_database["thermo_lib_sha256"] in blob


@pytest.mark.fast
def test_summary_warning_lists_the_ignored_files():
    ignored = [{"location": "cwd", "name": "thermo.lib", "path": "/x/thermo.lib"}]
    backend = SimpleNamespace(cea_database={"ignored_external_files": ignored, "thermo_lib_path": "/p/thermo.lib",
                                            "thermo_lib_sha256": "0" * 64})
    summary = {"warnings": []}
    postprocess.append_external_cea_database_warning(summary, backend)
    postprocess.append_external_cea_database_warning(summary, backend)
    assert [w["code"] for w in summary["warnings"]] == ["W_EXTERNAL_CEA_DATABASE_IGNORED"]
    assert summary["warnings"][0]["context"]["ignored_files"] == ignored
    quiet = {"warnings": []}
    postprocess.append_external_cea_database_warning(
        quiet, SimpleNamespace(cea_database={"ignored_external_files": []}))
    assert quiet["warnings"] == []
