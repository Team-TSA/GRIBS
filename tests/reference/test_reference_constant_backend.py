"""Phase A-5 acceptance tests for the constant-property reference backend."""

from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
MODULE_FILES = sorted(ROOT.glob("GRIBS_v*.py"))

if len(MODULE_FILES) != 1:
    raise RuntimeError(
        f"Expected exactly one GRIBS_v*.py file, found {len(MODULE_FILES)}: "
        f"{[path.name for path in MODULE_FILES]}"
    )

MODULE_PATH = MODULE_FILES[0]
SPEC = importlib.util.spec_from_file_location("gribs_a5_reference", MODULE_PATH)

if SPEC is None or SPEC.loader is None:
    raise ImportError(f"Could not load GRIBS module from {MODULE_PATH}")

gribs = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = gribs
SPEC.loader.exec_module(gribs)


@pytest.mark.fast
def test_reference_constant_backend_configuration_is_explicit_and_selectable():
    """The A-5 constant-property backend is selectable from JSON configuration."""
    source = ROOT / "gribs_config.json"
    supplied = json.loads(source.read_text(encoding="utf-8"))

    supplied["thermochemistry"]["backend"] = "reference_constant"
    supplied["thermochemistry"].pop("cea_python", None)
    supplied["thermochemistry"].pop("cea_legacy_executable", None)
    supplied["thermochemistry"]["reference_constant"] = {
        "temperature_K": 3000.0,
        "gamma": 1.2,
        "molecular_weight_kg_kmol": 25.0,
    }

    configuration = gribs.configuration_from_document(supplied)
    config = configuration.config

    assert config.thermo_backend == "reference_constant"
    assert config.reference_temperature_K == pytest.approx(3000.0)
    assert config.reference_gamma == pytest.approx(1.2)
    assert config.reference_molecular_weight_kg_kmol == pytest.approx(25.0)


def reference_configuration_document():
    """Return a fresh valid reference_constant configuration document."""
    source = ROOT / "gribs_config.json"
    supplied = json.loads(source.read_text(encoding="utf-8"))
    supplied["thermochemistry"]["backend"] = "reference_constant"
    supplied["thermochemistry"].pop("cea_python", None)
    supplied["thermochemistry"].pop("cea_legacy_executable", None)
    supplied["thermochemistry"]["reference_constant"] = {
        "temperature_K": 3000.0,
        "gamma": 1.2,
        "molecular_weight_kg_kmol": 25.0,
    }
    return supplied


@pytest.mark.fast
def test_reference_constant_backend_requires_its_configuration_section():
    """Selecting reference_constant requires its explicit property block."""
    supplied = reference_configuration_document()
    supplied["thermochemistry"].pop("reference_constant")

    with pytest.raises(
        gribs.ConfigurationError,
        match=r"thermochemistry\.reference_constant.*missing",
    ):
        gribs.configuration_from_document(supplied)


@pytest.mark.fast
def test_reference_constant_backend_rejects_unknown_property_key():
    """Unknown constant-property keys must not be silently ignored."""
    supplied = reference_configuration_document()
    supplied["thermochemistry"]["reference_constant"]["unknown_property"] = 1.0

    with pytest.raises(
        gribs.ConfigurationError,
        match=r"thermochemistry\.reference_constant.*unknown key",
    ):
        gribs.configuration_from_document(supplied)


@pytest.mark.fast
@pytest.mark.parametrize(
    ("key", "invalid_value"),
    [
        ("temperature_K", 0.0),
        ("temperature_K", -1.0),
        ("gamma", 1.0),
        ("gamma", 0.9),
        ("molecular_weight_kg_kmol", 0.0),
        ("molecular_weight_kg_kmol", -1.0),
    ],
)
def test_reference_constant_backend_rejects_nonphysical_properties(
    key,
    invalid_value,
):
    """The selectable reference backend accepts only physical constants."""
    supplied = reference_configuration_document()
    supplied["thermochemistry"]["reference_constant"][key] = invalid_value

    with pytest.raises(
        gribs.ConfigurationError,
        match=rf"thermochemistry\.reference_constant.*{key}",
    ):
        gribs.configuration_from_document(supplied)


@pytest.mark.integration
def test_build_backend_creates_reference_constant_backend(tmp_path):
    """The selectable reference backend must use the normal backend factory."""
    supplied = reference_configuration_document()
    configuration = gribs.configuration_from_document(supplied)
    config = configuration.config

    backend = gribs.build_backend(
        config,
        config.reactants_spec_list,
        config.packing_fraction,
        tmp_path,
    )

    assert backend.name == "reference_constant"
    assert backend.uses_cea is False
    assert backend.supports_two_phase is True
    assert backend.props(config.p0_init) == pytest.approx(
        (
            gribs.R_UNIVERSAL / config.reference_molecular_weight_kg_kmol,
            config.reference_temperature_K * config.eta_T0,
            config.reference_gamma,
        )
    )
    assert backend.phase_fractions(config.p0_init) == pytest.approx((1.0, 0.0))

    psi, dpsi_dp = backend.psi_and_derivative(config.p0_init)
    expected_psi = (
        gribs.R_UNIVERSAL
        / config.reference_molecular_weight_kg_kmol
        * config.reference_temperature_K
        * config.eta_T0
    )
    assert psi == pytest.approx(expected_psi)
    assert dpsi_dp == pytest.approx(0.0, abs=1.0e-12)

    metadata = backend.metadata()
    assert metadata["backend"] == "reference_constant"
    assert metadata["uses_cea"] is False
    assert metadata["constant_properties"] is True


@pytest.mark.integration
def test_initialize_reference_backend_uses_normal_product_path(
    tmp_path,
    monkeypatch,
):
    """Normal initialization registers the backend and updates propellant density."""
    supplied = reference_configuration_document()
    configuration = gribs.configuration_from_document(supplied)
    config = configuration.config

    original_backend = gribs._THERMO
    monkeypatch.setattr(gribs, "_THERMO", None)

    backend = gribs.initialize_thermochemistry(config, tmp_path)

    assert isinstance(backend, gribs.ReferenceConstantBackend)
    assert gribs.thermochemistry() is backend
    assert config.rho_p == pytest.approx(backend.bulk_mixture_density)
    assert backend.bulk_mixture_density == pytest.approx(
        gribs.ideal_mixture_density(config.reactants_spec_list)
        * config.packing_fraction
    )
    assert backend.cache_path is None
    assert backend.cache_used is False
    assert backend.cache_written is False
    assert list(tmp_path.iterdir()) == []

    assert backend.is_extrapolated(config.p_fit_min) is False
    assert backend.is_extrapolated(config.p_fit_max) is False
    assert backend.is_extrapolated(config.p_fit_min * 0.5) is True
    assert backend.is_extrapolated(config.p_fit_max * 2.0) is True

    monkeypatch.setattr(gribs, "_THERMO", original_backend)


@pytest.mark.integration
def test_level0_constant_property_blowdown_matches_analytic_solution(
    tmp_path,
    monkeypatch,
):
    """Level 0 chamber pressure matches the closed-form choked solution."""
    supplied = reference_configuration_document()

    burn_law = supplied["propellant"]["burn_law"]
    burn_law["coefficient_m_s"] = 1.0e-12
    burn_law["pressure_exponent"] = 0.0
    burn_law["temperature_sensitivity_1_K"] = 0.0
    burn_law["erosive_burning"]["enabled"] = False

    supplied["grain"]["burning_end_faces"] = 0
    supplied["nozzle"]["expansion_ratio"] = 1.0
    supplied["nozzle"]["discharge_coefficient"] = 1.0
    supplied["nozzle"]["thrust_efficiency"] = 1.0
    supplied["nozzle"]["flow_separation"]["enabled"] = False
    supplied["nozzle"]["throat_erosion"]["enabled"] = False

    supplied["environment"]["ambient_pressure_Pa"] = 1.0e5
    supplied["environment"]["initial_chamber_pressure_Pa"] = 5.0e6

    supplied["igniter"]["mass_flow_kg_s"] = 0.0
    supplied["igniter"]["duration_s"] = 0.0

    supplied["solver"]["method"] = "LSODA"
    supplied["solver"]["relative_tolerance"] = 1.0e-10
    supplied["solver"]["burning_time_limit_s"] = 5.0
    supplied["solver"]["maximum_burning_step_s"] = 0.02
    supplied["solver"]["unchoked_policy"] = "switch"
    supplied["solver"]["blowdown"]["enabled"] = False

    configuration = gribs.configuration_from_document(supplied)
    config = configuration.config

    monkeypatch.setattr(gribs, "_THERMO", None)
    backend = gribs.initialize_thermochemistry(config, tmp_path)

    result = gribs.run_model(config)
    solution = result["sol_burn"]

    assert solution.success
    assert solution.y.shape[0] == 7
    assert result["sol_blow"] is None
    assert result["init_regime"] == "choked"

    throat_area = gribs.math.pi * config.R_t0**2
    gas_constant, temperature, gamma = backend.props(config.p0_init)
    psi, dpsi_dp = backend.psi_and_derivative(config.p0_init)

    assert dpsi_dp == pytest.approx(0.0, abs=0.0)

    critical_factor = (
        2.0 / (gamma + 1.0)
    ) ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))
    flow_coefficient = (
        config.Cd
        * throat_area
        * gribs.math.sqrt(gamma / (gas_constant * temperature))
        * critical_factor
    )
    time_constant = config.V_g0 / (psi * flow_coefficient)

    initial_burning_area = gribs.geometry(0.0, config)[2]
    generated_mass_flow = (
        config.rho_p
        * initial_burning_area
        * config.a_burn
    )
    asymptotic_pressure = generated_mass_flow / flow_coefficient

    choking_pressure = gribs.choke_limit_pressure(
        config,
        throat_area,
        throat_area,
    )
    assert choking_pressure is not None

    time = solution.t
    numerical_pressure = solution.y[gribs.IP]
    choked_segment = numerical_pressure > 1.05 * choking_pressure

    assert gribs.np.count_nonzero(choked_segment) >= 10

    analytic_pressure = (
        asymptotic_pressure
        + (config.p0_init - asymptotic_pressure)
        * gribs.np.exp(-time[choked_segment] / time_constant)
    )
    maximum_relative_error = float(
        gribs.np.max(
            gribs.np.abs(
                numerical_pressure[choked_segment] / analytic_pressure - 1.0
            )
        )
    )

    assert maximum_relative_error <= 1.0e-8


QUASI_STEADY_MAXIMUM_RELATIVE_ERROR = 1.5e-3


@pytest.mark.integration
@pytest.mark.parametrize("burning_end_faces", [0, 1, 2])
def test_level1_quasi_steady_pressure_tracks_equilibrium(
    burning_end_faces,
    tmp_path,
    monkeypatch,
):
    """Low-exponent numerical pressure tracks quasi-steady equilibrium."""
    supplied = reference_configuration_document()

    burn_law = supplied["propellant"]["burn_law"]
    burn_law["pressure_exponent"] = 0.01
    burn_law["temperature_sensitivity_1_K"] = 0.0
    burn_law["erosive_burning"]["enabled"] = False

    supplied["grain"]["burning_end_faces"] = burning_end_faces

    supplied["nozzle"]["expansion_ratio"] = 1.0
    supplied["nozzle"]["flow_separation"]["enabled"] = False
    supplied["nozzle"]["throat_erosion"]["enabled"] = False

    supplied["igniter"]["mass_flow_kg_s"] = 0.0
    supplied["igniter"]["duration_s"] = 0.0

    supplied["solver"]["method"] = "LSODA"
    supplied["solver"]["relative_tolerance"] = 1.0e-10
    supplied["solver"]["unchoked_policy"] = "switch"
    supplied["solver"]["blowdown"]["enabled"] = False

    configuration = gribs.configuration_from_document(supplied)
    config = configuration.config

    monkeypatch.setattr(gribs, "_THERMO", None)
    gribs.initialize_thermochemistry(config, tmp_path)

    result = gribs.run_model(config)
    solution = result["sol_burn"]

    assert solution.success
    assert solution.y.shape[0] == 7
    assert result["stop_burn"] == "burnout"
    assert result["sol_blow"] is None

    index_start = max(1, solution.t.size // 4)
    index_stop = max(index_start + 1, 3 * solution.t.size // 4)
    indices = gribs.np.unique(
        gribs.np.linspace(
            index_start,
            index_stop - 1,
            min(60, index_stop - index_start),
        ).astype(int)
    )

    relative_errors = []
    equilibrium_pressures = []

    for index in indices:
        pressure = float(solution.y[gribs.IP, index])
        burn_depth = float(solution.y[gribs.IX, index])
        throat_radius = float(solution.y[gribs.IRT, index])

        bore_radius, unused_length, burning_area, unused_volume = (
            gribs.geometry(burn_depth, config)
        )
        throat_area = gribs.math.pi * throat_radius**2

        equilibrium = gribs.equilibrium_pressure(
            burning_area,
            bore_radius,
            throat_area,
            result["Ae"],
            config,
        )
        if equilibrium is None or not gribs.np.isfinite(equilibrium):
            continue

        equilibrium_pressures.append(float(equilibrium))
        relative_errors.append(
            abs(pressure / float(equilibrium) - 1.0)
        )

    assert len(relative_errors) >= 10

    equilibrium_pressures = gribs.np.asarray(
        equilibrium_pressures,
        dtype=float,
    )
    relative_errors = gribs.np.asarray(relative_errors, dtype=float)

    assert gribs.np.all(
        equilibrium_pressures >= config.p_fit_min
    )
    assert gribs.np.all(
        equilibrium_pressures <= config.p_fit_max
    )

    maximum_relative_error = float(gribs.np.max(relative_errors))
    assert (
        maximum_relative_error
        <= QUASI_STEADY_MAXIMUM_RELATIVE_ERROR
    )


REFERENCE_CASES = [
    (ROOT / "cases" / "reference" / "level1_end_faces_0.json", 0),
    (ROOT / "cases" / "reference" / "level1_end_faces_1.json", 1),
    (ROOT / "cases" / "reference" / "level1_end_faces_2.json", 2),
]


@pytest.mark.integration
@pytest.mark.parametrize(
    ("case_path", "expected_end_faces"),
    REFERENCE_CASES,
    ids=["end-faces-0", "end-faces-1", "end-faces-2"],
)
def test_level1_reference_case_file_runs_without_cea(
    case_path,
    expected_end_faces,
    tmp_path,
    monkeypatch,
):
    """Each committed Level 1 case loads and runs through the product path."""
    assert case_path.is_file()

    configuration = gribs.load_json_configuration(case_path)
    config = configuration.config

    assert config.thermo_backend == "reference_constant"
    assert config.n_end == pytest.approx(float(expected_end_faces))
    assert config.n_burn == pytest.approx(0.01)
    assert config.eps_nozzle == pytest.approx(1.0)
    assert config.use_separation is False
    assert config.ero_alpha == pytest.approx(0.0)
    assert config.ero_throat_c == pytest.approx(0.0)
    assert config.ign_mdot == pytest.approx(0.0)
    assert config.ign_time == pytest.approx(0.0)
    assert config.blowdown is False
    assert config.rtol == pytest.approx(1.0e-10)

    monkeypatch.setattr(gribs, "_THERMO", None)
    backend = gribs.initialize_thermochemistry(config, tmp_path)

    assert isinstance(backend, gribs.ReferenceConstantBackend)
    assert backend.uses_cea is False
    assert list(tmp_path.iterdir()) == []

    result = gribs.run_model(config)
    solution = result["sol_burn"]

    assert solution.success
    assert solution.y.shape[0] == 7
    assert result["stop_burn"] == "burnout"
    assert result["sol_blow"] is None

    midpoint = solution.t.size // 2
    pressure = float(solution.y[gribs.IP, midpoint])
    burn_depth = float(solution.y[gribs.IX, midpoint])
    throat_radius = float(solution.y[gribs.IRT, midpoint])

    bore_radius, unused_length, burning_area, unused_volume = (
        gribs.geometry(burn_depth, config)
    )
    throat_area = gribs.math.pi * throat_radius**2
    equilibrium = gribs.equilibrium_pressure(
        burning_area,
        bore_radius,
        throat_area,
        result["Ae"],
        config,
    )

    assert equilibrium is not None
    assert gribs.np.isfinite(equilibrium)
    assert config.p_fit_min <= equilibrium <= config.p_fit_max
    assert (
        abs(pressure / float(equilibrium) - 1.0)
        <= QUASI_STEADY_MAXIMUM_RELATIVE_ERROR
    )


@pytest.mark.integration
def test_reference_constant_cli_dumps_deterministic_properties(tmp_path):
    """--dump-thermo supports the selectable constant-property backend."""
    case_path = (
        ROOT
        / "cases"
        / "reference"
        / "level1_end_faces_0.json"
    )
    dump_path = tmp_path / "reference_constant_properties.json"

    result = subprocess.run(
        [
            sys.executable,
            str(MODULE_PATH),
            "--config",
            str(case_path),
            "--dump-thermo",
            str(dump_path),
        ],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, (
        f"stdout:\n{result.stdout}\n"
        f"stderr:\n{result.stderr}"
    )
    assert dump_path.is_file()

    payload = json.loads(
        dump_path.read_text(encoding="utf-8")
    )

    assert payload["backend"] == "reference_constant"
    assert payload["provenance"]["constant_properties"] is True
    assert payload["pressure_Pa"] == pytest.approx(
        [1.0e5, 8.0e6]
    )
    assert payload["temperature_K"] == pytest.approx(
        [3000.0, 3000.0]
    )
    assert payload["gamma_s"] == pytest.approx(
        [1.2, 1.2]
    )
    assert payload[
        "gas_phase_molecular_weight_kg_kmol"
    ] == pytest.approx(
        [25.0, 25.0]
    )
    assert payload["gas_mass_fraction"] == pytest.approx(
        [1.0, 1.0]
    )
    assert payload["condensed_mass_fraction"] == pytest.approx(
        [0.0, 0.0]
    )

    expected_gas_constant = gribs.R_UNIVERSAL / 25.0
    expected_psi = expected_gas_constant * 3000.0

    assert payload["gas_constant_J_kgK"] == pytest.approx(
        [expected_gas_constant, expected_gas_constant]
    )
    assert payload["psi_J_kg"] == pytest.approx(
        [expected_psi, expected_psi]
    )


@pytest.mark.integration
def test_reference_constant_banner_reports_analytic_properties(
    tmp_path,
):
    """The CLI banner describes constants, not a nonexistent property table."""
    case_path = REFERENCE_CASES[0][0]
    dump_path = tmp_path / "properties.json"

    result = subprocess.run(
        [
            sys.executable,
            str(MODULE_PATH),
            "--config",
            str(case_path),
            "--dump-thermo",
            str(dump_path),
        ],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, (
        f"stdout:\n{result.stdout}\n"
        f"stderr:\n{result.stderr}"
    )
    assert "thermochemistry backend : reference_constant" in result.stdout
    assert "constant analytic properties" in result.stdout
    assert "T0 = 3000 K" in result.stdout
    assert "gamma = 1.2" in result.stdout
    assert "M = 25 kg/kmol" in result.stdout
    assert "CEA used              : no" in result.stdout
    assert "property cache        : disabled" in result.stdout
    assert "81 points" not in result.stdout


@pytest.mark.fast
def test_reference_constant_parameters_are_documented():
    """Programmatic parameter documentation includes the A-5 backend."""
    backend_documentation = gribs.PARAM_DOC["thermo_backend"]

    assert "reference_constant" in backend_documentation
    assert (
        gribs.PARAM_DOC["reference_temperature_K"]
        == "reference_constant backend temperature [K]"
    )
    assert (
        gribs.PARAM_DOC["reference_gamma"]
        == "reference_constant backend specific-heat ratio [-]"
    )
    assert (
        gribs.PARAM_DOC["reference_molecular_weight_kg_kmol"]
        == "reference_constant backend molecular weight [kg/kmol]"
    )



@pytest.mark.integration
def test_cli_backend_override_uses_documented_reference_defaults(tmp_path):
    """CLI override selects reference constants without invoking CEA."""
    dump_path = tmp_path / "override_reference_properties.json"

    result = subprocess.run(
        [
            sys.executable,
            str(MODULE_PATH),
            "--config",
            str(ROOT / "gribs_config.json"),
            "--backend",
            "reference_constant",
            "--dump-thermo",
            str(dump_path),
        ],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, (
        f"stdout:\n{result.stdout}\n"
        f"stderr:\n{result.stderr}"
    )
    assert dump_path.is_file()

    expected_warning = (
        "thermochemistry.backend overridden on the command line "
        "to 'reference_constant'"
    )
    assert expected_warning in result.stderr
    assert "thermochemistry backend : reference_constant" in result.stdout
    assert "constant analytic properties" in result.stdout
    assert "T0 = 3000 K" in result.stdout
    assert "gamma = 1.2" in result.stdout
    assert "M = 25 kg/kmol" in result.stdout
    assert "official CEA package" not in result.stdout
    assert "CEA library version" not in result.stdout

    payload = json.loads(dump_path.read_text(encoding="utf-8"))

    assert payload["backend"] == "reference_constant"
    assert payload["provenance"]["uses_cea"] is False
    assert payload["provenance"]["constant_properties"] is True
    assert payload["temperature_K"] == pytest.approx([3000.0, 3000.0])
    assert payload["gamma_s"] == pytest.approx([1.2, 1.2])

    molecular_weight = payload[
        "gas_phase_molecular_weight_kg_kmol"
    ]
    assert molecular_weight == pytest.approx([25.0, 25.0])
