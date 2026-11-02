import importlib.util
import json
import re
import subprocess
import sys
from dataclasses import asdict
from types import SimpleNamespace
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
MODULE_FILES = sorted(ROOT.glob("GRIBS_v*.py"))

if len(MODULE_FILES) != 1:
    raise RuntimeError(
        f"Expected exactly one GRIBS_v*.py file, found {len(MODULE_FILES)}: "
        f"{[path.name for path in MODULE_FILES]}"
    )

MODULE_PATH = MODULE_FILES[0]
spec = importlib.util.spec_from_file_location("gribs", MODULE_PATH)

if spec is None or spec.loader is None:
    raise ImportError(f"Could not load GRIBS module from {MODULE_PATH}")

gribs = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = gribs
spec.loader.exec_module(gribs)


@pytest.mark.integration
def test_default_self_tests_pass():
    result = subprocess.run(
        [sys.executable, str(MODULE_PATH), "--selftest"],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, (
        f"Self-test failed with exit code {result.returncode}\n"
        f"stdout:\n{result.stdout}\n"
        f"stderr:\n{result.stderr}"
    )
    assert "-> all checks passed" in result.stdout


@pytest.mark.integration
def test_cd_nozzle_transition_classification(tmp_path):
    """C-D nozzle mass flow passes T4 while thrust jump remains diagnostic."""
    source = ROOT / "gribs_config.json"
    supplied = json.loads(source.read_text(encoding="utf-8"))
    supplied["nozzle"]["expansion_ratio"] = 4.0

    config_path = tmp_path / "gribs_config_eps4.json"
    config_path.write_text(
        json.dumps(supplied, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    result = subprocess.run(
        [
            sys.executable,
            str(MODULE_PATH),
            "--config",
            str(config_path),
            "--selftest",
        ],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, (
        "Ae/At = 4 self-test should pass after A-2 classification.\n"
        f"exit code: {result.returncode}\n"
        f"stdout:\n{result.stdout}\n"
        f"stderr:\n{result.stderr}"
    )

    diagnostic_lines = [
        line
        for line in result.stdout.splitlines()
        if "C-D nozzle transition mass-flow continuity" in line
    ]
    assert len(diagnostic_lines) == 1, (
        "Expected exactly one C-D nozzle transition diagnostic, "
        f"found {len(diagnostic_lines)}:\n"
        + "\n".join(diagnostic_lines)
    )

    diagnostic = diagnostic_lines[0]
    assert "[PASS]" in diagnostic
    assert "Ae/At = 4" in diagnostic
    assert "known jump-model limitation" in diagnostic

    match = re.search(
        r"mdot jump = ([0-9.eE+-]+), F jump = ([0-9.eE+-]+)",
        diagnostic,
    )
    assert match is not None, (
        f"Could not parse C-D nozzle transition diagnostics:\n{diagnostic}"
    )

    mdot_jump = float(match.group(1))
    thrust_jump = float(match.group(2))

    assert mdot_jump == pytest.approx(4.02e-08, rel=1.0e-2)
    assert thrust_jump == pytest.approx(1.04, rel=1.0e-2)
    assert "-> all checks passed" in result.stdout
    assert "[FAIL]" not in result.stdout


@pytest.mark.fast
def test_geometry_identity_finite_difference():
    c = gribs.Config()
    x = 0.5 * gribs.web_thickness(c)
    h = gribs.web_thickness(c) * 1e-6
    derivative = (
        gribs.geometry(x + h, c)[3]
        - gribs.geometry(x - h, c)[3]
    ) / (2 * h)
    area = gribs.geometry(x, c)[2]
    assert abs(derivative / area - 1.0) < 1e-7


@pytest.mark.fast
def test_example_configuration_keys_are_valid():
    path = ROOT / "examples" / "example_config.json"
    supplied = json.loads(path.read_text(encoding="utf-8"))
    # The example config must be loadable and valid
    configuration = gribs.configuration_from_document(supplied)
    assert configuration.config is not None


@pytest.mark.fast
def test_transition_policy_defaults_to_jump():
    """Omitting transition_policy preserves the current jump-model behavior."""
    source = ROOT / "gribs_config.json"
    supplied = json.loads(source.read_text(encoding="utf-8"))
    supplied["nozzle"].pop("transition_policy", None)

    configuration = gribs.configuration_from_document(supplied)

    assert configuration.config.transition_policy == "jump"


@pytest.mark.fast
@pytest.mark.parametrize("policy", ["jump", "shock"])
def test_transition_policy_accepts_known_values(policy):
    """The configuration schema accepts the documented transition policies."""
    source = ROOT / "gribs_config.json"
    supplied = json.loads(source.read_text(encoding="utf-8"))
    supplied["nozzle"]["transition_policy"] = policy

    configuration = gribs.configuration_from_document(supplied)

    assert configuration.config.transition_policy == policy


@pytest.mark.fast
def test_transition_policy_rejects_unknown_value():
    """Unknown policies must fail rather than silently selecting a model."""
    source = ROOT / "gribs_config.json"
    supplied = json.loads(source.read_text(encoding="utf-8"))
    supplied["nozzle"]["transition_policy"] = "interpolate"

    with pytest.raises(
        gribs.ConfigurationError,
        match=r"transition_policy.*jump.*shock",
    ):
        gribs.configuration_from_document(supplied)


@pytest.mark.fast
def test_shock_transition_policy_is_reserved_for_future_model():
    """The reserved shock policy must not silently use the jump model."""
    c = gribs.Config(transition_policy="shock")

    with pytest.raises(
        NotImplementedError,
        match=r"transition_policy.*shock.*not implemented",
    ):
        gribs.require_implemented_transition_policy(c)


@pytest.mark.fast
def test_jump_transition_policy_is_implemented():
    """The current direct-switch model remains available."""
    c = gribs.Config(transition_policy="jump")

    assert gribs.require_implemented_transition_policy(c) is None


@pytest.mark.fast
def test_run_model_rejects_unimplemented_shock_transition_policy():
    """The normal calculation path must reject the reserved shock model."""
    source = ROOT / "gribs_config.json"
    supplied = json.loads(source.read_text(encoding="utf-8"))
    supplied["nozzle"]["transition_policy"] = "shock"
    configuration = gribs.configuration_from_document(supplied)

    with pytest.raises(
        NotImplementedError,
        match=r"transition_policy.*shock.*not implemented",
    ):
        gribs.run_model(configuration.config)

@pytest.mark.fast
def test_masked_trapezoid_does_not_bridge_disjoint_regions():
    """Masked integration must not invent area across gaps."""
    t = gribs.np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    y = gribs.np.ones_like(t)
    mask = gribs.np.array([True, True, False, True, True])

    integral = gribs._masked_trapezoid(y, t, mask)

    assert integral == pytest.approx(2.0)


@pytest.mark.fast
def test_nozzle_transition_diagnostics_reports_contiguous_band_cost():
    """Transition diagnostics report time and impulse without bridging gaps."""
    c = gribs.Config(
        eps_nozzle=4.0,
        transition_policy="jump",
        use_separation=True,
    )
    t = gribs.np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    h = {
        "t": t,
        "p0": gribs.np.array([2.0, 3.0, 10.0, 3.0, 2.0]),
        "At": gribs.np.ones_like(t),
        "F": gribs.np.ones_like(t),
    }

    diagnostic = gribs.nozzle_transition_diagnostics(
        h,
        c,
        Ae=4.0,
        lower_pressure_Pa=1.0,
        upper_pressure_Pa=5.0,
    )

    assert diagnostic["applicable"] is True
    assert diagnostic["detected"] is True
    assert diagnostic["duration_s"] == pytest.approx(2.0)
    assert diagnostic["impulse_Ns"] == pytest.approx(2.0)
    assert diagnostic["interval_count"] == 2


@pytest.mark.fast
def test_nozzle_transition_summary_lines_include_warning_code():
    """Transition text output preserves units and the warning identifier."""
    transition = {
        "policy": "jump",
        "applicable": True,
        "detected": True,
        "warning_code": "W_NOZZLE_TRANSITION",
        "interval_count": 1,
        "duration_s": 0.02844,
        "duration_fraction": 0.012,
        "impulse_Ns": 0.6370,
        "impulse_fraction": 0.00041,
    }
    warnings = [
        {
            "code": "W_NOZZLE_TRANSITION",
            "severity": "warning",
            "message": "Synthetic transition warning for text testing.",
            "context": {},
        }
    ]

    lines = gribs.nozzle_transition_summary_lines(transition, warnings)
    text = "\n".join(lines)

    assert "[ nozzle transition diagnostic ]" in text
    assert "28.4400 ms" in text
    assert "0.041000 %" in text
    assert "WARNING [W_NOZZLE_TRANSITION]" in text

# =============================================================================
# Phase A-4: mass and energy ledger
# =============================================================================

@pytest.mark.fast
@pytest.mark.parametrize(
    ("simulation_end_s", "expected_mass_kg"),
    [
        (0.0, 0.0),
        (0.025, 0.0005),
        (0.1, 0.002),
        (0.5, 0.002),
    ],
)
def test_igniter_injected_mass_is_limited_by_simulation_duration(
    simulation_end_s,
    expected_mass_kg,
):
    """Actual igniter mass uses the simulated interval, not only commanded charge."""
    c = gribs.Config(
        ign_mdot=0.02,
        ign_time=0.1,
    )

    actual = gribs.igniter_injected_mass(c, simulation_end_s)

    assert actual == pytest.approx(expected_mass_kg)


@pytest.mark.fast
def test_igniter_injected_mass_is_zero_when_igniter_is_disabled():
    """Zero flow or zero duration produces no injected igniter mass."""
    assert gribs.igniter_injected_mass(
        gribs.Config(ign_mdot=0.0, ign_time=0.1),
        1.0,
    ) == pytest.approx(0.0)

    assert gribs.igniter_injected_mass(
        gribs.Config(ign_mdot=0.02, ign_time=0.0),
        1.0,
    ) == pytest.approx(0.0)


@pytest.mark.fast
def test_mass_energy_ledger_separates_igniter_from_propellant():
    """The product ledger must not report igniter charge as propellant error."""
    c = gribs.Config(
        rho_p=1000.0,
        R_i0=0.01,
        R_p=0.02,
        L_p0=0.1,
        n_end=0.0,
        ign_mdot=0.02,
        ign_time=0.1,
    )

    initial_propellant = gribs.propellant_mass(c)
    initial_chamber = 0.01
    igniter_injected = 0.002
    generated_propellant = initial_propellant
    generated_total = generated_propellant + igniter_injected
    final_chamber = 0.02
    discharged_total = (
        initial_chamber
        + generated_total
        - final_chamber
    )

    t = gribs.np.array([0.0, 0.05, 0.1])
    xw = gribs.web_thickness(c)
    h = {
        "t": t,
        "phase": gribs.np.array(
            ["burn", "burn", "blowdown"],
            dtype=object,
        ),
        "x": gribs.np.array([0.0, 0.5 * xw, xw]),
        "m_gen": gribs.np.array(
            [0.0, 0.5 * generated_total, generated_total]
        ),
        "m_igniter": gribs.np.array(
            [0.0, 0.5 * igniter_injected, igniter_injected]
        ),
        "m_out": gribs.np.array(
            [0.0, 0.5 * discharged_total, discharged_total]
        ),
        "m_total_eos": gribs.np.array(
            [initial_chamber, 0.03, final_chamber]
        ),
        "mdot_out_total": gribs.np.array(
            [discharged_total / 0.1] * 3
        ),
        "mdot_out_gas": gribs.np.array(
            [0.75 * discharged_total / 0.1] * 3
        ),
        "mdot_out_condensed": gribs.np.array(
            [0.25 * discharged_total / 0.1] * 3
        ),
    }
    res = {
        "m_total0": initial_chamber,
    }

    ledger = gribs.mass_energy_ledger(res, h, c)
    mass = ledger["mass"]

    assert mass["basis"] == "total_product_mass"
    assert mass["igniter_mass_commanded_kg"] == pytest.approx(0.002)
    assert mass["igniter_mass_injected_kg"] == pytest.approx(0.002)
    assert mass["generated_mass_total_kg"] == pytest.approx(
        generated_total
    )
    assert mass["generated_mass_propellant_kg"] == pytest.approx(
        initial_propellant
    )
    assert mass["initial_propellant_mass_kg"] == pytest.approx(
        initial_propellant
    )
    assert mass["final_unburned_propellant_kg"] == pytest.approx(0.0)
    assert mass["burned_propellant_geometry_kg"] == pytest.approx(
        initial_propellant
    )

    assert mass["total_system"]["signed_residual_kg"] == pytest.approx(
        0.0,
        abs=1.0e-14,
    )
    assert mass["total_system"]["relative_residual"] == pytest.approx(
        0.0,
        abs=1.0e-14,
    )
    assert mass["total_system"]["tolerance"] == pytest.approx(1.0e-6)
    assert mass["total_system"]["passed"] is True

    assert mass["propellant"]["signed_residual_kg"] == pytest.approx(
        0.0,
        abs=1.0e-14,
    )
    assert mass["propellant"]["relative_residual"] == pytest.approx(
        0.0,
        abs=1.0e-14,
    )
    assert mass["propellant"]["tolerance"] == pytest.approx(
        max(1.0e-9, 100.0 * c.rtol)
    )
    assert mass["propellant"]["passed"] is True

    phase = mass["phase_discharge_diagnostic"]
    assert phase["method"] == "sampled_history_trapezoid"
    assert phase["grid_dependent"] is True
    assert phase["gas_mass_kg"] == pytest.approx(
        0.75 * discharged_total
    )
    assert phase["condensed_mass_kg"] == pytest.approx(
        0.25 * discharged_total
    )
    assert phase["total_quadrature_mass_kg"] == pytest.approx(
        discharged_total
    )
    assert phase["ode_total_mass_kg"] == pytest.approx(
        discharged_total
    )
    assert phase["relative_quadrature_mismatch"] == pytest.approx(
        0.0,
        abs=1.0e-14,
    )

    assert mass["passed"] is True
    assert ledger["energy"]["available"] is False
    assert ledger["energy"]["closure_residual_available"] is False
    assert "independent energy state" in ledger["energy"]["reason"]
    assert ledger["passed"] is True


@pytest.mark.fast
def test_mass_energy_ledger_records_failed_acceptance_without_hiding_residual():
    """Ledger failures retain signed and normalized residual information."""
    c = gribs.Config(
        rho_p=1000.0,
        R_i0=0.01,
        R_p=0.02,
        L_p0=0.1,
        n_end=0.0,
        ign_mdot=0.0,
        ign_time=0.0,
    )

    initial_propellant = gribs.propellant_mass(c)
    total_residual = 2.0e-6
    propellant_residual = 2.0e-8 * initial_propellant

    h = {
        "t": gribs.np.array([0.0, 1.0]),
        "phase": gribs.np.array(["burn", "burn"], dtype=object),
        "x": gribs.np.array([0.0, 0.0]),
        "m_gen": gribs.np.array(
            [0.0, initial_propellant + propellant_residual]
        ),
        "m_igniter": gribs.np.array([0.0, 0.0]),
        "m_out": gribs.np.array(
            [0.0, initial_propellant - total_residual]
        ),
        "m_total_eos": gribs.np.array([0.0, 0.0]),
        "mdot_out_total": gribs.np.array(
            [initial_propellant - total_residual] * 2
        ),
        "mdot_out_gas": gribs.np.array(
            [initial_propellant - total_residual] * 2
        ),
        "mdot_out_condensed": gribs.np.array([0.0, 0.0]),
    }
    res = {
        "m_total0": 0.0,
    }

    ledger = gribs.mass_energy_ledger(res, h, c)
    mass = ledger["mass"]

    assert mass["total_system"]["signed_residual_kg"] > 0.0
    assert mass["total_system"]["relative_residual"] > 1.0e-6
    assert mass["total_system"]["passed"] is False

    assert mass["propellant"]["signed_residual_kg"] > 0.0
    assert mass["propellant"]["relative_residual"] > 1.0e-8
    assert mass["propellant"]["passed"] is False

    assert mass["passed"] is False
    assert ledger["passed"] is False


@pytest.mark.fast
def test_sampled_generation_rate_includes_igniter_during_blowdown():
    """Sampled mdot_gen_total must use the same igniter rule as the ODE."""
    class SyntheticSolution:
        def __init__(self, t0, t1, state):
            self.t = gribs.np.array([t0, t1], dtype=float)
            self._state = gribs.np.asarray(state, dtype=float)

        def sol(self, times):
            return gribs.np.repeat(
                self._state[:, None],
                len(times),
                axis=1,
            )

    c = gribs.Config(
        rho_p=1000.0,
        R_i0=0.01,
        R_p=0.02,
        L_p0=0.1,
        n_end=0.0,
        R_t0=0.0015,
        eps_nozzle=1.0,
        p_a=101325.0,
        p0_init=200000.0,
        ign_mdot=0.02,
        ign_time=1.0,
    )
    xw = gribs.web_thickness(c)
    state = gribs.np.array(
        [
            c.p0_init,
            xw,
            c.R_t0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
        dtype=float,
    )
    solution = SyntheticSolution(0.2, 0.3, state)
    result = {
        "Ae": c.eps_nozzle * gribs.math.pi * c.R_t0**2,
        "m_total0": 0.0,
        "sol_burn": solution,
        "sol_blow": solution,
    }

    original_phase_grid = gribs._phase_grid
    original_gas_props = gribs.gas_props
    original_phase_fractions = gribs.phase_fractions
    original_psi_and_deriv = gribs.psi_and_deriv
    original_nozzle_state = gribs.nozzle_state
    original_is_extrapolated = gribs.is_extrapolated
    original_equilibrium_pressure = gribs.equilibrium_pressure

    try:
        gribs._phase_grid = lambda sol, n_lin, log_head: sol.t
        gribs.gas_props = lambda p, config: (300.0, 1000.0, 1.2)
        gribs.phase_fractions = lambda p, config: (1.0, 0.0)
        gribs.psi_and_deriv = lambda p, config: (300000.0, 0.0)
        gribs.nozzle_state = lambda p, At, Ae, config: {
            "mdot_total": 0.0,
            "mdot_gas": 0.0,
            "mdot_condensed": 0.0,
            "F": 0.0,
            "Me": 0.0,
            "pe": config.p_a,
            "ve": 0.0,
            "regime": "subsonic",
        }
        gribs.is_extrapolated = lambda p, config: False
        gribs.equilibrium_pressure = lambda Ab, Ri, At, Ae, config: None

        history = gribs.sample(result, c)
    finally:
        gribs._phase_grid = original_phase_grid
        gribs.gas_props = original_gas_props
        gribs.phase_fractions = original_phase_fractions
        gribs.psi_and_deriv = original_psi_and_deriv
        gribs.nozzle_state = original_nozzle_state
        gribs.is_extrapolated = original_is_extrapolated
        gribs.equilibrium_pressure = original_equilibrium_pressure

    blowdown = history["phase"] == "blowdown"
    assert gribs.np.any(blowdown)
    assert history["mdot_gen_total"][blowdown] == pytest.approx(
        c.ign_mdot
    )


@pytest.mark.fast
def test_summarize_preserves_legacy_metrics_and_adds_ledger(monkeypatch):
    """A-4 adds results.ledger without changing the A-3 metric definitions."""
    sentinel_ledger = {
        "schema_version": "a4-ledger-1",
        "mass": {"passed": True},
        "energy": {
            "available": False,
            "closure_residual_available": False,
        },
        "passed": True,
    }

    monkeypatch.setattr(
        gribs,
        "mass_energy_ledger",
        lambda res, h, c: sentinel_ledger,
    )

    source = ROOT / "gribs_config.json"
    supplied = json.loads(source.read_text(encoding="utf-8"))
    configuration = gribs.configuration_from_document(supplied)
    c = configuration.config
    c.eta_T0 = 0.98

    class Backend:
        name = gribs.BACKEND_CEA_PYTHON
        pmin = 1.0e5
        pmax = 8.0e6
        supports_two_phase = True
        # An active two-phase backend must state its Yg range explicitly
        # (summarize() rejects a missing range with
        # E_PHASE_FRACTION_VALIDATION_INVALID instead of assuming Yg = 1).
        property_validation = {
            "gas_mass_fraction_min": 1.0,
            "gas_mass_fraction_max": 1.0,
        }
        cache_schema = 4
        mass_basis_semantics_version = 1

    monkeypatch.setattr(gribs, "thermochemistry", lambda: Backend())
    monkeypatch.setattr(
        gribs,
        "choke_limit_pressure",
        lambda config, At0, Ae: None,
    )
    monkeypatch.setattr(
        gribs,
        "nozzle_transition_diagnostics",
        lambda history, config, Ae: {
            "policy": "jump",
            "applicable": False,
            "detected": False,
            "warning_code": None,
            "interval_count": 0,
            "duration_s": 0.0,
            "duration_fraction": 0.0,
            "impulse_Ns": 0.0,
            "impulse_fraction": 0.0,
        },
    )
    monkeypatch.setattr(
        gribs,
        "nozzle_state",
        lambda p, At, Ae, config: {
            "mdot_total": 1.0,
        },
    )

    m_prop = gribs.propellant_mass(c)
    t = gribs.np.array([0.0, 1.0])
    At0 = gribs.math.pi * c.R_t0**2
    history = {
        "t": t,
        "phase": gribs.np.array(["burn", "burn"], dtype=object),
        "p0": gribs.np.array([2.0e5, 3.0e5]),
        "F": gribs.np.array([0.0, 1.0]),
        "impulse": gribs.np.array([0.0, 1.0]),
        "m_out": gribs.np.array([0.0, m_prop]),
        "m_gen": gribs.np.array([0.0, m_prop]),
        "At": gribs.np.array([At0, At0]),
        "Kn": gribs.np.array([1.0, 1.0]),
        "Rt": gribs.np.array([c.R_t0, c.R_t0]),
        "m_total_eos": gribs.np.array([1.0, 1.0]),
        "m_total_balance": gribs.np.array([1.0, 1.0]),
        "Y_gas": gribs.np.array([1.0, 1.0]),
        "Y_condensed": gribs.np.array([0.0, 0.0]),
        "extrap": gribs.np.array([False, False]),
        "regime": gribs.np.array(
            ["subsonic", "subsonic"],
            dtype=object,
        ),
    }
    result = {
        "m_total0": 1.0,
        "At0": At0,
        "Ae": At0,
        "stop_burn": "time-limit",
        "stop_blow": "not-run",
        "init_regime": "subsonic",
        "t_events": {},
    }

    summary = gribs.summarize(result, history, c)

    assert summary["ledger"] is sentinel_ledger
    assert any(
        warning["code"] == "W_THERMAL_EFFICIENCY_PARTIAL"
        for warning in summary["warnings"]
    )
    assert summary["propellant_mass_balance_error"] == pytest.approx(0.0)
    assert summary["total_mass_consistency_error"] == pytest.approx(0.0)
    assert summary["total_mass_consistency_error_invariant"] == pytest.approx(
        0.0
    )
    assert summary["total_mass_consistency_reference_mass_kg"] == pytest.approx(
        1.0
    )
    assert (
        summary["total_mass_consistency_metric_version"]
        == "run_wide_balance_mass_reference_v1"
    )
    assert summary["post_run_consistency"]["passed"] is True


@pytest.mark.fast
def test_summary_text_contains_mass_and_energy_ledger_sections():
    """Human-readable output exposes the A-4 ledger and its limitations."""
    c = gribs.Config()
    summary = {
        "web_mm": 1.0,
        "propellant_mass_g": 1.0,
        "thermochemistry_backend": "synthetic",
        "property_extrapolation": False,
        "property_extrapolation_points": 0,
        "property_extrapolation_fraction": 0.0,
        "property_range_min_Pa": 1.0e5,
        "property_range_max_Pa": 8.0e6,
        "pressure_min_Pa": 2.0e5,
        "pressure_max_Pa": 3.0e5,
        "two_phase_model": {
            "mode": "homogeneous_equilibrium",
            "active": True,
        },
        "maximum_condensed_mass_fraction_during_run": 0.0,
        "mean_condensed_mass_fraction_during_burn": 0.0,
        "condensed_mass_fraction_at_max_pressure": 0.0,
        "condensed_mass_fraction_table_min": 0.0,
        "condensed_mass_fraction_table_max": 0.0,
        "min_pressure_factor_1_minus_pPsi_over_Psi": 1.0,
        "nozzle_entrainment_assumption": "synthetic",
        "condensed_volume_assumption": "synthetic",
        "initial_regime": "subsonic",
        "choke_limit_pressure_Pa": None,
        "initial_mdot_gen_kg_s": 0.0,
        "initial_mdot_out_kg_s": 1.0,
        "event_times": {},
        "burn_stop_reason": "time-limit",
        "blowdown_stop_reason": "not-run",
        "burn_time_s": 1.0,
        "total_time_s": 1.0,
        "p_max_Pa": 3.0e5,
        "t_p_max_s": 1.0,
        "p_min_burn_Pa": 2.0e5,
        "t_p_min_s": 0.0,
        "p_mean_burn_Pa": 2.5e5,
        "F_max_N": 1.0,
        "t_F_max_s": 1.0,
        "F_mean_burn_N": 1.0,
        "impulse_burn_Ns": 1.0,
        "impulse_blowdown_Ns": 0.0,
        "impulse_total_Ns": 1.0,
        "Isp_s": 1.0,
        "cstar_eff_m_s": 1.0,
        "CF_eff": 1.0,
        "Kn_initial": 1.0,
        "Kn_max": 1.0,
        "Kn_final": 1.0,
        "throat_radius_final_mm": 1.5,
        "regimes_visited": ["subsonic"],
        "nozzle_transition": {
            "policy": "jump",
            "applicable": False,
            "detected": False,
            "warning_code": None,
            "interval_count": 0,
            "duration_s": 0.0,
            "duration_fraction": 0.0,
            "impulse_Ns": 0.0,
            "impulse_fraction": 0.0,
        },
        "warnings": [
            {
                "code": "W_THERMAL_EFFICIENCY_PARTIAL",
                "severity": "warning",
                "message": "Direct CEA comparisons are not like-for-like.",
                "context": {"temperature_efficiency_eta_T0": 0.98},
            }
        ],
        "ledger": {
            "schema_version": "a4-ledger-1",
            "mass": {
                "basis": "total_product_mass",
                "initial_chamber_inventory_kg": 0.01,
                "generated_mass_total_kg": 1.002,
                "generated_mass_propellant_kg": 1.0,
                "igniter_mass_injected_kg": 0.002,
                "igniter_mass_commanded_kg": 0.002,
                "discharged_mass_total_kg": 1.0,
                "final_chamber_inventory_kg": 0.012,
                "initial_propellant_mass_kg": 1.0,
                "final_unburned_propellant_kg": 0.0,
                "burned_propellant_geometry_kg": 1.0,
                "total_system": {
                    "signed_residual_kg": 0.0,
                    "relative_residual": 0.0,
                    "tolerance": 1.0e-6,
                    "passed": True,
                },
                "propellant": {
                    "signed_residual_kg": 0.0,
                    "relative_residual": 0.0,
                    "tolerance": 1.0e-9,
                    "passed": True,
                },
                "phase_discharge_diagnostic": {
                    "method": "sampled_history_trapezoid",
                    "grid_dependent": True,
                    "gas_mass_kg": 0.75,
                    "condensed_mass_kg": 0.25,
                    "relative_quadrature_mismatch": 0.0,
                },
                "passed": True,
            },
            "energy": {
                "available": False,
                "closure_residual_available": False,
                "reason": "No independent energy state is integrated.",
                "model": {
                    "equilibrium_constraint": (
                        "HP (assigned enthalpy and pressure)"
                    ),
                    "temperature_treatment": "T0 = eta_T0 * T_CEA",
                    "state_relation": "p0 * Vg = m_total * Psi",
                },
            },
            "passed": True,
        },
        "propellant_mass_balance_error": 0.002,
        "total_mass_consistency_error": 1.0e-10,
        "total_mass_consistency_error_invariant": 1.0e-12,
        "total_mass_consistency_reference_mass_kg": 1.0,
        "total_mass_consistency_metric_version": (
            "run_wide_balance_mass_reference_v1"
        ),
        "post_run_consistency": {
            "passed": True,
        },
    }

    summary["ledger"]["passed"] = False
    summary["ledger"]["mass"]["passed"] = False
    summary["ledger"]["mass"]["propellant"].update(
        {
            "passed": False,
            "signed_residual_kg": 2.0e-8,
            "relative_residual": 4.0e-8,
            "tolerance": 1.0e-9,
        }
    )
    gribs.append_run_verdict_warnings(summary)
    text = gribs.summary_text(summary, c)

    assert "[ verification warnings ]" in text
    assert "W_MASS_LEDGER_FAILED" in text
    assert "[ mass ledger ]" in text
    assert "propellant relative residual" in text
    assert "4.000e-08 (limit 1.0e-09, FAIL)" in text
    assert "mass-ledger result" in text
    assert "[ energy ledger ]" in text
    assert "W_THERMAL_EFFICIENCY_PARTIAL" in text
    assert "Direct CEA comparisons are not like-for-like." in text
    assert "closure residual available" in text
    assert "not a complete energy-conservation law" in text
    assert "legacy propellant balance metric" in text
    assert "includes igniter in the legacy numerator" in text
    assert "EOS-history invariant error" in text
    assert "legacy shrinking-mass error" in text
    assert "(deprecated)" in text


@pytest.mark.fast
def test_rhs_exposes_integrated_igniter_mass_as_appended_state(monkeypatch):
    """The appended ODE quadrature records the existing igniter mass source."""
    c = gribs.Config(
        rho_p=1000.0,
        R_i0=0.01,
        R_p=0.02,
        L_p0=0.1,
        n_end=0.0,
        R_t0=0.0015,
        eps_nozzle=1.0,
        p_ref=1.0e6,
        p_a=101325.0,
        ign_mdot=0.02,
        ign_time=0.1,
    )
    Ae = c.eps_nozzle * gribs.math.pi * c.R_t0**2

    monkeypatch.setattr(
        gribs,
        "nozzle_state",
        lambda p0, At, exit_area, config: {
            "mdot_total": 0.0,
            "F": 0.0,
        },
    )
    monkeypatch.setattr(
        gribs,
        "psi_and_deriv",
        lambda p0, config: (1.0, 0.0),
    )
    monkeypatch.setattr(
        gribs,
        "burn_rate",
        lambda p0, config, Ab, Ri: 0.0,
    )

    rhs = gribs.make_rhs(c, Ae, burning=True)
    y = gribs.np.array(
        [
            c.p0_init,
            0.0,
            c.R_t0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
        dtype=float,
    )

    before_cutoff = rhs(0.05, y)
    after_cutoff = rhs(0.2, y)

    assert len(before_cutoff) == 7
    assert len(after_cutoff) == 7
    assert before_cutoff[gribs.IMG] == pytest.approx(c.ign_mdot)
    assert before_cutoff[gribs.IMIGN] == pytest.approx(c.ign_mdot)
    assert after_cutoff[gribs.IMG] == pytest.approx(0.0)
    assert after_cutoff[gribs.IMIGN] == pytest.approx(0.0)


@pytest.mark.fast
def test_mass_energy_ledger_uses_integrated_igniter_state():
    """Ledger separation uses the ODE igniter quadrature, not analytic charge."""
    c = gribs.Config(
        rho_p=1000.0,
        R_i0=0.01,
        R_p=0.02,
        L_p0=0.1,
        n_end=0.0,
        ign_mdot=0.02,
        ign_time=0.1,
    )

    initial_propellant = gribs.propellant_mass(c)
    # Deliberately distinct from the commanded 0.002 kg charge so the
    # test detects whether the ledger uses the ODE quadrature.
    integrated_igniter = 0.0019
    generated_total = initial_propellant + integrated_igniter

    h = {
        "t": gribs.np.array([0.0, 0.1]),
        "phase": gribs.np.array(["burn", "burn"], dtype=object),
        "x": gribs.np.array([0.0, gribs.web_thickness(c)]),
        "m_gen": gribs.np.array([0.0, generated_total]),
        "m_igniter": gribs.np.array([0.0, integrated_igniter]),
        "m_out": gribs.np.array([0.0, generated_total]),
        "m_total_eos": gribs.np.array([0.0, 0.0]),
        "mdot_out_total": gribs.np.array(
            [generated_total / 0.1, generated_total / 0.1]
        ),
        "mdot_out_gas": gribs.np.array(
            [generated_total / 0.1, generated_total / 0.1]
        ),
        "mdot_out_condensed": gribs.np.array([0.0, 0.0]),
    }
    res = {"m_total0": 0.0}

    ledger = gribs.mass_energy_ledger(res, h, c)
    mass = ledger["mass"]

    assert mass["igniter_mass_commanded_kg"] == pytest.approx(0.002)
    assert mass["igniter_mass_injected_kg"] == pytest.approx(
        integrated_igniter
    )
    assert mass["generated_mass_propellant_kg"] == pytest.approx(
        initial_propellant
    )
    assert mass["propellant"]["relative_residual"] == pytest.approx(
        0.0,
        abs=1.0e-14,
    )
    assert mass["propellant"]["passed"] is True


@pytest.mark.fast
def test_default_solver_relative_tolerance_is_1e_10():
    """The product default supports the A-4 propellant-ledger criterion."""
    assert gribs.Config().rtol == pytest.approx(
        1.0e-10,
        rel=0.0,
        abs=0.0,
    )


@pytest.mark.fast
def test_example_configuration_uses_default_solver_relative_tolerance():
    """The distributed JSON example and runtime default stay aligned."""
    document = json.loads(
        (ROOT / "gribs_config.json").read_text(encoding="utf-8")
    )
    loaded = gribs.configuration_from_document(document)

    assert document["solver"]["relative_tolerance"] == pytest.approx(
        1.0e-10,
        rel=0.0,
        abs=0.0,
    )
    assert loaded.config.rtol == pytest.approx(
        gribs.Config().rtol,
        rel=0.0,
        abs=0.0,
    )


@pytest.mark.fast
def test_post_run_consistency_uses_time_invariant_reference_mass():
    """The A-5.1-2d gate must not shrink with final chamber inventory."""
    result = {
        "m_total0": 1.0,
    }
    history = {
        "m_gen": gribs.np.array([0.0, 0.0]),
        "m_out": gribs.np.array([0.0, 0.9989995]),
        "m_total_eos": gribs.np.array([1.0, 0.001]),
        "m_total_balance": gribs.np.array([1.0, 0.0010005]),
    }

    report = gribs.post_run_consistency_report(result, history)
    eos = report["eos_history"]

    assert eos["reference_mass_kg"] == pytest.approx(1.0)
    assert eos["invariant_relative_error"] == pytest.approx(5.0e-7)
    assert eos["legacy_instantaneous_mass_relative_error"] == pytest.approx(
        5.0e-4
    )
    assert eos["legacy_metric_deprecated"] is True
    assert eos["metric_version"] == "run_wide_balance_mass_reference_v1"
    assert eos["passed"] is True
    assert report["final_total_mass"]["passed"] is True
    assert report["passed"] is True


@pytest.mark.fast
def test_post_run_consistency_rejects_mismatched_history_shapes():
    result = {
        "m_total0": 1.0,
    }
    history = {
        "m_gen": gribs.np.array([0.0, 0.0]),
        "m_out": gribs.np.array([0.0, 0.0]),
        "m_total_eos": gribs.np.array([1.0, 1.0]),
        "m_total_balance": gribs.np.array([1.0]),
    }

    with pytest.raises(
        ValueError,
        match="m_total_balance must match m_total_eos",
    ):
        gribs.post_run_consistency_report(result, history)


@pytest.mark.fast
def test_run_verdict_warning_exposes_failed_ledger_context():
    summary = {
        "warnings": [],
        "post_run_consistency": {
            "passed": True,
        },
        "ledger": {
            "passed": False,
            "mass": {
                "passed": False,
                "igniter_mass_commanded_kg": 0.002,
                "igniter_mass_injected_kg": 0.002,
                "total_system": {
                    "passed": True,
                    "signed_residual_kg": 1.0e-10,
                    "relative_residual": 1.0e-9,
                    "tolerance": 1.0e-6,
                },
                "propellant": {
                    "passed": False,
                    "signed_residual_kg": 2.0e-8,
                    "relative_residual": 4.0e-8,
                    "tolerance": 1.0e-9,
                },
            },
        },
    }

    gribs.append_run_verdict_warnings(summary)
    gribs.append_run_verdict_warnings(summary)

    assert gribs.run_verdict_passed(summary) is False
    assert len(summary["warnings"]) == 1

    warning = summary["warnings"][0]
    assert warning["code"] == "W_MASS_LEDGER_FAILED"
    assert warning["severity"] == "error"
    assert warning["context"]["total_system_passed"] is True
    assert warning["context"]["total_system_relative_residual"] == pytest.approx(
        1.0e-9
    )
    assert warning["context"]["total_system_tolerance"] == pytest.approx(
        1.0e-6
    )
    assert warning["context"]["propellant_passed"] is False
    assert warning["context"]["propellant_relative_residual"] == pytest.approx(
        4.0e-8
    )
    assert warning["context"]["propellant_tolerance"] == pytest.approx(
        1.0e-9
    )
    assert warning["context"]["igniter_mass_commanded_kg"] == pytest.approx(
        0.002
    )
    assert warning["context"]["igniter_mass_injected_kg"] == pytest.approx(
        0.002
    )


@pytest.mark.fast
@pytest.mark.parametrize(
    "ledger_passed,post_run_passed,warnings,allow,expected",
    [
        pytest.param(True, True, [], False, 0, id="all-pass"),
        pytest.param(False, True, [], False, 5, id="ledger-only-fail"),
        pytest.param(False, True, [], True, 0, id="ledger-only-allowed"),
        pytest.param(True, False, [], True, 1, id="post-run-not-overridden"),
        pytest.param(False, False, [], True, 1, id="combined-failure-not-overridden"),
        pytest.param(
            True,
            True,
            [{"code": "E_OUT_OF_VALIDATED_RANGE"}],
            True,
            1,
            id="scope-error-not-overridden",
        ),
    ],
)
def test_run_verification_exit_code_obeys_q18_contract(
    ledger_passed,
    post_run_passed,
    warnings,
    allow,
    expected,
):
    summary = {
        "warnings": warnings,
        "post_run_consistency": {"passed": post_run_passed},
        "ledger": {"passed": ledger_passed},
    }

    assert gribs.run_verification_exit_code(
        summary,
        allow_verification_failure=allow,
    ) == expected


@pytest.mark.fast
@pytest.mark.parametrize(
    "allow,expected_exit",
    [(False, 5), (True, 0)],
    ids=["strict-ledger-exit", "explicit-ledger-override"],
)
def test_main_writes_diagnostics_before_ledger_exit(
    tmp_path,
    monkeypatch,
    capsys,
    allow,
    expected_exit,
):
    document = json.loads(
        (ROOT / "gribs_config.json").read_text(encoding="utf-8")
    )
    output_dir = tmp_path / "ledger-failure-output"
    document["output"]["directory"] = str(output_dir)
    config_path = tmp_path / "ledger-failure-config.json"
    config_path.write_text(
        json.dumps(document, indent=2) + "\n",
        encoding="utf-8",
    )

    backend = SimpleNamespace(
        name="synthetic-ledger-probe",
        metadata=lambda: {"backend": "synthetic-ledger-probe"},
        ideal_mixture_density=1000.0,
        packing_fraction=0.6,
        bulk_mixture_density=600.0,
    )
    summary = {
        "warnings": [],
        "post_run_consistency": {"passed": True},
        "ledger": {
            "passed": False,
            "mass": {
                "passed": False,
                "igniter_mass_commanded_kg": 0.002,
                "igniter_mass_injected_kg": 0.002,
                "total_system": {
                    "passed": True,
                    "signed_residual_kg": 1.0e-10,
                    "relative_residual": 1.0e-9,
                    "tolerance": 1.0e-6,
                },
                "propellant": {
                    "passed": False,
                    "signed_residual_kg": 2.0e-8,
                    "relative_residual": 4.0e-8,
                    "tolerance": 1.0e-9,
                },
            },
        },
        "two_phase_model": {"mode": "synthetic", "active": False},
    }

    monkeypatch.setattr(
        gribs, "initialize_thermochemistry", lambda *_args: backend
    )
    monkeypatch.setattr(gribs, "_print_backend_banner", lambda *_args: None)
    monkeypatch.setattr(gribs, "model_applicability_error", lambda *_args: None)
    monkeypatch.setattr(gribs, "run_model", lambda _config: object())
    monkeypatch.setattr(gribs, "sample", lambda *_args: {})
    monkeypatch.setattr(
        gribs,
        "summarize",
        lambda *_args, **_kwargs: summary,
    )
    monkeypatch.setattr(gribs, "environment_info", lambda: {})
    monkeypatch.setattr(
        gribs,
        "_cea_theoretical_diagnostics",
        lambda *_args: None,
    )
    monkeypatch.setattr(
        gribs,
        "summary_text",
        lambda item, *_args: (
            f"{item['warnings'][0]['code']} "
            f"{item['warnings'][0]['context']['propellant_relative_residual']} "
            f"{item['warnings'][0]['context']['propellant_tolerance']}"
        ),
    )
    monkeypatch.setattr(
        gribs,
        "write_csv",
        lambda _history, path: path.write_text("history\n", encoding="utf-8"),
    )
    monkeypatch.setattr(gribs, "make_figure", lambda *_args: None)
    monkeypatch.setattr(
        gribs,
        "self_tests_post_run",
        lambda *_args, **_kwargs: None,
    )

    argv = ["--config", str(config_path)]
    if allow:
        argv.append("--allow-verification-failure")
    exit_code = gribs.main(argv)

    assert exit_code == expected_exit
    assert (output_dir / "time_history.csv").exists()
    summary_text_path = output_dir / "summary.txt"
    summary_json_path = output_dir / "summary.json"
    assert summary_text_path.exists()
    assert summary_json_path.exists()
    assert "W_MASS_LEDGER_FAILED" in summary_text_path.read_text(
        encoding="utf-8"
    )
    record = json.loads(summary_json_path.read_text(encoding="utf-8"))
    warning = record["results"]["warnings"][0]
    assert warning["code"] == "W_MASS_LEDGER_FAILED"
    assert warning["context"]["propellant_relative_residual"] == pytest.approx(
        4.0e-8
    )
    captured = capsys.readouterr()
    if allow:
        assert "allowed by --allow-verification-failure" in captured.err
    else:
        assert "exit 5" not in captured.err
        assert "W_MASS_LEDGER_FAILED" in captured.err


@pytest.mark.fast
def test_run_verdict_warning_exposes_failed_post_run_context():
    summary = {
        "warnings": [],
        "post_run_consistency": {
            "passed": False,
            "tolerance": 1.0e-6,
            "final_total_mass": {
                "passed": True,
                "relative_error": 1.0e-9,
            },
            "eos_history": {
                "passed": False,
                "invariant_relative_error": 2.0e-6,
                "reference_mass_kg": 1.0,
                "legacy_instantaneous_mass_relative_error": 3.0e-4,
            },
        },
        "ledger": {
            "passed": True,
        },
    }

    gribs.append_run_verdict_warnings(summary)
    gribs.append_run_verdict_warnings(summary)

    assert gribs.run_verdict_passed(summary) is False
    assert len(summary["warnings"]) == 1
    assert summary["warnings"][0]["code"] == "W_POST_RUN_CONSISTENCY"
    assert (
        summary["warnings"][0]["context"][
            "eos_history_invariant_relative_error"
        ]
        == pytest.approx(2.0e-6)
    )


@pytest.mark.fast
def test_model_applicability_error_uses_the_out_of_scope_threshold():
    config = gribs.Config()
    backend = type("BackendProbe", (), {})()
    backend.supports_two_phase = True
    backend.property_validation = {
        "gas_mass_fraction_min": 0.65,
        "gas_mass_fraction_max": 1.0,
    }
    backend.pmin = 1.0e5
    backend.pmax = 8.0e6

    error = gribs.model_applicability_error(config, backend)

    assert error["code"] == "E_OUT_OF_VALIDATED_RANGE"
    assert error["severity"] == "error"
    assert error["context"]["thermochemistry_table_Yg_range"] == [0.65, 1.0]
    assert error["context"]["thermochemistry_pressure_range_Pa"] == [
        1.0e5,
        8.0e6,
    ]

    backend.property_validation["gas_mass_fraction_min"] = 0.6500001
    assert gribs.model_applicability_error(config, backend) is None

    config.two_phase_mode = gribs.TWO_PHASE_SINGLE_LEGACY
    backend.property_validation["gas_mass_fraction_min"] = 0.5
    assert gribs.model_applicability_error(config, backend) is None

    config.two_phase_mode = gribs.TWO_PHASE_HOMOGENEOUS
    backend.supports_two_phase = False
    assert gribs.model_applicability_error(config, backend) is None


@pytest.mark.fast
def test_main_fails_fast_for_out_of_scope_two_phase_table(
    tmp_path,
    monkeypatch,
    capsys,
):
    document = json.loads(
        (ROOT / "gribs_config.json").read_text(encoding="utf-8")
    )
    output_dir = tmp_path / "must_not_be_created"
    document["output"]["directory"] = str(output_dir)
    config_path = tmp_path / "out_of_scope.json"
    config_path.write_text(
        json.dumps(document, indent=2) + "\n",
        encoding="utf-8",
    )

    backend = type("BackendProbe", (), {})()
    backend.supports_two_phase = True
    backend.property_validation = {
        "gas_mass_fraction_min": 0.65,
        "gas_mass_fraction_max": 1.0,
    }
    backend.pmin = 1.0e5
    backend.pmax = 8.0e6
    monkeypatch.setattr(
        gribs,
        "initialize_thermochemistry",
        lambda _config, _cache_dir: backend,
    )
    monkeypatch.setattr(gribs, "_print_backend_banner", lambda *_args: None)
    monkeypatch.setattr(
        gribs,
        "run_model",
        lambda *_args: pytest.fail("out-of-scope run reached the solver"),
    )

    exit_code = gribs.main(["--config", str(config_path)])
    captured = capsys.readouterr()

    assert exit_code == 1
    assert "E_OUT_OF_VALIDATED_RANGE" in captured.err
    json_record = json.loads(
        next(line for line in captured.err.splitlines() if line.startswith("{"))
    )
    assert json_record["context"]["simulation_started"] is False
    assert not output_dir.exists()


@pytest.mark.fast
def test_out_of_scope_applicability_error_fails_run_verdict():
    summary = {
        "warnings": [],
        "two_phase_model": {
            "active": True,
            "mode": "homogeneous_equilibrium",
            "nozzle_entrainment": "complete",
        },
        "validation": {
            "model_status": "out_of_scope",
            "validated_Yg_range": [0.65, 1.0],
            "validated_pressure_range_Pa": [1.0e5, 8.0e6],
        },
        "gas_mass_fraction_table_min": 0.65,
        "gas_mass_fraction_table_max": 1.0,
        "condensed_mass_fraction_table_max": 0.35,
        "maximum_condensed_mass_fraction_during_run": 0.3,
        "post_run_consistency": {"passed": True},
        "ledger": {"passed": True},
    }

    gribs.append_model_applicability_warnings(summary)

    assert [warning["code"] for warning in summary["warnings"]] == [
        "E_OUT_OF_VALIDATED_RANGE",
        "W_TWO_PHASE_NOZZLE_APPROXIMATION",
    ]
    assert gribs.run_verdict_passed(summary) is False


@pytest.mark.fast
def test_model_applicability_warnings_cover_condensed_scope_and_nozzle_limits():
    summary = {
        "warnings": [],
        "two_phase_model": {
            "active": True,
            "mode": "homogeneous_equilibrium",
            "nozzle_entrainment": "complete",
        },
        "validation": {
            "model_status": "approximate_with_limitations",
        },
        "gas_mass_fraction_table_min": 0.72,
        "condensed_mass_fraction_table_max": 0.28,
        "maximum_condensed_mass_fraction_during_run": 0.25,
    }

    gribs.append_model_applicability_warnings(summary)
    warnings = {warning["code"]: warning for warning in summary["warnings"]}

    assert set(warnings) == {
        "W_HIGH_CONDENSATE_MODEL_LIMITATION",
        "W_TWO_PHASE_NOZZLE_APPROXIMATION",
    }
    high_condensate = warnings["W_HIGH_CONDENSATE_MODEL_LIMITATION"]
    assert high_condensate["severity"] == "warning"
    assert high_condensate["context"]["model_status"] == (
        "approximate_with_limitations"
    )
    assert high_condensate["context"]["gas_mass_fraction_table_min"] == (
        pytest.approx(0.72)
    )
    assert high_condensate["context"][
        "validated_in_scope_requires_Yg_min_greater_than"
    ] == pytest.approx(0.99)

    nozzle = warnings["W_TWO_PHASE_NOZZLE_APPROXIMATION"]
    assert nozzle["severity"] == "warning"
    assert nozzle["context"]["particle_slip_modelled"] is False
    assert nozzle["context"][
        "maximum_condensed_mass_fraction_during_run"
    ] == pytest.approx(0.25)

    # Emission is idempotent if a caller applies the reporting helper twice.
    gribs.append_model_applicability_warnings(summary)
    assert len(summary["warnings"]) == 2


@pytest.mark.fast
@pytest.mark.parametrize(
    "backend_name,eta_t0,expected_warning,direction",
    [
        (gribs.BACKEND_CEA_PYTHON, 0.98, True, "below"),
        (gribs.BACKEND_CEA_PYTHON, 1.02, True, "above"),
        (gribs.BACKEND_CEA_LEGACY_EXECUTABLE, 0.97, True, "below"),
        (gribs.BACKEND_CEA_PYTHON, 1.0, False, None),
        (gribs.BACKEND_REFERENCE_CONSTANT, 0.98, False, None),
    ],
    ids=["cea-scaled-down", "cea-scaled-up", "legacy-cea", "eta-one", "reference"],
)
def test_thermal_efficiency_partial_warning_scope_and_verdict(
    backend_name,
    eta_t0,
    expected_warning,
    direction,
):
    summary = {
        "warnings": [],
        "post_run_consistency": {"passed": True},
        "ledger": {"passed": True},
    }
    config = gribs.Config(eta_T0=eta_t0)
    backend = type("BackendProbe", (), {"name": backend_name})()

    gribs.append_thermal_efficiency_warning(summary, config, backend)

    if not expected_warning:
        assert summary["warnings"] == []
        return

    assert len(summary["warnings"]) == 1
    warning = summary["warnings"][0]
    assert warning["code"] == "W_THERMAL_EFFICIENCY_PARTIAL"
    assert warning["severity"] == "warning"
    assert warning["context"]["temperature_efficiency_eta_T0"] == eta_t0
    assert warning["context"]["cea_state_re_solved_after_scaling"] is False
    assert "density" in warning["context"][
        "direct_cea_comparisons_not_like_for_like"
    ]
    assert f"{direction} its equilibrium value" in warning["message"]
    assert gribs.run_verdict_passed(summary) is True

    gribs.append_thermal_efficiency_warning(summary, config, backend)
    assert len(summary["warnings"]) == 1


@pytest.mark.fast
@pytest.mark.parametrize(
    "two_phase_active,model_status,condensed_max,expected_codes",
    [
        (False, "approximate_with_limitations", 0.25, set()),
        (True, "validated_in_scope", 4.0e-16, set()),
        (
            True,
            "approximate_with_limitations",
            4.0e-16,
            {"W_HIGH_CONDENSATE_MODEL_LIMITATION"},
        ),
        (
            True,
            "out_of_scope",
            4.0e-16,
            {"E_OUT_OF_VALIDATED_RANGE"},
        ),
    ],
)
def test_model_applicability_warning_scope_and_roundoff_threshold(
    two_phase_active,
    model_status,
    condensed_max,
    expected_codes,
):
    yg_min = 0.64 if model_status == "out_of_scope" else 0.998
    summary = {
        "warnings": [],
        "two_phase_model": {
            "active": two_phase_active,
            "mode": "homogeneous_equilibrium",
            "nozzle_entrainment": "complete",
        },
        "validation": {
            "model_status": model_status,
            "validated_Yg_range": [yg_min, 1.0],
            "validated_pressure_range_Pa": [1.0e5, 8.0e6],
        },
        "gas_mass_fraction_table_min": yg_min,
        "condensed_mass_fraction_table_max": 1.0 - yg_min,
        "maximum_condensed_mass_fraction_during_run": condensed_max,
    }

    gribs.append_model_applicability_warnings(summary)

    assert {warning["code"] for warning in summary["warnings"]} == (
        expected_codes
    )


@pytest.mark.fast
def test_model_applicability_warnings_do_not_fail_run_verdict():
    summary = {
        "warnings": [],
        "two_phase_model": {
            "active": True,
            "mode": "homogeneous_equilibrium",
            "nozzle_entrainment": "complete",
        },
        "validation": {
            "model_status": "approximate_with_limitations",
        },
        "gas_mass_fraction_table_min": 0.72,
        "condensed_mass_fraction_table_max": 0.28,
        "maximum_condensed_mass_fraction_during_run": 0.25,
        "post_run_consistency": {"passed": True},
        "ledger": {"passed": True},
    }

    gribs.append_model_applicability_warnings(summary)

    assert len(summary["warnings"]) == 2
    assert gribs.run_verdict_passed(summary) is True


@pytest.mark.fast
def test_default_self_test_name_contract_is_exact_and_ordered():
    configuration = gribs.load_json_configuration(ROOT / "gribs_config.json")
    expected = gribs.expected_self_test_names(configuration.config)

    assert len(expected) == 15
    assert len(set(expected)) == 15
    assert (
        expected[3]
        == "converging-nozzle continuity (mdot, F)"
    )


@pytest.mark.fast
def test_cd_self_test_name_contract_uses_transition_test():
    supplied = json.loads(
        (ROOT / "gribs_config.json").read_text(encoding="utf-8")
    )
    supplied["nozzle"]["expansion_ratio"] = 4.0
    configuration = gribs.configuration_from_document(supplied)

    expected = gribs.expected_self_test_names(configuration.config)

    assert len(expected) == 15
    assert len(set(expected)) == 15
    assert (
        expected[3]
        == "C-D nozzle transition mass-flow continuity"
    )
    assert "converging-nozzle continuity (mdot, F)" not in expected


@pytest.mark.fast
def test_self_test_name_contract_detects_missing_unexpected_and_duplicate():
    c = gribs.Config()
    expected = list(gribs.expected_self_test_names(c))

    valid = gribs.self_test_name_contract(c, expected)
    assert valid["passed"] is True
    assert valid["expected_count"] == 15
    assert valid["actual_count"] == 15

    invalid_names = expected[:-1] + [expected[0], "unexpected test"]
    invalid = gribs.self_test_name_contract(c, invalid_names)

    assert invalid["passed"] is False
    assert expected[-1] in invalid["missing_names"]
    assert "unexpected test" in invalid["unexpected_names"]
    assert expected[0] in invalid["duplicate_names"]


@pytest.mark.fast
def test_reference_constant_phase_fraction_contract():
    supplied = json.loads(
        (
            ROOT
            / "cases"
            / "reference"
            / "level1_end_faces_0.json"
        ).read_text(encoding="utf-8")
    )
    configuration = gribs.configuration_from_document(supplied)
    c = configuration.config

    backend = gribs.ReferenceConstantBackend(
        c,
        c.reactants_spec_list,
        c.packing_fraction,
    )

    assert backend.phase_fractions(c.p0_init) == (1.0, 0.0)

    metadata = backend.metadata()
    assert metadata["gas_mass_fraction"] == pytest.approx(1.0)
    assert metadata["condensed_mass_fraction"] == pytest.approx(0.0)

    table = backend.thermo_dump_table()
    assert table["gas_mass_fraction"] == pytest.approx([1.0, 1.0])
    assert table["condensed_mass_fraction"] == pytest.approx([0.0, 0.0])


def _cache_reader_probe(tmp_path, table):
    """Build the minimal backend shell needed by _read_cache()."""
    backend = gribs.PressureTableBackend.__new__(
        gribs.PressureTableBackend
    )
    backend.cache_path = tmp_path / "cache-schema-probe.json"
    backend.key = "current-schema-probe-key"
    backend.name = "cache-schema-probe"
    backend.cache_path.write_text(
        json.dumps(table, indent=2) + "\n",
        encoding="utf-8",
    )
    return backend


@pytest.mark.fast
def test_cache_schema_identity_rejects_missing_schema(tmp_path):
    backend = _cache_reader_probe(tmp_path, {})

    with pytest.raises(
        gribs.ThermochemistryCacheError,
        match="no schema identity",
    ):
        backend._read_cache()


@pytest.mark.fast
def test_cache_schema_identity_rejects_noninteger_schema(tmp_path):
    backend = _cache_reader_probe(
        tmp_path,
        {"schema": "4"},
    )

    with pytest.raises(
        gribs.ThermochemistryCacheError,
        match="schema must be an integer",
    ):
        backend._read_cache()


@pytest.mark.fast
def test_cache_schema_identity_rejects_schema_3_mass_basis(tmp_path):
    backend = _cache_reader_probe(
        tmp_path,
        {"schema": 3},
    )

    with pytest.raises(
        gribs.ThermochemistryCacheError,
        match="schema 3 uses obsolete CEA mass-basis semantics",
    ) as exc_info:
        backend._read_cache()

    assert exc_info.value.code == "E_OLD_THERMO_CACHE"
    assert exc_info.value.context["schema_actual"] == 3
    assert exc_info.value.context["schema_expected"] == 4
    assert exc_info.value.context["cache_path"] == str(backend.cache_path)


class _CacheBuildReached(RuntimeError):
    pass


class _LegacyCacheProbeBackend(gribs.PressureTableBackend):
    name = "legacy-cache-probe"
    supports_two_phase = False
    table_stem = "thermo_legacy_cache_probe"
    equilibrium_formulation = "cache identity test"

    def _build_table(self):
        raise _CacheBuildReached("new-cache rebuild path reached")


def _legacy_cache_probe_inputs(tmp_path, **config_overrides):
    config = gribs.Config(
        thermo_backend=_LegacyCacheProbeBackend.name,
        two_phase_mode=gribs.TWO_PHASE_SINGLE_LEGACY,
        **config_overrides,
    )
    reactants = [gribs.ReactantSpec(
        name="TEST",
        wt_percent=100.0,
        temperature_K=298.15,
        density_kg_m3=1000.0,
    )]
    backend = _LegacyCacheProbeBackend.__new__(_LegacyCacheProbeBackend)
    backend.c = config
    backend.reactants = reactants
    backend.cache_root = tmp_path
    legacy_key = backend._legacy_schema3_cache_key()
    legacy_path = backend._legacy_schema3_cache_path()
    return config, reactants, legacy_path, legacy_key


def _write_legacy_schema3_cache(path, backend_name, cache_key):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps({
            "schema": 3,
            "backend": backend_name,
            "cache_key": cache_key,
        }),
        encoding="utf-8",
    )


@pytest.mark.fast
def test_matching_legacy_schema3_cache_raises_structured_error(tmp_path):
    config, reactants, legacy_path, legacy_key = _legacy_cache_probe_inputs(tmp_path)
    schema3_payload = {
        "cache_schema": 3,
        "backend": _LegacyCacheProbeBackend.name,
        "reactants": [[r.name, float(r.wt_percent), float(r.temperature_K)]
                      for r in reactants],
        "equilibrium": _LegacyCacheProbeBackend.equilibrium_formulation,
        "pressure_min_Pa": float(config.p_fit_min),
        "pressure_max_Pa": float(config.p_fit_max),
        "pressure_points": int(config.cea_pressure_points),
        "temperature_efficiency": float(config.eta_T0),
        "phase_fractions_tabulated": False,
        "molecular_weight_definition": (
            "gas_phase_M (gas-phase M; total MW is diagnostic only since "
            "v0.5.0-alpha)"
        ),
        "options": {},
    }
    expected_legacy_key = gribs.hashlib.sha256(
        json.dumps(schema3_payload, sort_keys=True, default=str).encode("utf-8")
    ).hexdigest()[:24]
    assert legacy_key == expected_legacy_key

    _write_legacy_schema3_cache(
        legacy_path,
        _LegacyCacheProbeBackend.name,
        legacy_key,
    )

    with pytest.raises(gribs.ThermochemistryCacheError) as exc_info:
        _LegacyCacheProbeBackend(config, reactants, 1.0, tmp_path)

    error = exc_info.value
    assert error.code == "E_OLD_THERMO_CACHE"
    assert error.context["backend"] == _LegacyCacheProbeBackend.name
    assert error.context["cache_path"] == str(legacy_path)
    assert error.context["schema_actual"] == 3
    assert error.context["schema_expected"] == 4
    assert error.context["detection"] == "matching_schema3_cache_key"
    assert error.context["cache_key_schema3"] == legacy_key
    assert error.context["cache_key_current"] != legacy_key
    assert "E_OLD_THERMO_CACHE" in str(error)
    assert "Expected schema: 4; actual schema: 3" in str(error)
    assert str(legacy_path) in str(error)


@pytest.mark.fast
def test_main_reports_old_cache_error_as_structured_stderr(
    tmp_path,
    monkeypatch,
    capsys,
):
    loaded = gribs.Configuration(
        config=gribs.Config(),
        output_dir=tmp_path,
        output_names={},
        document={},
        migration_notes=[],
    )
    monkeypatch.setattr(gribs, "load_json_configuration", lambda _path: loaded)
    monkeypatch.setattr(gribs, "validate", lambda _config: None)
    error = gribs.ThermochemistryCacheError(
        "E_OLD_THERMO_CACHE: obsolete test cache",
        code="E_OLD_THERMO_CACHE",
        context={
            "cache_path": str(tmp_path / "old-cache.json"),
            "schema_actual": 3,
            "schema_expected": 4,
        },
    )

    def raise_old_cache_error(_config, _cache_dir):
        raise error

    monkeypatch.setattr(
        gribs,
        "initialize_thermochemistry",
        raise_old_cache_error,
    )

    assert gribs.main(["--config", "ignored.json"]) == 3
    captured = capsys.readouterr()
    assert "Thermochemistry backend error [E_OLD_THERMO_CACHE]" in captured.err
    record = json.loads(captured.err.splitlines()[-1])
    assert record["code"] == "E_OLD_THERMO_CACHE"
    assert record["severity"] == "error"
    assert record["context"]["schema_actual"] == 3
    assert record["context"]["schema_expected"] == 4


@pytest.mark.fast
def test_unrelated_schema3_cache_does_not_block_rebuild(tmp_path):
    config, reactants, legacy_path, _ = _legacy_cache_probe_inputs(tmp_path)
    unrelated_key = "f" * 24
    unrelated_path = tmp_path / (
        f"{_LegacyCacheProbeBackend.table_stem}_{unrelated_key}.json")
    assert unrelated_path != legacy_path
    _write_legacy_schema3_cache(
        unrelated_path,
        _LegacyCacheProbeBackend.name,
        unrelated_key,
    )

    with pytest.raises(_CacheBuildReached, match="new-cache rebuild path reached"):
        _LegacyCacheProbeBackend(config, reactants, 1.0, tmp_path)


@pytest.mark.fast
@pytest.mark.parametrize(
    "config_overrides",
    [
        {"cea_rebuild_cache": True},
        {"cea_cache_enabled": False},
    ],
    ids=["explicit-rebuild", "cache-disabled"],
)
def test_legacy_schema3_cache_guard_respects_rebuild_escape(
    tmp_path,
    config_overrides,
):
    config, reactants, legacy_path, legacy_key = _legacy_cache_probe_inputs(
        tmp_path,
        **config_overrides,
    )
    _write_legacy_schema3_cache(
        legacy_path,
        _LegacyCacheProbeBackend.name,
        legacy_key,
    )

    with pytest.raises(_CacheBuildReached, match="new-cache rebuild path reached"):
        _LegacyCacheProbeBackend(config, reactants, 1.0, tmp_path)


@pytest.mark.fast
def test_cache_schema_identity_accepts_schema_4_identity(tmp_path):
    backend = _cache_reader_probe(
        tmp_path,
        {"schema": 4},
    )

    with pytest.raises(
        gribs.ThermochemistryCacheError,
        match="cache is missing",
    ) as exc_info:
        backend._read_cache()

    message = str(exc_info.value)
    assert "schema identity" not in message
    assert "incompatible with required schema" not in message
    assert "schema 3 uses obsolete" not in message


def _mass_basis_cache_probe(tmp_path, **table_overrides):
    backend = gribs.PressureTableBackend.__new__(gribs.PressureTableBackend)
    backend.name = "mass-basis-invariant-probe"
    backend.key = "mass-basis-invariant-key"
    backend.supports_two_phase = True
    backend.cache_path = tmp_path / "thermo_mass_basis_probe.json"
    backend.c = gribs.Config(
        cea_pressure_points=4,
        eta_T0=0.97,
        two_phase_mode=gribs.TWO_PHASE_HOMOGENEOUS,
    )
    pressure = [1.0e5, 2.0e5, 4.0e5, 8.0e5]
    temperature = [3000.0, 3100.0, 3200.0, 3300.0]
    molecular_weight = [25.0] * 4
    gas_constant = [gribs.R_UNIVERSAL / 25.0] * 4
    table = {
        "schema": 4,
        "backend": backend.name,
        "cache_key": backend.key,
        "mass_basis_semantics_version": 1,
        "pressure_Pa": pressure,
        "temperature_K": temperature,
        "gamma_s": [1.2] * 4,
        "equilibrium_molecular_weight_kg_kmol": molecular_weight,
        "n_gas_per_kg_total": [1.0 / 25.0] * 4,
        "effective_gas_constant_J_kgK": gas_constant,
        "gas_mass_fraction": [0.7, 0.8, 0.9, 1.0],
        "condensed_mass_fraction": [0.3, 0.2, 0.1, 0.0],
        "psi_J_kg": [
            r * t * backend.c.eta_T0
            for r, t in zip(gas_constant, temperature)
        ],
    }
    table.update(table_overrides)
    backend.cache_path.write_text(json.dumps(table), encoding="utf-8")
    return backend, table


@pytest.mark.fast
def test_mass_basis_cache_invariants_accept_exact_roundtrip(tmp_path):
    backend, _ = _mass_basis_cache_probe(tmp_path)
    table = backend._read_cache()
    assert table["mass_basis_semantics_version"] == 1


@pytest.mark.fast
def test_mass_basis_cache_invariants_allow_pre_tag_schema4_cache(tmp_path):
    backend, table = _mass_basis_cache_probe(tmp_path)
    table.pop("mass_basis_semantics_version")
    table.pop("n_gas_per_kg_total")
    backend.cache_path.write_text(json.dumps(table), encoding="utf-8")

    assert backend._read_cache()["cache_key"] == backend.key


@pytest.mark.fast
@pytest.mark.parametrize(
    "pressure, invariant",
    [
        ([1.0e5, 0.0, 4.0e5, 8.0e5], "pressure_grid_positive"),
        ([1.0e5, 2.0e5, 2.0e5, 8.0e5], "pressure_grid_strictly_increasing"),
        ([1.0e5, 2.0e5, float("nan"), 8.0e5], "pressure_grid_finite"),
        ([[1.0e5], [2.0e5], [4.0e5], [8.0e5]], "pressure_grid_shape"),
    ],
    ids=["nonpositive", "not-increasing", "nonfinite", "wrong-shape"],
)
def test_mass_basis_cache_invariants_reject_bad_pressure_grid(
    tmp_path,
    pressure,
    invariant,
):
    backend, table = _mass_basis_cache_probe(tmp_path)
    table["pressure_Pa"] = pressure
    backend.cache_path.write_text(json.dumps(table), encoding="utf-8")

    with pytest.raises(gribs.ThermochemistryInvariantError) as exc_info:
        backend._read_cache()

    assert exc_info.value.code == "E_THERMO_MASS_BASIS_INVARIANT_FAILED"
    assert exc_info.value.context["invariant"] == invariant


@pytest.mark.fast
@pytest.mark.parametrize(
    "column, index, multiplier, invariant",
    [
        ("effective_gas_constant_J_kgK", 1, 1.0 + 2.0e-12,
         "R_eff_times_M_CEA_equals_Ru"),
        ("n_gas_per_kg_total", 2, 1.0 + 2.0e-12,
         "M_CEA_times_n_gas_per_kg_total_equals_one"),
        ("psi_J_kg", 3, 1.0 + 2.0e-12,
         "psi_equals_R_eff_times_T_used"),
        ("gas_mass_fraction", 3, 1.0 + 2.0e-12,
         "phase_fraction_physical_ranges"),
        ("condensed_mass_fraction", 2, 1.0,
         "gas_plus_condensed_mass_fraction_equals_one"),
    ],
    ids=["R-times-M", "M-times-n", "psi", "phase-range", "phase-sum"],
)
def test_mass_basis_cache_invariants_reject_corrupted_columns(
    tmp_path,
    column,
    index,
    multiplier,
    invariant,
):
    backend, table = _mass_basis_cache_probe(tmp_path)
    values = list(table[column])
    if column == "condensed_mass_fraction":
        values[index] += 5.0e-13
    else:
        values[index] *= multiplier
    table[column] = values
    backend.cache_path.write_text(json.dumps(table), encoding="utf-8")

    with pytest.raises(gribs.ThermochemistryInvariantError) as exc_info:
        backend._read_cache()

    assert exc_info.value.code == "E_THERMO_MASS_BASIS_INVARIANT_FAILED"
    assert exc_info.value.context["invariant"] == invariant
    assert exc_info.value.context["cache_path"] == str(backend.cache_path)
    assert str(exc_info.value).startswith("E_THERMO_MASS_BASIS_INVARIANT_FAILED:")


@pytest.mark.fast
def test_gamma_s_below_unity_reports_invalid_pressure_intervals(tmp_path):
    backend, table = _mass_basis_cache_probe(tmp_path)
    table["gamma_s"] = [0.9985, 1.2, 0.999, 1.1]
    backend.cache_path.write_text(json.dumps(table), encoding="utf-8")

    with pytest.raises(gribs.GammaBelowUnityError) as exc_info:
        backend._read_cache()

    error = exc_info.value
    assert error.code == "E_GAMMA_S_BELOW_UNITY"
    assert error.context["table_source"] == "cache"
    assert error.context["invalid_point_count"] == 2
    assert error.context["invalid_pressure_intervals_Pa"] == [
        [1.0e5, 1.0e5],
        [4.0e5, 4.0e5],
    ]
    assert error.context["invalid_gamma_s_range"] == pytest.approx(
        [0.9985, 0.999]
    )
    assert "Do not clip these values" in str(error)


@pytest.mark.fast
def test_gamma_s_below_unity_generated_scan_reports_configured_range(tmp_path):
    backend, _ = _mass_basis_cache_probe(tmp_path)
    pressure = [1.0e5, 2.0e5, 4.0e5, 8.0e5]

    with pytest.raises(gribs.GammaBelowUnityError) as exc_info:
        backend._validate_gamma_s_above_unity(
            pressure,
            [1.2, 1.1, 0.9985, 0.9984],
            table_source="generated",
        )

    context = exc_info.value.context
    assert context["table_source"] == "generated"
    assert context["pressure_scan_range_Pa"] == [1.0e5, 8.0e5]
    assert context["configured_pressure_range_Pa"] == [
        backend.c.p_fit_min,
        backend.c.p_fit_max,
    ]
    assert context["invalid_pressure_intervals_Pa"] == [[4.0e5, 8.0e5]]


@pytest.mark.fast
def test_props_reports_gamma_below_unity_with_pressure_context(tmp_path):
    backend = gribs.PressureTableBackend.__new__(gribs.PressureTableBackend)
    backend.name = "interpolated-gamma-probe"
    backend.key = "interpolated-gamma-key"
    backend.cache_path = tmp_path / "interpolated-gamma.json"
    backend.c = gribs.Config(property_policy="clamp")
    backend.pmin = 1.0e5
    backend.pmax = 8.0e5
    backend.Ri = lambda _z: 300.0
    backend.Ti = lambda _z: 3000.0
    backend.gi = lambda _z: 0.9985

    with pytest.raises(gribs.GammaBelowUnityError) as exc_info:
        backend.props(2.0e5)

    assert exc_info.value.code == "E_GAMMA_S_BELOW_UNITY"
    assert exc_info.value.context["table_source"] == "interpolated"
    assert exc_info.value.context["invalid_pressure_intervals_Pa"] == [
        [2.0e5, 2.0e5]
    ]


@pytest.mark.fast
def test_cea_builder_raises_gamma_error_at_first_invalid_grid_point(monkeypatch):
    backend = gribs.CEAPythonBackend.__new__(gribs.CEAPythonBackend)
    backend.cea = SimpleNamespace(
        R=gribs.R_UNIVERSAL,
        HP="HP",
        __version__="test-cea",
        __file__="test-cea-module",
        lib_version=lambda: "test-cea-library",
    )
    backend.c = gribs.Config(
        p_fit_min=1.0e5,
        p_fit_max=8.0e5,
        cea_pressure_points=4,
        eta_T0=1.0,
    )
    backend.cache_path = Path("/tmp/gribs-gamma-builder-probe.json")
    backend.key = "gamma-builder-probe"
    solution = SimpleNamespace()

    class FakeSolver:
        def __init__(self):
            self.calls = 0

        def solve(self, result, _hp, _h_r, _pressure_bar, _weights):
            gamma_values = [1.2, 0.9985, 1.1, 1.15]
            result.converged = True
            result.last_error = 0
            result.T = 3000.0
            result.M = 25.0
            result.MW = 25.0
            result.gamma_s = gamma_values[self.calls]
            result.density = 1.0
            result.cp_eq = 1.0
            result.sonic_velocity = 1000.0
            result.enthalpy = 0.0
            result.mass_fractions = {"X": 1.0}
            result.n = 1.0 / 25.0
            self.calls += 1

    solver = FakeSolver()
    backend._build_mixtures = lambda: (
        None,
        None,
        solver,
        solution,
        gribs.np.array([1.0]),
        0.0,
    )
    monkeypatch.setattr(
        backend,
        "_phase_split",
        lambda _product, _solver, _solution: (0.7, 0.3, {}),
    )

    with pytest.raises(gribs.GammaBelowUnityError) as exc_info:
        backend._build_table()

    assert exc_info.value.context["table_source"] == "generated"
    intervals = exc_info.value.context["invalid_pressure_intervals_Pa"]
    assert len(intervals) == 1
    assert intervals[0] == pytest.approx([2.0e5, 2.0e5])
    assert exc_info.value.context["pressure_scan_range_Pa"] == pytest.approx(
        [1.0e5, 2.0e5]
    )
    assert solver.calls == 2


@pytest.mark.fast
def test_mass_basis_cache_semantics_version_mismatch_is_old_cache(tmp_path):
    backend, _ = _mass_basis_cache_probe(
        tmp_path,
        mass_basis_semantics_version=2,
    )

    with pytest.raises(gribs.ThermochemistryCacheError) as exc_info:
        backend._read_cache()

    assert exc_info.value.code == "E_OLD_THERMO_CACHE"
    assert exc_info.value.context["value_actual"] == 2
    assert exc_info.value.context["value_expected"] == 1
    assert str(backend.cache_path) in str(exc_info.value)


@pytest.mark.fast
def test_mass_basis_semantics_tag_requires_n_column(tmp_path):
    backend, table = _mass_basis_cache_probe(tmp_path)
    table.pop("n_gas_per_kg_total")
    backend.cache_path.write_text(json.dumps(table), encoding="utf-8")

    with pytest.raises(gribs.ThermochemistryInvariantError) as exc_info:
        backend._read_cache()

    assert exc_info.value.context["invariant"] == "column_present:n_gas_per_kg_total"


@pytest.mark.fast
@pytest.mark.parametrize(
    "error_class,code,context",
    [
        (
            gribs.ThermochemistryInvariantError,
            "E_THERMO_MASS_BASIS_INVARIANT_FAILED",
            {"invariant": "psi_equals_R_eff_times_T_used"},
        ),
        (
            gribs.GammaBelowUnityError,
            "E_GAMMA_S_BELOW_UNITY",
            {"invalid_pressure_intervals_Pa": [[4.0e5, 4.0e5]]},
        ),
    ],
    ids=["mass-basis-invariant", "gamma-below-unity"],
)
def test_main_reports_thermochemistry_error_as_structured_stderr(
    tmp_path,
    monkeypatch,
    capsys,
    error_class,
    code,
    context,
):
    loaded = gribs.Configuration(
        config=gribs.Config(),
        output_dir=tmp_path,
        output_names={},
        document={},
        migration_notes=[],
    )
    monkeypatch.setattr(gribs, "load_json_configuration", lambda _path: loaded)
    monkeypatch.setattr(gribs, "validate", lambda _config: None)
    error = error_class(f"{code}: test mismatch", context=context)

    def raise_thermochemistry_error(_config, _cache_dir):
        raise error

    monkeypatch.setattr(
        gribs,
        "initialize_thermochemistry",
        raise_thermochemistry_error,
    )

    assert gribs.main(["--config", "ignored.json"]) == 3
    captured = capsys.readouterr()
    assert f"Thermochemistry backend error [{code}]" in captured.err
    record = json.loads(captured.err.splitlines()[-1])
    assert record["code"] == code
    assert record["context"] == context


@pytest.mark.fast
def test_corrected_mass_basis_metadata_contract():
    backend = gribs.PressureTableBackend.__new__(
        gribs.PressureTableBackend
    )
    backend.name = "metadata-contract-probe"
    backend.reactants = []
    backend.supports_two_phase = True
    backend.equilibrium_formulation = "HP equilibrium"
    backend.c = type(
        "MetadataProbeConfig",
        (),
        {
            "p_fit_min": 1.0e5,
            "p_fit_max": 8.0e6,
            "cea_pressure_points": 4,
            "eta_T0": 1.0,
            "two_phase_mode": "homogeneous_equilibrium",
        },
    )()

    payload = backend.cache_key_payload()
    definition = payload["molecular_weight_definition"]

    assert payload["cache_schema"] == 4
    assert payload["mass_basis_semantics_version"] == 1
    assert "two_phase_mode" in payload
    assert payload["two_phase_mode"] == "homogeneous_equilibrium"
    assert "M_CEA = 1 / EqSolution.n" in definition
    assert "R_eff = Ru / M_CEA" in definition
    assert "no additional Yg factor" in definition
    assert "gas-phase M" not in definition

    source = MODULE_PATH.read_text(encoding="utf-8")
    assert "R_eff(p) = Ru / M_CEA(p) = Yg(p) * R_gas(p)" in source
    assert "no additional Yg factor" in source
    assert "p/(rho*T) = R_eff" in source

    forbidden_claims = (
        "EqSolution.M (gas-phase, kg/kmol)",
        "EqSolution.M - the GAS-PHASE molecular weight",
        "gas_phase_M (gas-phase M;",
    )
    assert all(claim not in source for claim in forbidden_claims)


def test_migration_model_change_warning_code_in_source():
    """The W_MIGRATION_MODEL_CHANGE warning code must be emitted when a
    v0.4.3-alpha configuration is migrated to homogeneous_equilibrium mode.
    This test verifies the code string is present and wired into main()."""
    source = MODULE_PATH.read_text(encoding="utf-8")
    assert "W_MIGRATION_MODEL_CHANGE" in source
    # The warning must be appended to the summary warnings list, not just
    # printed to stderr.
    assert 'summ.setdefault("warnings"' in source or 'summ["warnings"]' in source


# ---------------------------------------------------------------------------
# Phase-fraction validation for the two-phase applicability gate.
# A missing / non-numeric / non-finite / out-of-range Yg range on an ACTIVE
# two-phase backend is a structured E_PHASE_FRACTION_VALIDATION_INVALID error,
# distinct from a genuine E_OUT_OF_VALIDATED_RANGE classification.  Legacy mode
# and backends without two-phase support do not use the range.
# ---------------------------------------------------------------------------

def _phase_probe_backend(property_validation, supports_two_phase=True):
    backend = type("PhaseProbeBackend", (), {})()
    backend.name = "probe"
    backend.supports_two_phase = supports_two_phase
    backend.property_validation = property_validation
    backend.pmin = 1.0e5
    backend.pmax = 8.0e6
    return backend


@pytest.mark.fast
def test_reference_backend_states_its_analytic_phase_fractions(tmp_path):
    document = json.loads(
        (ROOT / "cases" / "golden" / "g1_reference_constant.json").read_text(encoding="utf-8")
    )
    # The shipped golden g1 uses single_phase_legacy; an explicit
    # homogeneous_equilibrium reference run reaches the two-phase gate.
    assert document["two_phase_model"]["mode"] == gribs.TWO_PHASE_SINGLE_LEGACY
    document["two_phase_model"]["mode"] = gribs.TWO_PHASE_HOMOGENEOUS
    document["output"]["directory"] = str(tmp_path / "out")
    config = gribs.configuration_from_document(document).config
    backend = gribs.initialize_thermochemistry(config, tmp_path / "cache")

    assert backend.name == gribs.BACKEND_REFERENCE_CONSTANT
    assert config.two_phase_mode == gribs.TWO_PHASE_HOMOGENEOUS
    assert backend.supports_two_phase is True
    pv = backend.property_validation
    for pressure in (backend.pmin, backend.pmax, 0.5 * (backend.pmin + backend.pmax)):
        yg, yc = backend.phase_fractions(pressure)
        assert pv["gas_mass_fraction_min"] == yg == pv["gas_mass_fraction_max"] == 1.0
        assert pv["condensed_mass_fraction_min"] == yc == pv["condensed_mass_fraction_max"] == 0.0
    assert backend.metadata()["gas_mass_fraction"] == pv["gas_mass_fraction_min"]
    assert gribs.phase_fraction_validation_error(config, backend) is None
    assert gribs.model_applicability_error(config, backend) is None


@pytest.mark.fast
@pytest.mark.parametrize(
    "property_validation, expected_problem",
    [
        ({}, "missing:gas_mass_fraction_min"),
        ({"gas_mass_fraction_max": 1.0}, "missing:gas_mass_fraction_min"),
        ({"gas_mass_fraction_min": 0.9}, "missing:gas_mass_fraction_max"),
        (None, "property_validation_missing"),
        ({"gas_mass_fraction_min": float("nan"), "gas_mass_fraction_max": 1.0}, "non_finite:gas_mass_fraction_min"),
        ({"gas_mass_fraction_min": 0.9, "gas_mass_fraction_max": float("inf")}, "non_finite:gas_mass_fraction_max"),
        ({"gas_mass_fraction_min": "0.9", "gas_mass_fraction_max": 1.0}, "not_a_number:gas_mass_fraction_min"),
        ({"gas_mass_fraction_min": True, "gas_mass_fraction_max": 1.0}, "not_a_number:gas_mass_fraction_min"),
        ({"gas_mass_fraction_min": 0.0, "gas_mass_fraction_max": 1.0}, "invalid_range"),
        ({"gas_mass_fraction_min": -0.1, "gas_mass_fraction_max": 1.0}, "invalid_range"),
        ({"gas_mass_fraction_min": 0.95, "gas_mass_fraction_max": 0.9}, "invalid_range"),
        ({"gas_mass_fraction_min": 0.9, "gas_mass_fraction_max": 1.2}, "invalid_range"),
        ({"gas_mass_fraction_min": 10**400, "gas_mass_fraction_max": 1.0}, "not_representable:gas_mass_fraction_min"),
        ({"gas_mass_fraction_min": 0.9, "gas_mass_fraction_max": -(10**400)}, "not_representable:gas_mass_fraction_max"),
    ],
)
def test_invalid_phase_fraction_validation_is_a_distinct_structured_error(
    property_validation, expected_problem,
):
    config = gribs.Config()
    assert config.two_phase_mode == gribs.TWO_PHASE_HOMOGENEOUS
    backend = _phase_probe_backend(property_validation)

    error = gribs.model_applicability_error(config, backend)

    assert error is not None
    assert error["code"] == gribs.PHASE_FRACTION_VALIDATION_ERROR_CODE
    assert error["code"] != "E_OUT_OF_VALIDATED_RANGE"
    assert error["severity"] == "error"
    assert expected_problem in error["context"]["problems"]
    assert error["context"]["model_status"] is None
    json.dumps(error, allow_nan=False)  # reported values are strings; record is strict JSON
    assert gribs.phase_fraction_validation_error(config, backend) == error


@pytest.mark.fast
def test_valid_phase_fraction_range_is_still_classified_by_the_thresholds():
    config = gribs.Config()
    in_scope = _phase_probe_backend({"gas_mass_fraction_min": 0.7, "gas_mass_fraction_max": 1.0})
    assert gribs.model_applicability_error(config, in_scope) is None
    out_of_scope = _phase_probe_backend({"gas_mass_fraction_min": 0.5, "gas_mass_fraction_max": 0.9})
    error = gribs.model_applicability_error(config, out_of_scope)
    assert error["code"] == "E_OUT_OF_VALIDATED_RANGE"
    assert error["context"]["model_status"] == "out_of_scope"


@pytest.mark.fast
@pytest.mark.parametrize("property_validation", [{}, None, {"gas_mass_fraction_min": float("nan")}])
def test_legacy_and_inactive_two_phase_do_not_use_the_phase_fraction_range(property_validation):
    legacy = gribs.Config()
    legacy.two_phase_mode = gribs.TWO_PHASE_SINGLE_LEGACY
    backend = _phase_probe_backend(property_validation)
    assert gribs.phase_fraction_validation_error(legacy, backend) is None
    assert gribs.model_applicability_error(legacy, backend) is None

    homogeneous = gribs.Config()
    no_two_phase = _phase_probe_backend(property_validation, supports_two_phase=False)
    assert gribs.phase_fraction_validation_error(homogeneous, no_two_phase) is None
    assert gribs.model_applicability_error(homogeneous, no_two_phase) is None


@pytest.mark.fast
def test_main_reports_invalid_phase_fraction_validation_before_the_solver(
    tmp_path,
    monkeypatch,
    capsys,
):
    document = json.loads((ROOT / "gribs_config.json").read_text(encoding="utf-8"))
    output_dir = tmp_path / "must_not_be_created"
    document["output"]["directory"] = str(output_dir)
    config_path = tmp_path / "missing_phase_fractions.json"
    config_path.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")

    backend = _phase_probe_backend({"passed": True})
    monkeypatch.setattr(gribs, "initialize_thermochemistry", lambda _config, _cache_dir: backend)
    monkeypatch.setattr(gribs, "_print_backend_banner", lambda *_args: None)
    monkeypatch.setattr(
        gribs, "run_model",
        lambda *_args: pytest.fail("run with unusable phase-fraction validation reached the solver"),
    )

    exit_code = gribs.main(["--config", str(config_path)])
    captured = capsys.readouterr()

    assert exit_code == 1
    assert gribs.PHASE_FRACTION_VALIDATION_ERROR_CODE in captured.err
    assert "E_OUT_OF_VALIDATED_RANGE" not in captured.err
    record = json.loads(captured.err.strip().splitlines()[-1])
    assert record["code"] == gribs.PHASE_FRACTION_VALIDATION_ERROR_CODE
    assert record["context"]["simulation_started"] is False
    assert "missing:gas_mass_fraction_min" in record["context"]["problems"]
    assert not output_dir.exists()


def _direct_summarize(monkeypatch, property_validation, two_phase_mode=None, supports_two_phase=True, table=None):
    c = gribs.configuration_from_document(
        json.loads((ROOT / "gribs_config.json").read_text(encoding="utf-8"))
    ).config
    if two_phase_mode is not None:
        c.two_phase_mode = two_phase_mode
    backend = _phase_probe_backend(property_validation, supports_two_phase=supports_two_phase)
    backend.name = gribs.BACKEND_CEA_PYTHON
    backend.cache_schema = 4
    backend.mass_basis_semantics_version = 1
    if table is not None:
        backend.table = table
    monkeypatch.setattr(gribs, "thermochemistry", lambda: backend)
    monkeypatch.setattr(gribs, "choke_limit_pressure", lambda config, At0, Ae: None)
    monkeypatch.setattr(
        gribs, "nozzle_transition_diagnostics",
        lambda history, config, Ae: {
            "policy": "jump", "applicable": False, "detected": False, "warning_code": None,
            "interval_count": 0, "duration_s": 0.0, "duration_fraction": 0.0,
            "impulse_Ns": 0.0, "impulse_fraction": 0.0,
        },
    )
    monkeypatch.setattr(gribs, "nozzle_state", lambda p, At, Ae, config: {"mdot_total": 1.0})
    monkeypatch.setattr(
        gribs, "mass_energy_ledger",
        lambda res, h, config: {"schema_version": "a4-ledger-1", "mass": {"passed": True},
                                "energy": {"available": False, "closure_residual_available": False},
                                "passed": True},
    )
    m_prop = gribs.propellant_mass(c)
    t = gribs.np.array([0.0, 1.0])
    At0 = gribs.math.pi * c.R_t0**2
    history = {
        "t": t,
        "phase": gribs.np.array(["burn", "burn"], dtype=object),
        "p0": gribs.np.array([2.0e5, 3.0e5]),
        "F": gribs.np.array([0.0, 1.0]),
        "impulse": gribs.np.array([0.0, 1.0]),
        "m_out": gribs.np.array([0.0, m_prop]),
        "m_gen": gribs.np.array([0.0, m_prop]),
        "At": gribs.np.array([At0, At0]),
        "Kn": gribs.np.array([1.0, 1.0]),
        "Rt": gribs.np.array([c.R_t0, c.R_t0]),
        "m_total_eos": gribs.np.array([1.0, 1.0]),
        "m_total_balance": gribs.np.array([1.0, 1.0]),
        "Y_gas": gribs.np.array([1.0, 1.0]),
        "Y_condensed": gribs.np.array([0.0, 0.0]),
        "extrap": gribs.np.array([False, False]),
        "regime": gribs.np.array(["subsonic", "subsonic"], dtype=object),
    }
    result = {
        "m_total0": 1.0, "At0": At0, "Ae": At0, "stop_burn": "time-limit",
        "stop_blow": "not-run", "init_regime": "subsonic", "t_events": {},
    }
    return gribs.summarize(result, history, c)


@pytest.mark.fast
@pytest.mark.parametrize("property_validation", [
    {},
    {"gas_mass_fraction_min": float("nan"), "gas_mass_fraction_max": 1.0},
    {"gas_mass_fraction_min": 10**400, "gas_mass_fraction_max": 1.0},
])
def test_direct_summarize_rejects_unusable_phase_fraction_validation(monkeypatch, property_validation):
    with pytest.raises(gribs.PhaseFractionValidationError) as info:
        _direct_summarize(monkeypatch, property_validation)
    assert info.value.code == gribs.PHASE_FRACTION_VALIDATION_ERROR_CODE
    assert info.value.record["code"] == gribs.PHASE_FRACTION_VALIDATION_ERROR_CODE


@pytest.mark.fast
def test_direct_summarize_legacy_and_inactive_two_phase_keep_previous_behavior(monkeypatch):
    # P1 (K-02): legacy -> regression_comparison_only, inactive two-phase ->
    # in_scope_not_validated; both map back to the former validated_in_scope.
    legacy = _direct_summarize(monkeypatch, {}, two_phase_mode=gribs.TWO_PHASE_SINGLE_LEGACY)
    assert legacy["validation"]["model_status"] == "regression_comparison_only"
    assert legacy["validation"]["status_previous_tokens"]["model_status"] == "validated_in_scope"
    inactive = _direct_summarize(monkeypatch, {}, supports_two_phase=False)
    assert inactive["validation"]["model_status"] == "in_scope_not_validated"
    assert inactive["validation"]["status_previous_tokens"]["model_status"] == "validated_in_scope"
    valid = _direct_summarize(monkeypatch, {"gas_mass_fraction_min": 0.8, "gas_mass_fraction_max": 1.0})
    assert valid["validation"]["model_status"] == "approximate_with_limitations"


# Reported values must never make the phase-fraction check itself raise.
class _ReprRaises:
    def __repr__(self):
        raise RuntimeError("repr refuses")


class _FloatReprRaises(float):
    def __repr__(self):
        raise RuntimeError("float repr refuses")


@pytest.mark.fast
@pytest.mark.parametrize("raw, expected", [
    (0.95, "0.95"),
    ("0.9", "'0.9'"),
    (float("nan"), "nan"),
    (True, "True"),
    (10**400, repr(10**400)[:200]),
], ids=["float", "str", "nan", "bool", "int_1e400"])
def test_phase_fraction_report_keeps_normal_repr(raw, expected):
    error = gribs.model_applicability_error(
        gribs.Config(), _phase_probe_backend({"gas_mass_fraction_min": raw, "gas_mass_fraction_max": 0.9})
    )
    assert error["code"] == gribs.PHASE_FRACTION_VALIDATION_ERROR_CODE
    assert error["context"]["reported_values"]["gas_mass_fraction_min"] == expected


@pytest.mark.fast
@pytest.mark.parametrize("raw", [10**4300, 10**5000, -(10**5000)], ids=["1e4300", "1e5000", "-1e5000"])  # explicit ids: str() of these ints exceeds the digit limit
def test_huge_integer_beyond_str_digit_limit_is_a_structured_error(raw):
    error = gribs.model_applicability_error(
        gribs.Config(), _phase_probe_backend({"gas_mass_fraction_min": raw, "gas_mass_fraction_max": 1.0})
    )
    assert error["code"] == gribs.PHASE_FRACTION_VALIDATION_ERROR_CODE
    assert "not_representable:gas_mass_fraction_min" in error["context"]["problems"]
    text = error["context"]["reported_values"]["gas_mass_fraction_min"]
    try:
        expected = repr(raw)[:200]          # within the interpreter's digit limit
    except ValueError:
        expected = None                     # beyond it (default limit: 4300 digits)
    if expected is None:
        assert text.startswith(f"<int with {raw.bit_length()} bits")
    else:
        assert text == expected
    assert len(text) <= 200
    json.dumps(error, allow_nan=False)


@pytest.mark.fast
def test_lowered_str_digit_limit_is_restored_after_the_check():
    original = sys.get_int_max_str_digits()
    try:
        sys.set_int_max_str_digits(640)  # smallest allowed non-zero limit
        raw = 10**700                    # repr fails under this limit
        error = gribs.model_applicability_error(
            gribs.Config(), _phase_probe_backend({"gas_mass_fraction_min": raw, "gas_mass_fraction_max": 1.0})
        )
        assert error["code"] == gribs.PHASE_FRACTION_VALIDATION_ERROR_CODE
        assert error["context"]["reported_values"]["gas_mass_fraction_min"].startswith(
            f"<int with {raw.bit_length()} bits"
        )
    finally:
        sys.set_int_max_str_digits(original)
    assert sys.get_int_max_str_digits() == original


@pytest.mark.fast
def test_value_whose_repr_raises_is_a_structured_error():
    error = gribs.model_applicability_error(
        gribs.Config(), _phase_probe_backend({"gas_mass_fraction_min": _ReprRaises(), "gas_mass_fraction_max": 1.0})
    )
    assert error["code"] == gribs.PHASE_FRACTION_VALIDATION_ERROR_CODE
    assert "not_a_number:gas_mass_fraction_min" in error["context"]["problems"]
    assert error["context"]["reported_values"]["gas_mass_fraction_min"] == "<_ReprRaises value; repr failed (RuntimeError)>"

    # a numeric value whose repr raises: classified on its value, reported safely
    error = gribs.model_applicability_error(
        gribs.Config(), _phase_probe_backend({"gas_mass_fraction_min": _FloatReprRaises(float("inf")), "gas_mass_fraction_max": 1.0})
    )
    assert "non_finite:gas_mass_fraction_min" in error["context"]["problems"]
    assert error["context"]["reported_values"]["gas_mass_fraction_min"] == "<_FloatReprRaises value; repr failed (RuntimeError)>"
    json.dumps(error, allow_nan=False)


@pytest.mark.fast
def test_direct_summarize_range_fields_tolerate_incomplete_yg_when_two_phase_is_not_active(monkeypatch):
    """Legacy / two-phase-inactive runs do not validate the Yg range; a missing
    Yg maximum must not break the (P3b) table_Yg_range field, and the former
    validated_Yg_range value is kept."""
    table = {"pressure_Pa": [1.0e5, 8.0e6]}
    pv = {"gas_mass_fraction_min": 0.9}            # maximum missing
    legacy = _direct_summarize(monkeypatch, pv, two_phase_mode=gribs.TWO_PHASE_SINGLE_LEGACY, table=table)
    inactive = _direct_summarize(monkeypatch, pv, supports_two_phase=False, table=table)
    for summary in (legacy, inactive):
        v = summary["validation"]
        assert v["validated_Yg_range"] == [0.9, 1.0]
        assert v["table_Yg_range"] is None
        assert v["table_pressure_range_Pa"] == [1.0e5, 8.0e6]


# float() conversion of a reported Yg value: ordinary exceptions become a
# structured error; KeyboardInterrupt / SystemExit are not swallowed.
def _float_raising(exc_type):
    class _FloatRaises(float):
        def __float__(self):
            raise exc_type("from __float__")
    return _FloatRaises(0.9)


@pytest.mark.fast
@pytest.mark.parametrize("exc_type", [RuntimeError, ArithmeticError, LookupError], ids=["RuntimeError", "ArithmeticError", "LookupError"])
def test_float_conversion_error_of_a_yg_value_is_a_structured_error(exc_type):
    backend = _phase_probe_backend({"gas_mass_fraction_min": _float_raising(exc_type), "gas_mass_fraction_max": 1.0})
    error = gribs.model_applicability_error(gribs.Config(), backend)
    assert error["code"] == gribs.PHASE_FRACTION_VALIDATION_ERROR_CODE
    assert "not_representable:gas_mass_fraction_min" in error["context"]["problems"]
    json.dumps(error, allow_nan=False)


@pytest.mark.fast
@pytest.mark.parametrize("exc_type", [KeyboardInterrupt, SystemExit], ids=["KeyboardInterrupt", "SystemExit"])
def test_float_conversion_does_not_swallow_interrupts(exc_type):
    backend = _phase_probe_backend({"gas_mass_fraction_min": _float_raising(exc_type), "gas_mass_fraction_max": 1.0})
    with pytest.raises(exc_type):
        gribs.phase_fraction_validation_error(gribs.Config(), backend)
