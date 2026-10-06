import importlib.util
import json
import re
import subprocess
import sys
from dataclasses import asdict
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
    defaults = asdict(gribs.Config())
    assert set(supplied).issubset(defaults)


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
    assert mass["propellant"]["tolerance"] == pytest.approx(1.0e-9)
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
    propellant_residual = 2.0e-9 * initial_propellant

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
    assert mass["propellant"]["relative_residual"] > 1.0e-9
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

    class Backend:
        name = "synthetic"
        pmin = 1.0e5
        pmax = 8.0e6
        supports_two_phase = True
        property_validation = {}

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
        "m_out": gribs.np.array([0.0, 0.5]),
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
        "At0": At0,
        "Ae": At0,
        "stop_burn": "time-limit",
        "stop_blow": "not-run",
        "init_regime": "subsonic",
        "t_events": {},
    }

    summary = gribs.summarize(result, history, c)

    assert summary["ledger"] is sentinel_ledger
    assert summary["propellant_mass_balance_error"] == pytest.approx(0.0)
    assert summary["total_mass_consistency_error"] == pytest.approx(0.0)


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
        "warnings": [],
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
    }

    text = gribs.summary_text(summary, c)

    assert "[ mass ledger ]" in text
    assert "propellant relative residual" in text
    assert "mass-ledger result" in text
    assert "[ energy ledger ]" in text
    assert "closure residual available" in text
    assert "not a complete energy-conservation law" in text
    assert "legacy propellant balance metric" in text


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
