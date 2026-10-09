"""A-5.1-3 independent red tests for the corrected CEA mass basis.

These tests state the A-5.1-4 acceptance contract before product correction.
They are intentionally isolated from normal pytest discovery.

The synthetic state is not a substitute for measured CEA evidence.  It is a
small structural fixture with a deliberately large condensed fraction so that
an extra factor of Yg is unambiguous.
"""

from __future__ import annotations

import importlib.util
import json
import math
import sys
from dataclasses import asdict
from pathlib import Path
from types import SimpleNamespace

import pytest


ROOT = Path(__file__).resolve().parents[2]
MODULE_FILES = sorted(ROOT.glob("GRIBS_v*.py"))

if len(MODULE_FILES) != 1:
    raise RuntimeError(
        "Expected exactly one GRIBS_v*.py file, found "
        f"{len(MODULE_FILES)}: {[path.name for path in MODULE_FILES]}"
    )

MODULE_PATH = MODULE_FILES[0]
SPEC = importlib.util.spec_from_file_location("gribs_a5_1_3_red", MODULE_PATH)

if SPEC is None or SPEC.loader is None:
    raise ImportError(f"Could not load GRIBS module from {MODULE_PATH}")

gribs = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = gribs
SPEC.loader.exec_module(gribs)


YG = 0.7
YC = 0.3
M_CEA_KG_KMOL = 25.0
M_GAS_KG_KMOL = YG * M_CEA_KG_KMOL
R_EFF_J_KGK = gribs.R_UNIVERSAL / M_CEA_KG_KMOL
R_GAS_J_KGK = gribs.R_UNIVERSAL / M_GAS_KG_KMOL
T_USED_K = 3000.0
GAMMA = 1.2
PRESSURE_PA = 5.0e6


class IndependentMassBasisBackend:
    """Minimal corrected-semantics backend owned by the red tests."""

    name = "a5_1_3_independent_mass_basis"
    uses_cea = False
    supports_two_phase = True
    pmin = 1.0e5
    pmax = 8.0e6

    def props(self, pressure_pa):
        del pressure_pa
        return R_EFF_J_KGK, T_USED_K, GAMMA

    def phase_fractions(self, pressure_pa):
        del pressure_pa
        return YG, YC

    def psi_and_derivative(self, pressure_pa):
        del pressure_pa
        return R_EFF_J_KGK * T_USED_K, 0.0

    def theta_and_derivative(self, pressure_pa):
        del pressure_pa
        return R_EFF_J_KGK * T_USED_K, 0.0

    def is_extrapolated(self, pressure_pa):
        return not self.pmin <= pressure_pa <= self.pmax

    @property
    def outside_range_policy(self):
        return "extrapolate"

    def metadata(self):
        return {
            "backend": self.name,
            "synthetic_structural_fixture": True,
            "molecular_weight_definition": (
                "M_CEA = 1 / EqSolution.n on total equilibrium-mixture mass basis"
            ),
            "effective_gas_constant_definition": "R_eff = Ru / M_CEA",
        }


@pytest.fixture
def corrected_backend(monkeypatch):
    backend = IndependentMassBasisBackend()
    monkeypatch.setattr(gribs, "_THERMO", backend)
    return backend


@pytest.fixture
def product_psi_backend(monkeypatch):
    """Load an independent table through the product Psi construction path."""

    backend = gribs.PressureTableBackend.__new__(
        gribs.PressureTableBackend
    )
    backend.name = "a5_1_3_product_psi_probe"
    backend.supports_two_phase = True
    backend.c = gribs.Config(
        two_phase_mode=gribs.TWO_PHASE_HOMOGENEOUS,
        eta_T0=1.0,
        property_policy="extrapolate",
    )

    pressures = [1.0e5, 2.0e5, 4.0e5, 8.0e5]
    table = {
        "pressure_Pa": pressures,
        "temperature_K": [T_USED_K] * len(pressures),
        "gamma_s": [GAMMA] * len(pressures),
        "equilibrium_molecular_weight_kg_kmol": [
            M_CEA_KG_KMOL
        ] * len(pressures),
        "effective_gas_constant_J_kgK": [
            R_EFF_J_KGK
        ] * len(pressures),
        "gas_mass_fraction": [YG] * len(pressures),
        "condensed_mass_fraction": [YC] * len(pressures),
    }

    backend.table = table
    backend._load_table(table)
    monkeypatch.setattr(gribs, "_THERMO", backend)
    return backend


@pytest.fixture
def homogeneous_config():
    return gribs.Config(
        two_phase_mode=gribs.TWO_PHASE_HOMOGENEOUS,
        p_a=1.0e5,
        Cd=0.97,
        eta_thrust=0.93,
        use_separation=False,
        eps_nozzle=1.0,
    )


def choked_mass_flow(config, pressure_pa, throat_area):
    critical_factor = (
        2.0 / (GAMMA + 1.0)
    ) ** ((GAMMA + 1.0) / (2.0 * (GAMMA - 1.0)))
    return (
        config.Cd
        * throat_area
        * pressure_pa
        * math.sqrt(GAMMA / (R_EFF_J_KGK * T_USED_K))
        * critical_factor
    )


def test_cea_molecular_weight_and_gas_constant_relations():
    """CEA M is total-mixture based; pure-gas M is smaller by Yg."""

    assert M_CEA_KG_KMOL == pytest.approx(M_GAS_KG_KMOL / YG)
    assert R_EFF_J_KGK == pytest.approx(YG * R_GAS_J_KGK)


def test_density_identity_uses_effective_gas_constant():
    """The homogeneous-mixture EOS gives p/(rho*T) = R_eff."""

    rho_total = PRESSURE_PA / (R_EFF_J_KGK * T_USED_K)

    assert PRESSURE_PA / (rho_total * T_USED_K) == pytest.approx(
        R_EFF_J_KGK
    )


def test_chamber_psi_is_reff_times_temperature(
    product_psi_backend,
    homogeneous_config,
):
    """Product table loading must construct Psi without a second Yg."""

    actual, derivative = gribs.psi_and_deriv(
        PRESSURE_PA,
        homogeneous_config,
    )
    expected = R_EFF_J_KGK * T_USED_K

    assert derivative == pytest.approx(0.0, abs=1.0e-12)
    assert actual == pytest.approx(expected), (
        "A-5.1-3 expected product Psi = R_eff*T_used; "
        f"actual={actual:.16g}, expected={expected:.16g}, "
        f"actual/expected={actual / expected:.16g}, Yg={YG:.16g}"
    )


def test_chamber_psi_rejects_legacy_extra_yg(
    product_psi_backend,
    homogeneous_config,
):
    """Product table loading must reject Yg*R_eff*T when Yg != 1."""

    actual, _ = gribs.psi_and_deriv(PRESSURE_PA, homogeneous_config)
    legacy = YG * R_EFF_J_KGK * T_USED_K

    assert actual != pytest.approx(legacy), (
        "A-5.1-3 rejects the product's legacy extra Yg in Psi: "
        f"actual={actual:.16g}, legacy={legacy:.16g}, Yg={YG:.16g}"
    )


def test_nozzle_formula_returns_total_flow_without_extra_yg(
    corrected_backend,
    homogeneous_config,
):
    """The isentropic formula directly returns total-mixture mass flow."""

    throat_area = 1.0e-5
    state = gribs.nozzle_state(
        PRESSURE_PA,
        throat_area,
        throat_area,
        homogeneous_config,
    )
    expected_total = choked_mass_flow(
        homogeneous_config,
        PRESSURE_PA,
        throat_area,
    )

    assert state["mdot_total"] == pytest.approx(expected_total), (
        "A-5.1-3 forbids an additional /Yg in nozzle flow: "
        f"actual_total={state['mdot_total']:.16g}, "
        f"expected_total={expected_total:.16g}, "
        f"ratio={state['mdot_total'] / expected_total:.16g}, "
        f"Yg={YG:.16g}"
    )
    assert state["mdot"] == pytest.approx(expected_total)


def test_nozzle_phase_decomposition_is_conservative(
    corrected_backend,
    homogeneous_config,
):
    """Phase flows decompose, rather than rescale, mdot_total."""

    throat_area = 1.0e-5
    state = gribs.nozzle_state(
        PRESSURE_PA,
        throat_area,
        throat_area,
        homogeneous_config,
    )

    assert state["mdot_gas"] == pytest.approx(YG * state["mdot_total"])
    assert state["mdot_condensed"] == pytest.approx(
        YC * state["mdot_total"]
    )
    assert (
        state["mdot_gas"] + state["mdot_condensed"]
    ) == pytest.approx(state["mdot_total"])


def test_momentum_thrust_uses_total_flow_once(
    corrected_backend,
    homogeneous_config,
):
    """Momentum thrust is eta_F*mdot_total*ve plus pressure thrust."""

    throat_area = 1.0e-5
    state = gribs.nozzle_state(
        PRESSURE_PA,
        throat_area,
        throat_area,
        homogeneous_config,
    )
    expected = (
        homogeneous_config.eta_thrust
        * state["mdot_total"]
        * state["ve"]
        + (state["pe"] - homogeneous_config.p_a) * state["A_eff"]
    )

    assert state["F"] == pytest.approx(expected)


def test_single_phase_limit_is_unchanged(monkeypatch, homogeneous_config):
    """At Yg=1, corrected and legacy mass-basis formulas coincide."""

    backend = IndependentMassBasisBackend()
    monkeypatch.setattr(backend, "phase_fractions", lambda pressure_pa: (1.0, 0.0))
    monkeypatch.setattr(gribs, "_THERMO", backend)

    throat_area = 1.0e-5
    state = gribs.nozzle_state(
        PRESSURE_PA,
        throat_area,
        throat_area,
        homogeneous_config,
    )
    expected_total = choked_mass_flow(
        homogeneous_config,
        PRESSURE_PA,
        throat_area,
    )

    assert state["mdot_total"] == pytest.approx(expected_total)
    assert state["mdot_gas"] == pytest.approx(expected_total)
    assert state["mdot_condensed"] == pytest.approx(0.0)


def test_corrected_and_legacy_total_flow_relation(
    corrected_backend,
    homogeneous_config,
):
    """Corrected total divided by legacy total is 1/sqrt(Yg)."""

    throat_area = 1.0e-5
    corrected = gribs.nozzle_state(
        PRESSURE_PA,
        throat_area,
        throat_area,
        homogeneous_config,
    )["mdot_total"]

    legacy_config = gribs.Config(
        **{
            **asdict(homogeneous_config),
            "two_phase_mode": gribs.TWO_PHASE_SINGLE_LEGACY,
        }
    )

    monkey_backend = SimpleNamespace(
        props=lambda pressure_pa: (R_GAS_J_KGK, T_USED_K, GAMMA),
        phase_fractions=lambda pressure_pa: (1.0, 0.0),
        psi_and_derivative=lambda pressure_pa: (
            R_GAS_J_KGK * T_USED_K,
            0.0,
        ),
        theta_and_derivative=lambda pressure_pa: (
            R_GAS_J_KGK * T_USED_K,
            0.0,
        ),
        is_extrapolated=lambda pressure_pa: False,
    )

    saved_backend = gribs._THERMO
    try:
        gribs._THERMO = monkey_backend
        legacy = gribs.nozzle_state(
            PRESSURE_PA,
            throat_area,
            throat_area,
            legacy_config,
        )["mdot_total"]
    finally:
        gribs._THERMO = saved_backend

    assert corrected / legacy == pytest.approx(1.0 / math.sqrt(YG)), (
        "A-5.1-3 corrected/legacy total-flow relation failed: "
        f"actual_ratio={corrected / legacy:.16g}, "
        f"expected_ratio={1.0 / math.sqrt(YG):.16g}, Yg={YG:.16g}"
    )


def test_cache_schema_3_is_rejected(tmp_path):
    """A complete legacy schema-3 cache must be rejected by schema identity."""

    backend = gribs.PressureTableBackend.__new__(
        gribs.PressureTableBackend
    )
    backend.name = "schema_rejection_probe"
    backend.key = "a5-1-3-schema-probe"
    backend.c = SimpleNamespace(
        two_phase_mode=gribs.TWO_PHASE_HOMOGENEOUS,
        cea_pressure_points=4,
    )
    backend.cache_path = tmp_path / "schema3.json"

    table = {
        "schema": 3,
        "backend": backend.name,
        "cache_key": backend.key,
        "pressure_Pa": [1.0e5, 2.0e5, 4.0e5, 8.0e5],
        "temperature_K": [T_USED_K] * 4,
        "gamma_s": [GAMMA] * 4,
        "equilibrium_molecular_weight_kg_kmol": [M_CEA_KG_KMOL] * 4,
        "effective_gas_constant_J_kgK": [R_EFF_J_KGK] * 4,
        "gas_mass_fraction": [YG] * 4,
        "condensed_mass_fraction": [YC] * 4,
        "psi_J_kg": [R_EFF_J_KGK * T_USED_K] * 4,
    }
    backend.cache_path.write_text(
        json.dumps(table, indent=2) + "\n",
        encoding="utf-8",
    )

    with pytest.raises(
        gribs.ThermochemistryCacheError,
        match=r"schema 3|mass.basis|incompatible",
    ):
        backend._read_cache()


def test_metadata_does_not_claim_eqsolution_m_is_pure_gas_basis():
    """CEA metadata must describe EqSolution.M on total-mixture mass basis."""

    source = MODULE_PATH.read_text(encoding="utf-8")
    forbidden_claims = (
        "EqSolution.M (gas-phase, kg/kmol)",
        "EqSolution.M - the GAS-PHASE molecular weight",
        "gas_phase_M (gas-phase M;",
    )

    present = [claim for claim in forbidden_claims if claim in source]

    assert present == [], (
        "A-5.1-3 found obsolete pure-gas-basis claims for EqSolution.M: "
        f"{present}"
    )


@pytest.mark.xfail(
    reason=(
        "Public v0.4.3-alpha default gas_phase_M comparison artifact "
        "has not yet been obtained; this is an external comparison gap, "
        "not evidence for or against the corrected equations."
    ),
    strict=False,
)
def test_public_v043_default_gas_phase_m_comparison_is_available():
    """Reserve the external public-v0.4.3 comparison acceptance point."""

    artifact = (
        ROOT
        / "verification"
        / "mass_basis"
        / "public_v0.4.3-alpha_default_gas_phase_M.json"
    )

    assert artifact.is_file(), (
        "External comparison artifact is not present: "
        f"{artifact.relative_to(ROOT)}"
    )
