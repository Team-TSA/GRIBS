"""Post-processing (A-6 unit 6). Moved verbatim from section 8 of the entry: the sampling of the solution on the
phase grid (sample), the nozzle-transition diagnostics, the post-run consistency report, the model and semantic
status, the phase-fraction validation and model applicability checks, the run warnings, summarize, the run
verdict and exit codes, and the mass / energy ledger.

Tests patch _phase_grid, sample, summarize, mass_energy_ledger, model_applicability_error and
nozzle_transition_diagnostics here, e.g. "gribs.simulation.postprocess.mass_energy_ledger"; summarize and
sample read them as module globals of this module, and the entry re-exports every name, so a patch on the
entry would not reach them. Names of other package modules that tests patch are read as <module>.<name>
(ADR decision 1). The state indices and igniter_mdot come from gribs.simulation.integrator (U-57).
"""
from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import brentq

import gribs.chamber.equilibrium
import gribs.nozzle
import gribs.propellant.burning
import gribs.thermo
from gribs.config import Config, geometry, propellant_mass, web_thickness
from gribs.constants import (
    BACKEND_CEA_LEGACY_EXECUTABLE,
    BACKEND_CEA_PYTHON,
    G0,
    HIGH_CONDENSATE_G4B_UNVALIDATED_YG_THRESHOLD,
    MODEL_STATUS_OUT_OF_SCOPE_YG_THRESHOLD,
    MODEL_STATUS_VALIDATED_YG_THRESHOLD,
    POST_RUN_MASS_TOLERANCE,
    TWO_PHASE_HOMOGENEOUS,
    TWO_PHASE_SINGLE_LEGACY,
    _PHASE_FRACTION_ROUNDING_TOL,
    _TRAPZ,
)
from gribs.errors import ThermochemistryError
from gribs.nozzle import choke_margin, mach_from_area
from gribs.propellant.burning import base_burn_rate
from gribs.simulation.integrator import IIMP, IMG, IMIGN, IMO, IP, IRT, IX, igniter_mdot
from gribs.thermo import ThermochemistryBackend
from gribs.util import _report_value_text


def _phase_grid(sol, n_lin: int, log_head: bool) -> np.ndarray:
    t0, t1 = float(sol.t[0]), float(sol.t[-1])
    if t1 <= t0:
        return np.array([t0])
    g = [np.linspace(t0, t1, n_lin)]
    span = t1 - t0
    if log_head:
        g.append(t0 + np.geomspace(max(1e-7, span * 1e-7), span, 900))
    g.append(np.asarray(sol.t, dtype=float))
    for te in sol.t_events:
        if len(te):
            g.append(np.asarray(te, dtype=float))
    t = np.unique(np.concatenate(g))
    return t[(t >= t0) & (t <= t1)]


def sample(res: dict, c: Config) -> dict:
    """Dense post-processing of the ODE solution.

    Column meanings (v0.5.0-alpha two-phase model B):
      mdot_gen_total        generated TOTAL product mass flow [kg/s]
      mdot_out_total/gas/condensed  nozzle outflow decomposition [kg/s]
      m_gen, m_out          integrated TOTAL product masses [kg] (ODE states)
      m_igniter           integrated igniter-source mass [kg] (ODE state)
      m_total_eos           TOTAL product mass from the EOS: p0*Vg/Psi(p0) [kg]
                            (replaces the v0.4.3 column m_gas_eos, whose meaning
                            changed - the old name is deliberately NOT kept)
      m_total_balance       mt0 + m_gen - m_out [kg] (replaces m_gas_bal)
      m_gas_equilibrium     Yg(p0)*m_total_eos [kg] (gas part of m_total_eos)
      m_condensed           Yc(p0)*m_total_eos [kg] (condensed part)
      Y_gas, Y_condensed    equilibrium phase mass fractions [-]
      Psi                   R_eff*T0 [J/kg]
    """
    Ae = res["Ae"]
    xw = web_thickness(c)
    segs = [(res["sol_burn"], "burn", 4000, True)]
    if res["sol_blow"] is not None:
        segs.append((res["sol_blow"], "blowdown", 1600, True))

    rec = {k: [] for k in
           ("t", "phase", "p0", "x", "Rt", "At", "Kn", "Ri", "Lp", "Ab", "Vg",
            "r", "mdot_gen_total", "mdot_out_total", "mdot_out_gas",
            "mdot_out_condensed", "F", "Me", "pe", "ve", "regime",
            "R", "T0", "gamma", "Y_gas", "Y_condensed", "Psi",
            "m_out", "impulse", "m_gen", "m_igniter",
            "m_total_eos", "m_total_balance",
            "m_gas_equilibrium", "m_condensed", "extrap", "cstar", "CF")}

    first = True
    for sol, phase, n_lin, head in segs:
        tg = _phase_grid(sol, n_lin, head)
        if not first:
            tg = tg[1:]
        first = False
        ys = sol.sol(tg)
        for k, t in enumerate(tg):
            p0 = max(float(ys[IP, k]), 1.0)
            x = min(max(float(ys[IX, k]), 0.0), xw)
            Rt = float(ys[IRT, k])
            At = math.pi * Rt * Rt
            Ri, Lp, Ab, Vg = geometry(x, c)
            R, T0, g = gribs.thermo.gas_props(p0, c)
            Yg, Yc = gribs.thermo.phase_fractions(p0, c)
            psi, _ = gribs.thermo.psi_and_deriv(p0, c)
            nz = gribs.nozzle.nozzle_state(p0, At, Ae, c)
            if phase == "burn":
                r = gribs.propellant.burning.burn_rate(p0, c, Ab, Ri)
            else:
                r, Ab = 0.0, 0.0
            # Match make_rhs(): igniter injection is part of the total
            # generated-product rate in every integrated phase.
            m_gen_rate = c.rho_p * r * Ab + igniter_mdot(t, c)
            m_total_eos = p0 * Vg / psi
            m_total_bal = res["m_total0"] + float(ys[IMG, k]) - float(ys[IMO, k])

            rec["t"].append(float(t));          rec["phase"].append(phase)
            rec["p0"].append(p0);               rec["x"].append(x)
            rec["Rt"].append(Rt);               rec["At"].append(At)
            rec["Kn"].append(Ab / At);          rec["Ri"].append(Ri)
            rec["Lp"].append(Lp);               rec["Ab"].append(Ab)
            rec["Vg"].append(Vg);               rec["r"].append(r)
            rec["mdot_gen_total"].append(m_gen_rate)
            rec["mdot_out_total"].append(nz["mdot_total"])
            rec["mdot_out_gas"].append(nz["mdot_gas"])
            rec["mdot_out_condensed"].append(nz["mdot_condensed"])
            rec["F"].append(nz["F"]);           rec["Me"].append(nz["Me"])
            rec["pe"].append(nz["pe"]);         rec["ve"].append(nz["ve"])
            rec["regime"].append(nz["regime"]); rec["R"].append(R)
            rec["T0"].append(T0);               rec["gamma"].append(g)
            rec["Y_gas"].append(Yg);            rec["Y_condensed"].append(Yc)
            rec["Psi"].append(psi)
            rec["m_out"].append(float(ys[IMO, k]))
            rec["impulse"].append(float(ys[IIMP, k]))
            rec["m_gen"].append(float(ys[IMG, k]))
            rec["m_igniter"].append(float(ys[IMIGN, k]))
            rec["m_total_eos"].append(m_total_eos)
            rec["m_total_balance"].append(m_total_bal)
            rec["m_gas_equilibrium"].append(Yg * m_total_eos)
            rec["m_condensed"].append(Yc * m_total_eos)
            rec["extrap"].append(gribs.thermo.is_extrapolated(p0, c))
            rec["cstar"].append(p0 * At / nz["mdot_total"]
                                if nz["mdot_total"] > 0 else np.nan)
            rec["CF"].append(nz["F"] / (p0 * At) if p0 * At > 0 else np.nan)

    out = {k: (np.asarray(v) if k not in ("phase", "regime") else np.asarray(v, dtype=object))
           for k, v in rec.items()}

    # quasi-steady equilibrium pressure (diagnostic, burning phase only)
    # S-7: default OFF for performance (each point does a pressure sweep + brentq)
    t = out["t"]
    p_eq = np.full_like(t, np.nan)
    if c.equilibrium_pressure_diagnostic == "sampled":
        burn_idx = np.nonzero(out["phase"] == "burn")[0]
        if burn_idx.size:
            pick = burn_idx[np.unique(np.linspace(0, burn_idx.size - 1,
                                                  min(200, burn_idx.size)).astype(int))]
            vals = np.array([gribs.chamber.equilibrium.equilibrium_pressure(out["Ab"][i], out["Ri"][i],
                                                  out["At"][i], Ae, c) or np.nan
                             for i in pick])
            ok = np.isfinite(vals)
            if ok.sum() >= 2:
                p_eq[burn_idx] = np.interp(t[burn_idx], t[pick][ok], vals[ok],
                                           left=np.nan, right=np.nan)
    out["p_eq"] = p_eq
    return out


def _masked_trapezoid(
    y: np.ndarray,
    x: np.ndarray,
    mask: np.ndarray,
) -> float:
    """Integrate only intervals whose two endpoints are inside the mask."""
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    mask = np.asarray(mask, dtype=bool)

    if y.ndim != 1 or x.ndim != 1 or mask.ndim != 1:
        raise ValueError("y, x, and mask must be one-dimensional.")
    if not (y.size == x.size == mask.size):
        raise ValueError("y, x, and mask must have the same length.")
    if x.size < 2:
        return 0.0

    valid_intervals = mask[:-1] & mask[1:]
    dx = np.diff(x)
    interval_areas = 0.5 * (y[:-1] + y[1:]) * dx

    return float(np.sum(interval_areas[valid_intervals]))


def nozzle_transition_diagnostics(
    h: dict,
    c: Config,
    Ae: float,
    *,
    lower_pressure_Pa: Optional[float] = None,
    upper_pressure_Pa: Optional[float] = None,
) -> dict:
    """Quantify time and impulse spent in the jump-model transition band.

    Explicit pressure bounds are primarily intended for isolated tests. During
    normal calculations, the lower boundary is the choking threshold and the
    upper boundary is the onset of the configured flow-separation criterion.
    """
    result = {
        "policy": c.transition_policy,
        "applicable": False,
        "detected": False,
        "warning_code": None,
        "interval_count": 0,
        "duration_s": 0.0,
        "duration_fraction": 0.0,
        "impulse_Ns": 0.0,
        "impulse_fraction": 0.0,
    }

    if c.transition_policy != "jump":
        return result
    if c.eps_nozzle <= 1.0 + 1.0e-12:
        return result
    if not c.use_separation:
        return result

    t = np.asarray(h["t"], dtype=float)
    p0 = np.asarray(h["p0"], dtype=float)
    At = np.asarray(h["At"], dtype=float)
    thrust = np.asarray(h["F"], dtype=float)

    if not (t.size == p0.size == At.size == thrust.size):
        raise ValueError(
            "Transition diagnostic histories t, p0, At, and F "
            "must have the same length."
        )

    result["applicable"] = True

    if t.size < 2:
        return result

    if lower_pressure_Pa is not None or upper_pressure_Pa is not None:
        if lower_pressure_Pa is None or upper_pressure_Pa is None:
            raise ValueError(
                "lower_pressure_Pa and upper_pressure_Pa must be supplied together."
            )
        if not lower_pressure_Pa < upper_pressure_Pa:
            raise ValueError(
                "Require lower_pressure_Pa < upper_pressure_Pa."
            )

        transition_mask = (
            (p0 > float(lower_pressure_Pa))
            & (p0 < float(upper_pressure_Pa))
        )
    else:
        transition_mask = np.zeros(t.size, dtype=bool)

        for i in range(t.size):
            if not np.isfinite(p0[i]) or not np.isfinite(At[i]):
                continue
            if p0[i] <= c.p_a or At[i] <= 0.0:
                continue

            try:
                if choke_margin(p0[i], At[i], Ae, c) < 0.0:
                    continue

                _, _, gamma = gribs.thermo.gas_props(p0[i], c)
                eps = Ae / At[i]
                if eps <= 1.0 + 1.0e-12:
                    continue

                exit_mach = mach_from_area(eps, gamma, supersonic=True)
                exit_pressure = p0[i] * (
                    1.0 + 0.5 * (gamma - 1.0) * exit_mach * exit_mach
                ) ** (-gamma / (gamma - 1.0))

                transition_mask[i] = (
                    exit_pressure < c.sep_ratio * c.p_a
                )
            except (ValueError, RuntimeError, ThermochemistryError):
                continue

    starts = transition_mask & np.concatenate(
        (np.array([True]), ~transition_mask[:-1])
    )
    interval_count = int(np.count_nonzero(starts))

    duration = _masked_trapezoid(
        np.ones_like(t),
        t,
        transition_mask,
    )
    transition_impulse = _masked_trapezoid(
        thrust,
        t,
        transition_mask,
    )

    total_duration = float(t[-1] - t[0])
    total_impulse = float(_TRAPZ(thrust, t))

    detected = bool(interval_count > 0 and duration > 0.0)

    result.update(
        detected=detected,
        warning_code="W_NOZZLE_TRANSITION" if detected else None,
        interval_count=interval_count,
        duration_s=duration,
        duration_fraction=(
            duration / total_duration if total_duration > 0.0 else 0.0
        ),
        impulse_Ns=transition_impulse,
        impulse_fraction=(
            transition_impulse / total_impulse
            if abs(total_impulse) > 0.0
            else 0.0
        ),
    )
    return result


def post_run_consistency_report(res: dict, h: dict) -> dict:
    """Return one authoritative post-simulation mass-consistency report.

    The invariant EOS-history metric uses one run-wide reference mass and
    therefore does not amplify a fixed absolute error as the chamber
    inventory approaches zero during blowdown.

    The former instantaneous-mass-normalized metric remains available as a
    deprecated diagnostic for historical comparison, but it is not used as
    the post-run gate.
    """
    m_total0 = float(res["m_total0"])
    m_gen = np.asarray(h["m_gen"], dtype=float)
    m_out = np.asarray(h["m_out"], dtype=float)
    m_eos = np.asarray(h["m_total_eos"], dtype=float)
    m_balance = np.asarray(h["m_total_balance"], dtype=float)

    if m_gen.size == 0:
        raise ValueError("Post-run consistency requires a non-empty history.")

    expected_shape = m_eos.shape
    for name, values in (
        ("m_gen", m_gen),
        ("m_out", m_out),
        ("m_total_balance", m_balance),
    ):
        if values.shape != expected_shape:
            raise ValueError(
                f"Post-run history array {name} must match m_total_eos."
            )

    lhs = m_total0 + float(m_gen[-1])
    rhs = float(m_out[-1]) + float(m_eos[-1])
    final_relative_error = abs(lhs - rhs) / max(abs(lhs), 1.0e-300)

    absolute_difference = np.abs(m_balance - m_eos)

    legacy_denominator = np.maximum(np.abs(m_eos), 1.0e-15)
    legacy_relative_error = float(
        np.max(absolute_difference / legacy_denominator)
    )

    reference_mass = max(
        float(np.max(np.abs(m_balance))),
        1.0e-30,
    )
    invariant_relative_error = float(
        np.max(absolute_difference) / reference_mass
    )

    final_passed = bool(
        final_relative_error < POST_RUN_MASS_TOLERANCE
    )
    eos_history_passed = bool(
        invariant_relative_error < POST_RUN_MASS_TOLERANCE
    )

    return {
        "schema_version": "a5-1-2d-post-run-1",
        "tolerance": POST_RUN_MASS_TOLERANCE,
        "final_total_mass": {
            "lhs_initial_plus_generated_kg": lhs,
            "rhs_discharged_plus_remaining_kg": rhs,
            "relative_error": final_relative_error,
            "passed": final_passed,
        },
        "eos_history": {
            "metric_version": "run_wide_balance_mass_reference_v1",
            "reference_mass_kg": reference_mass,
            "invariant_relative_error": invariant_relative_error,
            "legacy_instantaneous_mass_relative_error": (
                legacy_relative_error
            ),
            "legacy_metric_deprecated": True,
            "passed": eos_history_passed,
        },
        "passed": bool(final_passed and eos_history_passed),
    }


def model_status_from_yg(
    two_phase_active: bool,
    two_phase_mode: str,
    gas_mass_fraction_min: float,
) -> str:
    """Classify the modeled Yg range using the A-5.1-6 scope thresholds.

    The thresholds are code constants that are not yet linked to external or
    experimental validation evidence (A-5.1 G4a/G4b open), so no branch claims
    'validated'.  ``single_phase_legacy`` is a regression-comparison mode.
    Former token 'validated_in_scope' -> MODEL_STATUS_PREVIOUS_TOKENS.
    """
    if two_phase_mode == TWO_PHASE_SINGLE_LEGACY:
        return "regression_comparison_only"
    if not two_phase_active:
        return "in_scope_not_validated"
    if gas_mass_fraction_min <= MODEL_STATUS_OUT_OF_SCOPE_YG_THRESHOLD:
        return "out_of_scope"
    if gas_mass_fraction_min <= MODEL_STATUS_VALIDATED_YG_THRESHOLD:
        return "approximate_with_limitations"
    return "in_scope_not_validated"


def _finite_number_or_none(raw: Any) -> Optional[float]:
    """float(raw) for a finite real number (not bool); None otherwise. Never raises."""
    if isinstance(raw, bool) or not isinstance(raw, (int, float, np.integer, np.floating)):
        return None
    try:
        value = float(raw)
    except Exception:  # noqa: BLE001 - e.g. OverflowError, or a __float__ that raises
        return None
    return value if math.isfinite(value) else None


def _validation_range_fields(backend: Any, pv_table: Any, has_table: bool) -> Dict[str, Any]:
    """validation.* range fields of the run summary (cleanup P3b).

    ``validated_pressure_range_Pa`` / ``validated_Yg_range`` keep their former
    values and expressions for existing readers.  ``table_*`` are filled only
    for a backend with a chamber-property table and a usable (both keys
    present, finite) Yg range; otherwise None, so a legacy / two-phase-inactive
    backend without a complete Yg range never fails here.
    ``validated_range_basis`` states what the values actually are.
    """
    yg_lo = _finite_number_or_none(pv_table.get("gas_mass_fraction_min")) if isinstance(pv_table, dict) else None
    yg_hi = _finite_number_or_none(pv_table.get("gas_mass_fraction_max")) if isinstance(pv_table, dict) else None
    yg_complete = yg_lo is not None and yg_hi is not None
    if has_table:
        basis = ("validated_pressure_range_Pa / validated_Yg_range are the "
                 "chamber-property table range, not an externally validated "
                 "range")
    elif yg_complete:
        basis = ("no chamber-property table: validated_pressure_range_Pa is the "
                 "configured pressure range and validated_Yg_range is the "
                 f"backend's stated analytic phase-fraction range [{yg_lo:g}, {yg_hi:g}] "
                 "(e.g. reference_constant: a single gas phase); neither is a "
                 "measured or validated range")
    else:
        basis = ("no chamber-property table: validated_pressure_range_Pa is the "
                 "configured pressure range and validated_Yg_range uses the "
                 "defaults 0 / 1 for the entries the backend does not state; "
                 "neither is a measured or validated range")
    return {
        "validated_pressure_range_Pa": [float(backend.pmin), float(backend.pmax)],
        "validated_Yg_range": [
            float(pv_table.get("gas_mass_fraction_min", 0.0)),
            float(pv_table.get("gas_mass_fraction_max", 1.0)),
        ],
        # The two 'validated_*' keys above are kept for existing readers;
        # their values are property-table (or stated analytic) ranges.
        "table_pressure_range_Pa": (
            [float(backend.pmin), float(backend.pmax)] if has_table else None),
        "table_Yg_range": [yg_lo, yg_hi] if (has_table and yg_complete) else None,
        "validated_range_basis": basis,
    }


def model_status_basis(two_phase_active: bool, two_phase_mode: str) -> str:
    """Cleanup P3c: explanation matching the branch taken in
    model_status_from_yg (thresholds unchanged)."""
    if two_phase_mode == TWO_PHASE_SINGLE_LEGACY:
        return ("two_phase_mode single_phase_legacy: regression-comparison "
                "mode; Yg thresholds are not applied")
    if not two_phase_active:
        return ("two-phase model inactive (the backend has no two-phase "
                "support); Yg thresholds are not applied; not linked to "
                "external or experimental validation evidence")
    return ("Yg_min threshold classification (code constants 0.99 / 0.65); "
            "not linked to external or experimental validation evidence")


# Cleanup 2026-11-01 (K-02/K-03): former status tokens, kept for readers that
# still compare against them.  Thresholds and classification are unchanged.
MODEL_STATUS_PREVIOUS_TOKENS = {
    "regression_comparison_only": "validated_in_scope",
    "in_scope_not_validated": "validated_in_scope",
    "approximate_with_limitations": "approximate_with_limitations",
    "out_of_scope": "out_of_scope",
}
# For not_applicable the former token depended on the backend (reference
# constant: 'fail'; synthetic table: 'pass'); the per-run value is reported
# in summary validation.status_previous_tokens.semantic_status.
SEMANTIC_STATUS_PREVIOUS_TOKENS = {
    "schema_ok_semantic_probe_not_linked": "pass",
    "not_applicable_no_cea_table": "pass_or_fail_by_backend",
    "fail": "fail",
}


def semantic_status_for_backend(backend: Any) -> Tuple[str, str, str]:
    """Cleanup P3c: (token, basis, previous_token) for summary validation.

    * backends that do not use CEA (analytic reference_constant, the
      synthetic self-test table) -> 'not_applicable_no_cea_table'; the
      inherited PressureTableBackend.cache_schema of the synthetic backend is
      deliberately ignored because it carries no CEA semantics.
    * CEA table backends with cache schema >= 4 ->
      'schema_ok_semantic_probe_not_linked' (schema number only).
    * CEA table backends with schema < 4 -> 'fail'.
    * CEA table backends with an unknown schema (not a finite real number,
      e.g. None, a string or a bool) -> 'fail', and the basis says so.
    previous_token reproduces the pre-cleanup rule
    'pass' if getattr(backend, 'cache_schema', 0) >= 4 else 'fail'
    for numeric schemas; for an unknown schema (where that expression could
    not be evaluated) it is 'fail'.  Never raises.
    """
    raw_schema = getattr(backend, "cache_schema", 0)
    schema = _finite_number_or_none(raw_schema)
    previous = "pass" if (schema is not None and schema >= 4) else "fail"
    name = str(getattr(backend, "name", type(backend).__name__))
    if not getattr(backend, "uses_cea", False):
        return ("not_applicable_no_cea_table",
                f"backend '{name}' does not use CEA; the CEA mass-basis "
                "semantics check does not apply",
                previous)
    if schema is None:
        return ("fail",
                f"backend '{name}': unknown cache schema "
                f"{_report_value_text(raw_schema, 80)}; treated as fail",
                previous)
    if schema >= 4:
        return ("schema_ok_semantic_probe_not_linked",
                f"backend '{name}': cache schema number {raw_schema} only; no CEA "
                "semantic probe result is linked to this run",
                previous)
    return ("fail",
            f"backend '{name}': CEA table schema {raw_schema} < 4 (obsolete mass "
            "basis semantics)",
            previous)


def out_of_validated_range_error(
    two_phase_mode: str,
    gas_mass_fraction_range: Sequence[float],
    pressure_range_Pa: Sequence[float],
) -> dict:
    """Build the stable error record for an out-of-scope two-phase table."""
    yg_min = float(gas_mass_fraction_range[0])
    return {
        "code": "E_OUT_OF_VALIDATED_RANGE",
        "severity": "error",
        "message": (
            "The two-phase thermochemistry table reaches a gas mass fraction "
            f"Yg_min={yg_min:.6g}, at or below the out-of-scope threshold "
            f"{MODEL_STATUS_OUT_OF_SCOPE_YG_THRESHOLD:.6g}. This range is not "
            "validated for the current model; do not run the ballistics "
            "simulation as a validated case."
        ),
        "context": {
            "model_status": "out_of_scope",
            "two_phase_mode": two_phase_mode,
            "thermochemistry_table_Yg_range": [
                float(value) for value in gas_mass_fraction_range
            ],
            "thermochemistry_pressure_range_Pa": [
                float(value) for value in pressure_range_Pa
            ],
            "validated_in_scope_requires_Yg_min_greater_than": (
                MODEL_STATUS_VALIDATED_YG_THRESHOLD
            ),
            "out_of_scope_if_Yg_min_at_or_below": (
                MODEL_STATUS_OUT_OF_SCOPE_YG_THRESHOLD
            ),
        },
    }


PHASE_FRACTION_VALIDATION_ERROR_CODE = "E_PHASE_FRACTION_VALIDATION_INVALID"
_PHASE_FRACTION_KEYS = ("gas_mass_fraction_min", "gas_mass_fraction_max")


class PhaseFractionValidationError(ThermochemistryError):
    """An active two-phase backend did not provide a usable Yg range.

    Distinct from ``E_OUT_OF_VALIDATED_RANGE``: that code means a valid Yg
    range was classified as out of scope; this one means the range itself is
    missing, non-numeric, non-finite or outside 0 < Yg_min <= Yg_max <= 1, so
    no applicability classification is possible.
    """

    code = PHASE_FRACTION_VALIDATION_ERROR_CODE

    def __init__(self, record: Dict[str, Any]):
        super().__init__(record["message"])
        self.record = dict(record)
        self.context = dict(record.get("context", {}))


def _two_phase_is_active(config: Config, backend: Any) -> bool:
    return (
        config.two_phase_mode == TWO_PHASE_HOMOGENEOUS
        and bool(getattr(backend, "supports_two_phase", False))
    )


def phase_fraction_validation_error(
    config: Config,
    backend: Any,
) -> Optional[dict]:
    """Return a structured error if an ACTIVE two-phase backend lacks a usable Yg range.

    Returns None when two-phase is not active (legacy mode or a backend without
    two-phase support): the Yg range is not used there.  No default is ever
    substituted for a missing value.
    """
    if not _two_phase_is_active(config, backend):
        return None
    table = getattr(backend, "property_validation", None)
    problems: List[str] = []
    values: Dict[str, float] = {}
    reported: Dict[str, str] = {}
    if not isinstance(table, dict):
        problems.append("property_validation_missing")
        table = {}
    for key in _PHASE_FRACTION_KEYS:
        if key not in table:
            problems.append(f"missing:{key}")
            reported[key] = "<missing>"
            continue
        raw = table[key]
        reported[key] = _report_value_text(raw)
        if isinstance(raw, bool) or not isinstance(raw, (int, float, np.integer, np.floating)):
            problems.append(f"not_a_number:{key}")
            continue
        try:
            value = float(raw)
        except Exception:  # noqa: BLE001 - KeyboardInterrupt / SystemExit are not caught
            # e.g. a Python int beyond the float range (10**400), or a float
            # subclass whose __float__ raises
            problems.append(f"not_representable:{key}")
            continue
        if not math.isfinite(value):
            problems.append(f"non_finite:{key}")
            continue
        values[key] = value
    if len(values) == len(_PHASE_FRACTION_KEYS):
        lo, hi = values["gas_mass_fraction_min"], values["gas_mass_fraction_max"]
        if not (0.0 < lo <= hi <= 1.0):
            problems.append("invalid_range")
    if not problems:
        return None
    return {
        "code": PHASE_FRACTION_VALIDATION_ERROR_CODE,
        "severity": "error",
        "message": (
            "The active two-phase thermochemistry backend did not provide a "
            "usable gas mass fraction range (" + ", ".join(problems) + "). "
            "Model applicability cannot be classified; this is not an "
            "out-of-scope result and no default range is assumed."
        ),
        "context": {
            "problems": problems,
            "reported_values": reported,
            "required": "finite numbers with 0 < Yg_min <= Yg_max <= 1",
            "two_phase_mode": config.two_phase_mode,
            "backend": str(getattr(backend, "name", type(backend).__name__)),
            "model_status": None,
        },
    }


def model_applicability_error(
    config: Config,
    backend: ThermochemistryBackend,
) -> Optional[dict]:
    """Return a fail-fast error if an active two-phase table is out of scope.

    A missing / non-finite / invalid Yg range is reported first, with its own
    code (E_PHASE_FRACTION_VALIDATION_INVALID), never as out-of-scope.
    """
    if not _two_phase_is_active(config, backend):
        return None
    invalid = phase_fraction_validation_error(config, backend)
    if invalid is not None:
        return invalid

    table = backend.property_validation
    yg_range = [
        float(table["gas_mass_fraction_min"]),
        float(table["gas_mass_fraction_max"]),
    ]
    two_phase_active = True
    status = model_status_from_yg(
        two_phase_active,
        config.two_phase_mode,
        yg_range[0],
    )
    if status != "out_of_scope":
        return None

    return out_of_validated_range_error(
        config.two_phase_mode,
        yg_range,
        [float(backend.pmin), float(backend.pmax)],
    )


def append_model_applicability_warnings(summary: dict) -> None:
    """Append machine-readable two-phase warnings and scope errors."""
    warnings = summary.setdefault("warnings", [])
    existing_codes = {
        warning.get("code")
        for warning in warnings
        if isinstance(warning, dict)
    }
    two_phase = summary.get("two_phase_model", {})
    if not isinstance(two_phase, dict) or not two_phase.get("active"):
        return

    validation = summary.get("validation", {})
    model_status = validation.get("model_status")
    if (
        model_status == "out_of_scope"
        and "E_OUT_OF_VALIDATED_RANGE" not in existing_codes
    ):
        yg_range = validation.get(
            "validated_Yg_range",
            [
                summary.get("gas_mass_fraction_table_min", 1.0),
                summary.get("gas_mass_fraction_table_max", 1.0),
            ],
        )
        pressure_range = validation.get(
            "validated_pressure_range_Pa", [0.0, 0.0]
        )
        warnings.append(
            out_of_validated_range_error(
                two_phase.get("mode", TWO_PHASE_HOMOGENEOUS),
                yg_range,
                pressure_range,
            )
        )
        existing_codes.add("E_OUT_OF_VALIDATED_RANGE")

    if (
        model_status == "approximate_with_limitations"
        and "W_HIGH_CONDENSATE_MODEL_LIMITATION" not in existing_codes
    ):
        scope_description = "in the approximate-with-limitations Yg scope"
        warnings.append(
            {
                "code": "W_HIGH_CONDENSATE_MODEL_LIMITATION",
                "severity": "warning",
                "message": (
                    "The thermochemistry table is "
                    f"{scope_description} (model_status={model_status!r}). "
                    "The current homogeneous-equilibrium model provides "
                    "equilibrium phase bookkeeping only; it does not model "
                    "particle slip, finite-rate condensation, finite "
                    "condensed-phase volume, or phase-density lag. Treat "
                    "performance in this range as approximate, not as "
                    "validated high-condensate performance."
                ),
                "context": {
                    "model_status": model_status,
                    "gas_mass_fraction_table_min": summary.get(
                        "gas_mass_fraction_table_min"
                    ),
                    "condensed_mass_fraction_table_max": summary.get(
                        "condensed_mass_fraction_table_max"
                    ),
                    "validated_in_scope_requires_Yg_min_greater_than": (
                        MODEL_STATUS_VALIDATED_YG_THRESHOLD
                    ),
                    "out_of_scope_if_Yg_min_at_or_below": (
                        MODEL_STATUS_OUT_OF_SCOPE_YG_THRESHOLD
                    ),
                    "unmodelled_effects": [
                        "particle_slip",
                        "finite_rate_condensation",
                        "finite_condensed_phase_volume",
                        "phase_density_lag",
                    ],
                },
            }
        )
        existing_codes.add("W_HIGH_CONDENSATE_MODEL_LIMITATION")

    # Validation-evidence warning for the high-condensate band.  It is emitted
    # in addition to W_HIGH_CONDENSATE_MODEL_LIMITATION (which describes the
    # unmodelled physics) and does not change model_status, the run verdict or
    # the exit code.  Out-of-scope tables (Yg_min <= 0.65) report
    # E_OUT_OF_VALIDATED_RANGE instead and do not get this warning.
    g4b_yg_min = _finite_number_or_none(summary.get("gas_mass_fraction_table_min"))
    if (
        model_status == "approximate_with_limitations"
        and g4b_yg_min is not None
        and g4b_yg_min <= HIGH_CONDENSATE_G4B_UNVALIDATED_YG_THRESHOLD
        and "W_HIGH_CONDENSATE_G4B_UNVALIDATED" not in existing_codes
    ):
        warnings.append(
            {
                "code": "W_HIGH_CONDENSATE_G4B_UNVALIDATED",
                "severity": "warning",
                "message": (
                    "Validation evidence: the thermochemistry table reaches "
                    f"Yg_min={g4b_yg_min:.6g}, inside the high-condensate band "
                    f"(Yg_min <= {HIGH_CONDENSATE_G4B_UNVALIDATED_YG_THRESHOLD:.6g}). "
                    "Reproduction of reported or experimental data (gate G4b) "
                    "has not been done for this band, so results here are not "
                    "validated predictions. This warning concerns missing "
                    "validation evidence; the unmodelled physics is described by "
                    "W_HIGH_CONDENSATE_MODEL_LIMITATION. It does not change the "
                    "run verdict or the exit code."
                ),
                "context": {
                    "model_status": model_status,
                    "gas_mass_fraction_table_min": g4b_yg_min,
                    "high_condensate_band_if_Yg_min_at_or_below": (
                        HIGH_CONDENSATE_G4B_UNVALIDATED_YG_THRESHOLD
                    ),
                    "out_of_scope_if_Yg_min_at_or_below": (
                        MODEL_STATUS_OUT_OF_SCOPE_YG_THRESHOLD
                    ),
                    "g4b_reported_or_experimental_reproduction": "not_done",
                    "related_warning": "W_HIGH_CONDENSATE_MODEL_LIMITATION",
                },
            }
        )
        existing_codes.add("W_HIGH_CONDENSATE_G4B_UNVALIDATED")

    condensed_max = float(
        summary.get("maximum_condensed_mass_fraction_during_run", 0.0)
    )
    if (
        two_phase.get("nozzle_entrainment") == "complete"
        and condensed_max > _PHASE_FRACTION_ROUNDING_TOL
        and "W_TWO_PHASE_NOZZLE_APPROXIMATION" not in existing_codes
    ):
        warnings.append(
            {
                "code": "W_TWO_PHASE_NOZZLE_APPROXIMATION",
                "severity": "warning",
                "message": (
                    "The two-phase nozzle assumes complete entrainment and "
                    "equal gas/condensed-phase exit velocity. Particle slip "
                    "is not modelled, so momentum thrust may be overestimated "
                    "when condensed particles have significant relative "
                    "velocity."
                ),
                "context": {
                    "nozzle_entrainment": two_phase.get("nozzle_entrainment"),
                    "particle_slip_modelled": False,
                    "maximum_condensed_mass_fraction_during_run": condensed_max,
                    "roundoff_threshold": _PHASE_FRACTION_ROUNDING_TOL,
                },
            }
        )


def append_thermal_efficiency_warning(
    summary: dict,
    config: Config,
    backend: ThermochemistryBackend,
) -> None:
    """Report that eta_T0 scales CEA temperature without re-solving HP state."""
    if getattr(backend, "name", None) not in (
        BACKEND_CEA_PYTHON,
        BACKEND_CEA_LEGACY_EXECUTABLE,
    ):
        return
    eta_t0 = float(config.eta_T0)
    if not math.isfinite(eta_t0) or eta_t0 == 1.0:
        return

    warnings = summary.setdefault("warnings", [])
    code = "W_THERMAL_EFFICIENCY_PARTIAL"
    if any(
        isinstance(warning, dict) and warning.get("code") == code
        for warning in warnings
    ):
        return

    direction = "above" if eta_t0 > 1.0 else "below"
    warnings.append(
        {
            "code": code,
            "severity": "warning",
            "message": (
                f"temperature_efficiency eta_T0={eta_t0:.6g} scales the HP CEA "
                f"temperature {direction} its equilibrium value without "
                "re-solving the thermochemical state. gamma_s, M_CEA, species "
                "fractions and phase fractions remain from the unscaled CEA "
                "solution; direct comparisons of density, sonic velocity or "
                "CEA theoretical performance are not like-for-like. Use "
                "eta_T0=1 for direct CEA validation comparisons."
            ),
            "context": {
                "temperature_efficiency_eta_T0": eta_t0,
                "temperature_relation": "T_used = eta_T0 * T_CEA",
                "cea_state_re_solved_after_scaling": False,
                "unscaled_equilibrium_outputs": [
                    "gamma_s",
                    "M_CEA",
                    "species_fractions",
                    "phase_fractions",
                ],
                "direct_cea_comparisons_not_like_for_like": [
                    "density",
                    "sonic_velocity",
                    "theoretical_rocket_performance",
                ],
                "recommended_eta_T0_for_direct_cea_validation": 1.0,
            },
        }
    )


def append_external_cea_database_warning(
    summary: dict,
    backend: ThermochemistryBackend,
) -> None:
    """Report thermo.lib / trans.lib files that GRIBS ignored (P5-1).

    The cea_python backend always loads the database of the official package;
    files in the working directory or in CEA_DATA_DIR, which a plain
    ``import cea`` would have loaded instead, are listed here.
    """
    identity = getattr(backend, "cea_database", None)
    if not isinstance(identity, dict):
        return
    ignored = identity.get("ignored_external_files") or []
    if not ignored:
        return
    warnings = summary.setdefault("warnings", [])
    code = "W_EXTERNAL_CEA_DATABASE_IGNORED"
    if any(isinstance(w, dict) and w.get("code") == code for w in warnings):
        return
    warnings.append(
        {
            "code": code,
            "severity": "warning",
            "message": (
                "CEA database files were found in the working directory or in "
                "CEA_DATA_DIR and were ignored; the thermo.lib of the official "
                "CEA package was used. Results do not depend on these files."
            ),
            "context": {
                "ignored_files": [dict(f) for f in ignored],
                "used_thermo_lib": identity.get("thermo_lib_path"),
                "used_thermo_lib_sha256": identity.get("thermo_lib_sha256"),
            },
        }
    )


def summarize(res: dict, h: dict, c: Config,
              backend: Optional[ThermochemistryBackend] = None) -> dict:
    t = h["t"]
    burn = h["phase"] == "burn"
    m_prop = propellant_mass(c)
    At0, Ae = res["At0"], res["Ae"]
    t_bo = float(t[burn][-1])

    i_pmax = int(np.argmax(h["p0"]))
    i_fmax = int(np.argmax(h["F"]))
    I_burn = float(h["impulse"][burn][-1])
    I_tot = float(h["impulse"][-1])
    m_out_tot = float(h["m_out"][-1])
    m_gen_tot = float(h["m_gen"][burn][-1])

    m_out_burn = float(h["m_out"][burn][-1])
    pAt_int = float(_TRAPZ(h["p0"][burn] * h["At"][burn], t[burn]))
    # cstar_eff_m_s (and legacy_mixed_interval_cstar_m_s) keep the mixed definition: the burning-phase integral of
    # p0*At divided by the mass discharged over the whole run (burn + blowdown).
    cstar_eff = pAt_int / m_out_tot if m_out_tot > 0 else float("nan")
    # P5-6 (F-11): cstar_burn_interval_m_s is burn over burn: the burning-phase integral of p0*At divided by the mass
    # discharged up to burnout.  Equal to cstar_eff_m_s when there is no blowdown.
    cstar_burn = pAt_int / m_out_burn if m_out_burn > 0 else float("nan")
    CF_eff = I_burn / pAt_int if pAt_int > 0 else float("nan")
    p_mean = float(pAt_int / (At0 * (t_bo - t[0]))) if t_bo > t[0] else float("nan")

    # --- c* indicators (A-5.1-5 c* reorganization) ---
    # local_model_cstar: median of p0*At/mdot_total at each time step during
    # the burn phase (a per-step diagnostic of the nozzle model).
    cstar_local = np.asarray(h.get("cstar", []), dtype=float)
    cstar_local_burn = cstar_local[burn] if cstar_local.size > 0 else np.array([])
    cstar_local_median = (float(np.median(cstar_local_burn))
                          if cstar_local_burn.size > 0 else float("nan"))
    # full-interval c*  (burn + blowdown)
    pAt_full = float(_TRAPZ(h["p0"] * h["At"], t))
    cstar_full = pAt_full / m_out_tot if m_out_tot > 0 else float("nan")

    # One authoritative post-run consistency report.  The legacy metric is
    # retained under its existing field name for output compatibility.
    post_run = post_run_consistency_report(res, h)
    mass_err = post_run["eos_history"][
        "legacy_instantaneous_mass_relative_error"
    ]

    # two-phase statistics over the computed history
    Yg_h = np.asarray(h["Y_gas"], dtype=float)
    Yc_h = np.asarray(h["Y_condensed"], dtype=float)

    p_star = gribs.nozzle.choke_limit_pressure(c, At0, Ae)
    transition = nozzle_transition_diagnostics(h, c, Ae)

    mg0 = c.rho_p * base_burn_rate(c.p0_init, c) * geometry(0.0, c)[2]
    mo0 = gribs.nozzle.nozzle_state(c.p0_init, At0, Ae, c)["mdot_total"]

    if backend is None:
        backend = gribs.thermo.thermochemistry()
    pv_table = getattr(backend, "property_validation", {})
    ledger = mass_energy_ledger(res, h, c)
    two_phase_active = (c.two_phase_mode == TWO_PHASE_HOMOGENEOUS
                        and bool(getattr(backend, "supports_two_phase", False)))
    # An active two-phase run needs a usable Yg range; never classify a
    # missing / non-finite / invalid range with a default (direct callers of
    # summarize() bypass the pre-integration gate in main()).
    _pf_error = phase_fraction_validation_error(c, backend)
    if _pf_error is not None:
        raise PhaseFractionValidationError(_pf_error)
    # A-5.1-6: model_status determined by Yg range in the property table
    _yg_table_min = (float(pv_table["gas_mass_fraction_min"]) if two_phase_active
                     else float(pv_table.get("gas_mass_fraction_min", 1.0)))
    _model_status = model_status_from_yg(
        two_phase_active,
        c.two_phase_mode,
        _yg_table_min,
    )
    _semantic_status, _semantic_basis, _semantic_previous = (
        semantic_status_for_backend(backend))
    _has_table = isinstance(getattr(backend, "table", None), dict) and \
        "pressure_Pa" in getattr(backend, "table", {})
    s = dict(
        web_mm=web_thickness(c) * 1e3,
        propellant_density_kg_m3=c.rho_p,
        propellant_mass_g=m_prop * 1e3,
        burn_time_s=t_bo,
        total_time_s=float(t[-1]),
        burn_stop_reason=res["stop_burn"],
        blowdown_stop_reason=res["stop_blow"],
        initial_regime=res["init_regime"],
        initial_mdot_gen_kg_s=mg0,
        initial_mdot_out_kg_s=mo0,
        choke_limit_pressure_Pa=p_star,
        p_max_Pa=float(h["p0"][i_pmax]), t_p_max_s=float(t[i_pmax]),
        p_min_burn_Pa=float(np.min(h["p0"][burn])),
        t_p_min_s=float(t[burn][int(np.argmin(h["p0"][burn]))]),
        p_mean_burn_Pa=p_mean,
        F_max_N=float(h["F"][i_fmax]), t_F_max_s=float(t[i_fmax]),
        F_mean_burn_N=float(I_burn / (t_bo - t[0])) if t_bo > t[0] else float("nan"),
        impulse_burn_Ns=I_burn,
        impulse_total_Ns=I_tot,
        impulse_blowdown_Ns=I_tot - I_burn,
        # Isp_s: total impulse per initial propellant mass (unchanged).  P5-6 (F-12): the same impulse per mass
        # discharged over the whole run is Isp_discharged_mass_s.
        Isp_s=I_tot / (m_prop * G0) if m_prop > 0 else float("nan"),
        Isp_discharged_mass_s=I_tot / (m_out_tot * G0) if m_out_tot > 0 else float("nan"),
        cstar_eff_m_s=cstar_eff,
        local_model_cstar_m_s=cstar_local_median,
        cstar_burn_interval_m_s=cstar_burn,
        cstar_full_interval_m_s=cstar_full,
        legacy_mixed_interval_cstar_m_s=cstar_eff,
        CF_eff=CF_eff,
        Kn_initial=float(h["Kn"][0]), Kn_max=float(np.max(h["Kn"][burn])),
        Kn_final=float(h["Kn"][burn][-1]),
        throat_radius_final_mm=float(h["Rt"][-1]) * 1e3,
        propellant_mass_balance_error=abs(m_gen_tot / m_prop - 1.0) if m_prop > 0 else float("nan"),
        total_mass_consistency_error=mass_err,
        total_mass_consistency_error_invariant=(
            post_run["eos_history"]["invariant_relative_error"]
        ),
        total_mass_consistency_reference_mass_kg=(
            post_run["eos_history"]["reference_mass_kg"]
        ),
        total_mass_consistency_metric_version=(
            post_run["eos_history"]["metric_version"]
        ),
        post_run_consistency=post_run,
        # ---------------- Phase A-4 ledger ----------------
        ledger=ledger,
        # ---------------- two-phase model ----------------
        two_phase_model={
            "mode": c.two_phase_mode,
            "condensed_volume": c.two_phase_condensed_volume,
            "nozzle_entrainment": c.two_phase_nozzle_entrainment,
            "active": bool(two_phase_active),
        },
        gas_mass_fraction_min=float(np.min(Yg_h)),
        gas_mass_fraction_max=float(np.max(Yg_h)),
        condensed_mass_fraction_min=float(np.min(Yc_h)),
        condensed_mass_fraction_max=float(np.max(Yc_h)),
        maximum_condensed_mass_fraction_during_run=float(np.max(Yc_h)),
        mean_condensed_mass_fraction_during_burn=float(np.mean(Yc_h[burn]))
        if np.any(burn) else float("nan"),
        condensed_mass_fraction_at_max_pressure=float(Yc_h[int(np.argmax(h["p0"]))]),
        gas_mass_fraction_table_min=float(pv_table.get("gas_mass_fraction_min", float("nan"))),
        gas_mass_fraction_table_max=float(pv_table.get("gas_mass_fraction_max", float("nan"))),
        condensed_mass_fraction_table_min=float(pv_table.get("condensed_mass_fraction_min", float("nan"))),
        condensed_mass_fraction_table_max=float(pv_table.get("condensed_mass_fraction_max", float("nan"))),
        min_pressure_factor_1_minus_pPsi_over_Psi=float(
            pv_table.get("min_pressure_factor_1_minus_pPsi_over_Psi", float("nan"))),
        nozzle_entrainment_assumption=(
            "complete entrainment: the isentropic formula evaluated with R_eff "
            "returns mdot_total directly (no additional /Yg); both phases "
            "leave with the same exit velocity (particle slip NOT modelled; "
            "momentum thrust may be overestimated for large condensed "
            "fractions)" if two_phase_active else
            "single-phase legacy: total flow = gas reference flow (Yg := 1)"),
        condensed_volume_assumption=(
            "condensed-phase volume neglected against the chamber free volume; "
            "the gas phase occupies Vg" if two_phase_active else
            "not applicable (single-phase legacy mode)"),
        property_extrapolation=bool(np.any(h["extrap"])),
        property_extrapolation_points=int(np.count_nonzero(h["extrap"])),
        property_extrapolation_fraction=float(np.mean(h["extrap"])),
        property_policy=c.property_policy,
        thermochemistry_backend=backend.name,
        property_range_min_Pa=float(backend.pmin),
        property_range_max_Pa=float(backend.pmax),
        pressure_min_Pa=float(np.min(h["p0"])),
        pressure_max_Pa=float(np.max(h["p0"])),
        regimes_visited=sorted(set(h["regime"].tolist())),
        nozzle_transition=transition,
        # ---------------- validation status (A-5.1-5 / A-5.1-6) ----------------
        validation={
            "semantic_status": _semantic_status,
            "semantic_status_basis": _semantic_basis,
            "model_status_basis": model_status_basis(
                two_phase_active, c.two_phase_mode),
            "status_previous_tokens": {
                "semantic_status": _semantic_previous,
                "model_status": MODEL_STATUS_PREVIOUS_TOKENS.get(
                    _model_status, _model_status),
            },
            "numerical_status": (
                "pass" if (post_run["passed"] and ledger["passed"]) else "fail"
            ),
            "model_status": _model_status,
            **_validation_range_fields(backend, pv_table, _has_table),
            "mass_basis_semantics_version": getattr(
                backend, "mass_basis_semantics_version", None),
            "notes": (
                "Two-phase model provides equilibrium phase bookkeeping only. "
                "No particle slip, finite-rate condensation, condensed volume, "
                "or phase density lag.  See two_phase_capabilities in "
                "summary.thermochemistry."
                if two_phase_active else
                "Single-phase legacy mode (v0.4.3-alpha regression anchor)."
            ),
        },
        warnings=(
            [
                {
                    "code": "W_NOZZLE_TRANSITION",
                    "severity": "warning",
                    "message": (
                        "The direct-switch nozzle model traversed the "
                        "C-D nozzle transition band. Thrust and impulse in "
                        "this band are model-dependent because an internal-"
                        "shock solution is not implemented."
                    ),
                    "context": {
                        "policy": c.transition_policy,
                        "duration_s": transition["duration_s"],
                        "duration_fraction": transition["duration_fraction"],
                        "impulse_Ns": transition["impulse_Ns"],
                        "impulse_fraction": transition["impulse_fraction"],
                        "interval_count": transition["interval_count"],
                    },
                }
            ]
            if transition["warning_code"] == "W_NOZZLE_TRANSITION"
            else []
        ) + (
            # S-8: Warn when blowdown was truncated by time limit.
            # The total impulse and Isp may be underestimated.
            [{
                "code": "W_INTEGRATION_TRUNCATED",
                "severity": "warning",
                "message": (
                    "The blowdown phase was terminated by the time limit "
                    f"({c.blowdown_tmax} s), not by reaching ambient pressure.  "
                    "Total impulse and Isp may be underestimated."
                ),
                "context": {
                    "phase": "blowdown",
                    "time_limit_s": c.blowdown_tmax,
                    "total_time_s": float(t[-1]),
                    "impulse_blowdown_Ns": float(I_tot - I_burn),
                },
            }]
            if res["stop_blow"] == "time-limit" and c.blowdown
            else []
        ) + (
            # Also warn if the burning phase was truncated by time limit.
            [{
                "code": "W_INTEGRATION_TRUNCATED",
                "severity": "warning",
                "message": (
                    "The burning phase was terminated by the time limit "
                    f"({c.t_max} s), not by burnout.  Results are incomplete."
                ),
                "context": {
                    "phase": "burn",
                    "time_limit_s": c.t_max,
                    "total_time_s": float(t[-1]),
                },
            }]
            if res["stop_burn"] == "time-limit"
            else []
        ),
        event_times=dict(res["t_events"]),
    )
    append_throat_exceeds_exit_area_warning(s, res, h, c)
    append_model_applicability_warnings(s)
    append_thermal_efficiency_warning(s, c, backend)
    append_external_cea_database_warning(s, backend)
    append_property_table_crossing_warning(s, res)
    append_run_state_warnings(s, h)
    return s


def append_property_table_crossing_warning(summary: dict, res: dict) -> None:
    """Warn when the run crossed p_fit_min / p_fit_max (property policy "clamp").

    Outside the table the chamber properties are held at the nearest limit and
    their pressure derivative is zero, so the results there rest on a clamped
    property model.  Diagnostic only: histories and the summary schema are
    unchanged apart from one entry in ``warnings``.
    """
    budget = (res.get("integration_budget") if isinstance(res, dict) else None) or {}
    burn = list(budget.get("table_restarts_burn") or [])
    blow = list(budget.get("table_restarts_blowdown") or [])
    if not (burn or blow):
        return
    warnings = summary.setdefault("warnings", [])
    if any(isinstance(w, dict) and w.get("code") == "W_PROPERTY_TABLE_RANGE_CROSSED"
           for w in warnings):
        return
    warnings.append({
        "code": "W_PROPERTY_TABLE_RANGE_CROSSED",
        "severity": "warning",
        "message": (
            "The chamber pressure crossed the property-table limit "
            "(p_fit_min / p_fit_max) during the run.  With property_policy "
            "'clamp' the properties outside the table are held at the nearest "
            "limit and their pressure dependence is suppressed; results in "
            "that pressure range are not validated by the table."),
        "context": {
            "restarts_burn": burn,
            "restarts_blowdown": blow,
            "restart_count": len(burn) + len(blow),
        },
    })


W_THROAT_EXCEEDS_EXIT_AREA = "W_THROAT_EXCEEDS_EXIT_AREA"
THROAT_EXCEEDS_EXIT_EVENT = "throat_exceeds_exit_s"


def throat_exceeds_exit_area_time(res: dict, h: dict, c: Config) -> Optional[float]:
    """P5-3 (F-06): first time at which the eroding throat area exceeds the fixed exit area, or None.

    Found after the integration, which is not changed: the first sample of the history with At > Ae
    (gribs.nozzle.throat_exceeds_exit_area) and the sample before it bracket the crossing, and the crossing is
    located on the dense solution of that phase.  If the first sample is already beyond Ae, its time is returned.
    The throat radius does not decrease (the erosion rate is >= 0), so the first crossing is the only one.
    """
    t = np.asarray(h["t"], dtype=float)
    over = np.flatnonzero(gribs.nozzle.throat_exceeds_exit_area(
        np.asarray(h["Rt"], dtype=float), c.R_t0, c.eps_nozzle))
    if over.size == 0:
        return None
    k = int(over[0])
    if k == 0:
        return float(t[0])
    limit = gribs.nozzle.exit_area_throat_radius(c.R_t0, c.eps_nozzle)
    t_burn_end = float(res["sol_burn"].t[-1])

    def excess(s: float) -> float:
        sol = res["sol_burn"] if (s <= t_burn_end or res["sol_blow"] is None) else res["sol_blow"]
        return float(sol.sol(s)[IRT]) - limit

    a, b = float(t[k - 1]), float(t[k])
    if excess(a) >= 0.0:
        return a
    if not excess(b) > 0.0:
        return b
    return float(brentq(excess, a, b, xtol=np.finfo(float).tiny,
                        rtol=4.0 * np.finfo(float).eps, maxiter=200))


def append_throat_exceeds_exit_area_warning(summary: dict, res: dict, h: dict, c: Config) -> None:
    """P5-3 (F-06): record when the throat area exceeds the fixed exit area.

    Adds event_times["throat_exceeds_exit_s"] (empty when the run stays at At <= Ae, as the other events) and,
    when the crossing occurs, the warning W_THROAT_EXCEEDS_EXIT_AREA.  Nothing else is changed: the run continues,
    the nozzle model and all computed values are unchanged, and the warning does not change the run verdict.
    """
    t_cross = throat_exceeds_exit_area_time(res, h, c)
    summary["event_times"][THROAT_EXCEEDS_EXIT_EVENT] = [] if t_cross is None else [t_cross]
    if t_cross is None:
        return
    Rt = np.asarray(h["Rt"], dtype=float)
    k = int(np.flatnonzero(gribs.nozzle.throat_exceeds_exit_area(Rt, c.R_t0, c.eps_nozzle))[0])
    summary.setdefault("warnings", []).append({
        "code": W_THROAT_EXCEEDS_EXIT_AREA,
        "severity": "warning",
        "message": (
            f"The eroding throat area exceeded the fixed exit area (At > Ae) at t = {t_cross:.6f} s. "
            "This is outside the model domain: the exit area is fixed, and from this time the nozzle is "
            "computed with the exit as the smallest section. Thrust and mass flow after this time are not "
            "validated; at an unchoking with At > Ae they change discontinuously by the factor At/Ae."
        ),
        "context": {
            "time_s": t_cross,
            "phase": str(h["phase"][k]),
            "eps_nozzle": c.eps_nozzle,
            "initial_throat_radius_m": c.R_t0,
            "exit_area_throat_radius_m": gribs.nozzle.exit_area_throat_radius(c.R_t0, c.eps_nozzle),
            "final_throat_radius_m": float(Rt[-1]),
            "final_throat_to_exit_area_ratio": float((Rt[-1] / c.R_t0) ** 2 / c.eps_nozzle),
        },
    })


W_NON_NOMINAL_TERMINATION = "W_NON_NOMINAL_TERMINATION"
W_SUBSONIC_BURNING = "W_SUBSONIC_BURNING"
W_PROPERTY_TABLE_EXTRAPOLATED = "W_PROPERTY_TABLE_EXTRAPOLATED"
RUN_STATE_KEYS = (
    "solver_completed",
    "mass_consistency_passed",
    "nominal_termination_reached",
    "model_domain_passed",
    "warnings_present",
    "results_usable",
)


def nominal_termination(burn_stop_reason: Any, blowdown_stop_reason: Any) -> bool:
    """P5-5: the burning phase ended by burnout and the blowdown, when it was run, reached ambient pressure.

    After a burnout the blowdown is "not-run" only when it is disabled, which is a nominal end.
    """
    return burn_stop_reason == "burnout" and blowdown_stop_reason in ("ambient", "not-run")


def append_run_state_warnings(summary: dict, h: dict) -> None:
    """P5-5 (F-16): warnings for run states that the run-state keys report.

    W_NON_NOMINAL_TERMINATION: the run did not end nominally (see nominal_termination).
    W_SUBSONIC_BURNING: the nozzle was subsonic at one or more samples of the burning phase.
    W_PROPERTY_TABLE_EXTRAPOLATED: the chamber pressure left the pressure range of the chamber-property table.
    The last two do not make the results unusable.  None of them changes the run verdict, the exit code or a
    computed value.  A code already present is not added again.
    """
    warnings = summary.setdefault("warnings", [])
    existing = {w.get("code") for w in warnings if isinstance(w, dict)}
    burn_reason = summary.get("burn_stop_reason")
    blow_reason = summary.get("blowdown_stop_reason")
    if not nominal_termination(burn_reason, blow_reason) and W_NON_NOMINAL_TERMINATION not in existing:
        warnings.append({
            "code": W_NON_NOMINAL_TERMINATION,
            "severity": "warning",
            "message": (
                f"The run did not end nominally: the burning phase ended by {burn_reason!r} and the blowdown "
                f"by {blow_reason!r} (nominal: burnout, then ambient pressure or no blowdown). The results "
                "describe a truncated or interrupted run and are reported as not usable "
                "(results_usable = false)."
            ),
            "context": {
                "burn_stop_reason": burn_reason,
                "blowdown_stop_reason": blow_reason,
                "burn_time_s": summary.get("burn_time_s"),
                "total_time_s": summary.get("total_time_s"),
            },
        })
    phase = np.asarray(h["phase"])
    regime = np.asarray(h["regime"])
    burn = phase == "burn"
    subsonic = np.flatnonzero(burn & (regime == "subsonic"))
    if subsonic.size and W_SUBSONIC_BURNING not in existing:
        warnings.append({
            "code": W_SUBSONIC_BURNING,
            "severity": "warning",
            "message": (
                "The nozzle was not choked (subsonic flow) during part of the burning phase. This is a valid "
                "branch of the nozzle model, but burn-phase c* indicators lose their meaning there. It does not "
                "make the results unusable."
            ),
            "context": {
                "subsonic_burn_samples": int(subsonic.size),
                "burn_samples": int(np.count_nonzero(burn)),
                "first_subsonic_burn_time_s": float(np.asarray(h["t"], dtype=float)[subsonic[0]]),
                "initial_regime": summary.get("initial_regime"),
            },
        })
    if summary.get("property_extrapolation") and W_PROPERTY_TABLE_EXTRAPOLATED not in existing:
        policy = summary.get("property_policy")
        fraction = _finite_number_or_none(summary.get("property_extrapolation_fraction"))
        warnings.append({
            "code": W_PROPERTY_TABLE_EXTRAPOLATED,
            "severity": "warning",
            "message": (
                "The chamber pressure left the pressure range of the chamber-property table at "
                f"{summary.get('property_extrapolation_points')} samples"
                + ("" if fraction is None else f" ({100.0 * fraction:.3g} % of the samples)")
                + "; outside it the properties were "
                + ("held at the nearest table limit" if policy == "clamp" else "extrapolated")
                + f" (outside_range_policy = {policy!r}). This is reported as a warning; it is not a "
                "model-domain violation and does not make the results unusable. To cover the run, widen the "
                "table range (thermochemistry.pressure_range.minimum_Pa / maximum_Pa)."
            ),
            "context": {
                "policy": policy,
                "points": summary.get("property_extrapolation_points"),
                "fraction": summary.get("property_extrapolation_fraction"),
                "property_range_min_Pa": summary.get("property_range_min_Pa"),
                "property_range_max_Pa": summary.get("property_range_max_Pa"),
                "pressure_min_Pa": summary.get("pressure_min_Pa"),
                "pressure_max_Pa": summary.get("pressure_max_Pa"),
                # added by the integration of A-9 (the only two A-9 keys without an equivalent above)
                "points_total": int(np.asarray(h["extrap"]).size),
                "basis": "sampled_history_all_phases",
            },
        })


def add_run_state_keys(summary: dict) -> None:
    """P5-5 (F-16): set the six run-state keys at the top level of the summary of a completed run.

    Called after all warnings are appended (after append_run_verdict_warnings).  A run that raises has no summary,
    and no keys are made for it.  The keys do not change run_passed or the exit code.
      solver_completed: both phases were integrated and the chamber pressure is finite (a run that reaches this
        point has completed the integration; the pressure extremes are NaN or infinite otherwise).
      mass_consistency_passed: the post-run mass consistency and the mass ledger passed.
      nominal_termination_reached: see nominal_termination.
      model_domain_passed: no At > Ae (event_times.throat_exceeds_exit_s is empty) and no
        E_OUT_OF_VALIDATED_RANGE.  Leaving the property-table range is not a domain violation (a warning).
      warnings_present: summary["warnings"] is not empty.
      results_usable: the first four together.
    """
    warnings = summary.get("warnings") or []
    codes = {w.get("code") for w in warnings if isinstance(w, dict)}
    solver = (_finite_number_or_none(summary.get("pressure_min_Pa")) is not None
              and _finite_number_or_none(summary.get("pressure_max_Pa")) is not None)
    mass = bool(summary["post_run_consistency"]["passed"] and summary["ledger"]["passed"])
    nominal = nominal_termination(summary.get("burn_stop_reason"), summary.get("blowdown_stop_reason"))
    domain = (not (summary.get("event_times") or {}).get(THROAT_EXCEEDS_EXIT_EVENT)
              and "E_OUT_OF_VALIDATED_RANGE" not in codes)
    summary.update({
        "solver_completed": bool(solver),
        "mass_consistency_passed": mass,
        "nominal_termination_reached": bool(nominal),
        "model_domain_passed": bool(domain),
        "warnings_present": bool(warnings),
        "results_usable": bool(solver and mass and nominal and domain),
    })


def append_run_verdict_warnings(summary: dict) -> None:
    """Append machine-readable warnings for failed run-level verdicts.

    This function changes only the external diagnostic contract.  It does not
    modify simulation histories, ledger equations, or acceptance tolerances.
    """
    warnings = summary.setdefault("warnings", [])
    existing_codes = {
        warning.get("code")
        for warning in warnings
        if isinstance(warning, dict)
    }

    post_run = summary["post_run_consistency"]
    if not post_run["passed"] and "W_POST_RUN_CONSISTENCY" not in existing_codes:
        final = post_run["final_total_mass"]
        eos = post_run["eos_history"]
        warnings.append(
            {
                "code": "W_POST_RUN_CONSISTENCY",
                "severity": "error",
                "message": (
                    "The post-run mass-consistency contract failed. "
                    "Diagnostic output was preserved before returning a "
                    "nonzero process status."
                ),
                "context": {
                    "final_total_mass_passed": final["passed"],
                    "final_total_mass_relative_error": (
                        final["relative_error"]
                    ),
                    "eos_history_passed": eos["passed"],
                    "eos_history_invariant_relative_error": (
                        eos["invariant_relative_error"]
                    ),
                    "eos_history_reference_mass_kg": (
                        eos["reference_mass_kg"]
                    ),
                    "legacy_instantaneous_mass_relative_error": (
                        eos[
                            "legacy_instantaneous_mass_relative_error"
                        ]
                    ),
                    "tolerance": post_run["tolerance"],
                },
            }
        )

    ledger = summary["ledger"]
    if (
        not ledger["passed"]
        and "W_MASS_LEDGER_FAILED" not in existing_codes
    ):
        mass = ledger["mass"]
        total = mass["total_system"]
        propellant = mass["propellant"]
        warnings.append(
            {
                "code": "W_MASS_LEDGER_FAILED",
                "severity": "error",
                "message": (
                    "The mass-ledger acceptance contract failed. "
                    "Diagnostic output was preserved before returning a "
                    "nonzero process status."
                ),
                "context": {
                    "mass_passed": mass["passed"],
                    "total_system_passed": total["passed"],
                    "total_system_signed_residual_kg": (
                        total["signed_residual_kg"]
                    ),
                    "total_system_relative_residual": (
                        total["relative_residual"]
                    ),
                    "total_system_tolerance": total["tolerance"],
                    "propellant_passed": propellant["passed"],
                    "propellant_signed_residual_kg": (
                        propellant["signed_residual_kg"]
                    ),
                    "propellant_relative_residual": (
                        propellant["relative_residual"]
                    ),
                    "propellant_tolerance": propellant["tolerance"],
                    "igniter_mass_commanded_kg": (
                        mass["igniter_mass_commanded_kg"]
                    ),
                    "igniter_mass_injected_kg": (
                        mass["igniter_mass_injected_kg"]
                    ),
                },
            }
        )


def run_verdict_passed(summary: dict) -> bool:
    """Return the combined post-run, ledger, and applicability verdict."""
    scope_error = any(
        isinstance(warning, dict)
        and warning.get("code") == "E_OUT_OF_VALIDATED_RANGE"
        for warning in (summary.get("warnings", []) or [])
    )
    return bool(
        summary["post_run_consistency"]["passed"]
        and summary["ledger"]["passed"]
        and not scope_error
    )


EXIT_CODE_SUCCESS = 0
EXIT_CODE_VERIFICATION_FAILED = 1
EXIT_CODE_MASS_LEDGER_FAILED = 5


def run_verification_exit_code(
    summary: dict,
    *,
    allow_verification_failure: bool = False,
) -> int:
    """Return the CLI status for post-run, ledger, and scope gates.

    Code 5 is reserved for a mass-ledger-only failure. The explicit override
    may downgrade that one gate to 0 for sweeps; it never masks a failed
    post-run consistency or an out-of-validated-range error.
    """
    if run_verdict_passed(summary):
        return EXIT_CODE_SUCCESS

    ledger_failed = not summary["ledger"]["passed"]
    post_run_failed = not summary["post_run_consistency"]["passed"]
    scope_error = any(
        isinstance(warning, dict)
        and warning.get("code") == "E_OUT_OF_VALIDATED_RANGE"
        for warning in (summary.get("warnings", []) or [])
    )
    if ledger_failed and not post_run_failed and not scope_error:
        return (
            EXIT_CODE_SUCCESS
            if allow_verification_failure
            else EXIT_CODE_MASS_LEDGER_FAILED
        )
    return EXIT_CODE_VERIFICATION_FAILED


MASS_LEDGER_TOTAL_RELATIVE_TOLERANCE = 1.0e-6
MASS_LEDGER_PROPELLANT_RELATIVE_TOLERANCE = 1.0e-9
# A-5.1-6 / S-2: ledger propellant tol scales with rtol so that tightening
# the ODE tolerance does not break the gate.  The floor (1e-9) handles the
# numerical residual floor observed in the30-point sweep (§11-Q17).
MASS_LEDGER_PROPELLANT_TOLERANCE_FLOOR = 1.0e-9
MASS_LEDGER_PROPELLANT_TOLERANCE_RTL_FACTOR = 100.0


def mass_energy_ledger(res: dict, h: dict, c: Config) -> dict:
    """Build the Phase A-4 mass ledger and energy-availability metadata."""
    tiny = 1.0e-30

    # A-5.1-6: compute rtol-linked propellant tolerance
    propellant_tol = max(
        MASS_LEDGER_PROPELLANT_TOLERANCE_FLOOR,
        MASS_LEDGER_PROPELLANT_TOLERANCE_RTL_FACTOR * float(c.rtol),
    )

    t = np.asarray(h["t"], dtype=float)
    x = np.asarray(h["x"], dtype=float)
    if t.size == 0 or x.size == 0:
        raise ValueError("The mass ledger requires a non-empty history.")
    if t.size != x.size:
        raise ValueError("Ledger history arrays t and x must have equal length.")

    igniter_commanded = max(c.ign_mdot, 0.0) * max(c.ign_time, 0.0)
    igniter_history = np.asarray(h["m_igniter"], dtype=float)
    if igniter_history.shape != t.shape:
        raise ValueError(
            "Ledger history array m_igniter must match the time array."
        )
    igniter_injected = float(igniter_history[-1])

    initial_chamber = float(res["m_total0"])
    generated_total = float(np.asarray(h["m_gen"], dtype=float)[-1])
    generated_propellant = generated_total - igniter_injected
    discharged_total = float(np.asarray(h["m_out"], dtype=float)[-1])
    final_chamber = float(np.asarray(h["m_total_eos"], dtype=float)[-1])

    initial_propellant = propellant_mass(c)
    xw = web_thickness(c)
    x_final = min(max(float(x[-1]), 0.0), xw)
    Ri = min(c.R_i0 + x_final, c.R_p)
    Lp = max(c.L_p0 - c.n_end * x_final, 0.0)
    remaining_volume = (
        math.pi
        * max(c.R_p * c.R_p - Ri * Ri, 0.0)
        * Lp
    )
    final_unburned = c.rho_p * remaining_volume
    burned_geometry = initial_propellant - final_unburned

    total_signed = (
        initial_chamber
        + generated_total
        - discharged_total
        - final_chamber
    )
    total_scale = max(
        abs(initial_chamber + generated_total),
        tiny,
    )
    total_relative = abs(total_signed) / total_scale

    propellant_signed = generated_propellant - burned_geometry
    propellant_scale = max(abs(initial_propellant), tiny)
    propellant_relative = abs(propellant_signed) / propellant_scale

    mdot_total = np.asarray(h["mdot_out_total"], dtype=float)
    mdot_gas = np.asarray(h["mdot_out_gas"], dtype=float)
    mdot_condensed = np.asarray(
        h["mdot_out_condensed"],
        dtype=float,
    )
    expected_shape = t.shape
    for name, values in (
        ("mdot_out_total", mdot_total),
        ("mdot_out_gas", mdot_gas),
        ("mdot_out_condensed", mdot_condensed),
    ):
        if values.shape != expected_shape:
            raise ValueError(
                f"Ledger history array {name} must match the time array."
            )

    total_quadrature = float(_TRAPZ(mdot_total, t))
    gas_quadrature = float(_TRAPZ(mdot_gas, t))
    condensed_quadrature = float(_TRAPZ(mdot_condensed, t))
    quadrature_scale = max(abs(discharged_total), tiny)
    quadrature_mismatch = (
        abs(total_quadrature - discharged_total) / quadrature_scale
    )

    total_passed = bool(
        total_relative <= MASS_LEDGER_TOTAL_RELATIVE_TOLERANCE
    )
    propellant_passed = bool(
        propellant_relative <= propellant_tol
    )
    mass_passed = bool(total_passed and propellant_passed)

    mass = {
        "basis": "total_product_mass",
        "initial_chamber_inventory_kg": initial_chamber,
        "generated_mass_total_kg": generated_total,
        "generated_mass_propellant_kg": generated_propellant,
        "igniter_mass_commanded_kg": igniter_commanded,
        "igniter_mass_injected_kg": igniter_injected,
        "discharged_mass_total_kg": discharged_total,
        "final_chamber_inventory_kg": final_chamber,
        "initial_propellant_mass_kg": initial_propellant,
        "final_unburned_propellant_kg": final_unburned,
        "burned_propellant_geometry_kg": burned_geometry,
        "total_system": {
            "signed_residual_kg": total_signed,
            "relative_residual": total_relative,
            "tolerance": MASS_LEDGER_TOTAL_RELATIVE_TOLERANCE,
            "passed": total_passed,
        },
        "propellant": {
            "signed_residual_kg": propellant_signed,
            "relative_residual": propellant_relative,
            "tolerance": propellant_tol,
            "tolerance_formula": "max(1e-9, 100*rtol)",
            "rtol": float(c.rtol),
            "passed": propellant_passed,
        },
        "phase_discharge_diagnostic": {
            "method": "sampled_history_trapezoid",
            "grid_dependent": True,
            "gas_mass_kg": gas_quadrature,
            "condensed_mass_kg": condensed_quadrature,
            "total_quadrature_mass_kg": total_quadrature,
            "ode_total_mass_kg": discharged_total,
            "relative_quadrature_mismatch": quadrature_mismatch,
        },
        "passed": mass_passed,
    }

    energy = {
        "available": False,
        "closure_residual_available": False,
        "reason": (
            "No independent energy state is integrated by the current "
            "GRIBS model."
        ),
        "model": {
            "equilibrium_constraint": "HP (assigned enthalpy and pressure)",
            "temperature_treatment": (
                "T0 is the CEA equilibrium temperature multiplied by eta_T0."
            ),
            "psi_definition": "Psi = R_eff * T0",
            "state_relation": "p0 * Vg = m_total * Psi",
            "state_relation_is_energy_closure": False,
        },
    }

    return {
        "schema_version": "a4-ledger-1",
        "mass": mass,
        "energy": energy,
        "passed": mass_passed,
    }
