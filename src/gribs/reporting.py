"""File output and text reporting (A-6 unit 7). Moved verbatim from section 11 of the entry (the time-history CSV
columns and units, write_csv, the nozzle-transition summary lines, summary_text), followed by the CEA
theoretical diagnostics display helpers of section 13 (cea_diagnostics_status, _cea_diag_display and their
field lists), which the text summary and the figure read (U-35 placement).

Tests patch summary_text and write_csv here, e.g. "gribs.reporting.write_csv"; the entry re-exports every
name and reads them as gribs.reporting.<name>, so a patch on the entry would not reach them.
"""
from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Any, List, Optional, Sequence, Tuple

import numpy as np

from gribs.config import Config
from gribs.constants import PROGRAM_VERSION, SCHEMA_VERSION


#: v0.5.0-alpha column set.  MEANING CHANGES vs v0.4.3-alpha are RENAMED, not
#: silently reused:  mdot_gen -> mdot_gen_total, mdot_out -> mdot_out_total,
#: m_gas_eos -> m_total_eos, m_gas_bal -> m_total_balance (TOTAL product mass
#: basis: gas + condensed).  m_gen / m_out are the integrated TOTAL masses.
CSV_COLUMNS = ["t", "phase", "p0", "x", "Ri", "Lp", "Ab", "Vg", "Kn", "Rt", "At",
               "r", "mdot_gen_total", "mdot_out_total", "mdot_out_gas",
               "mdot_out_condensed", "m_gen", "m_igniter", "m_out",
               "m_total_eos",
               "m_total_balance", "m_gas_equilibrium", "m_condensed",
               "Y_gas", "Y_condensed", "Psi",
               "regime", "Me", "pe", "ve", "F", "impulse",
               "cstar", "CF", "R", "T0", "gamma", "extrap"]
CSV_UNITS = {"t": "s", "p0": "Pa", "x": "m", "Ri": "m", "Lp": "m", "Ab": "m2",
             "Vg": "m3", "Kn": "-", "Rt": "m", "At": "m2", "r": "m/s",
             "mdot_gen_total": "kg/s", "mdot_out_total": "kg/s",
             "mdot_out_gas": "kg/s", "mdot_out_condensed": "kg/s",
             "m_gen": "kg", "m_igniter": "kg", "m_out": "kg",
             "m_total_eos": "kg", "m_total_balance": "kg",
             "m_gas_equilibrium": "kg", "m_condensed": "kg",
             "Y_gas": "-", "Y_condensed": "-", "Psi": "J/kg",
             "Me": "-", "pe": "Pa", "ve": "m/s", "F": "N", "impulse": "N.s",
             "cstar": "m/s", "CF": "-",
             "R": "J/(kg.K)", "T0": "K", "gamma": "-"}


def write_csv(h: dict, path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([f"{k}[{CSV_UNITS[k]}]" if k in CSV_UNITS else k
                    for k in CSV_COLUMNS])
        n = len(h["t"])
        for i in range(n):
            row = []
            for k in CSV_COLUMNS:
                v = h[k][i]
                row.append(v if isinstance(v, (str, np.str_, bool, np.bool_))
                           else f"{float(v):.9e}")
            w.writerow(row)


def nozzle_transition_summary_lines(
    transition: dict,
    warnings: Sequence[dict],
) -> List[str]:
    """Format the nozzle-transition diagnostic for the text summary."""
    lines = [
        "[ nozzle transition diagnostic ]",
        f"  policy                            : {transition['policy']}",
        f"  applicable                        : {transition['applicable']}",
        f"  transition band detected          : {transition['detected']}",
        f"  contiguous intervals              : {transition['interval_count']}",
        (
            "  residence time                    : "
            f"{1.0e3 * transition['duration_s']:.4f} ms"
        ),
        (
            "  fraction of computed duration     : "
            f"{100.0 * transition['duration_fraction']:.6f} %"
        ),
        (
            "  model impulse in transition band  : "
            f"{transition['impulse_Ns']:.6f} N.s"
        ),
        (
            "  fraction of total model impulse   : "
            f"{100.0 * transition['impulse_fraction']:.6f} %"
        ),
    ]

    for warning in warnings:
        if warning.get("code") == "W_NOZZLE_TRANSITION":
            lines.append(f"  WARNING [{warning['code']}]")
            lines.append(f"    {warning['message']}")

    return lines


def summary_text(s: dict, c: Config, thermo_md: Optional[dict] = None,
                 cea_diag: Optional[dict] = None,
                 migration_notes: Optional[Sequence[str]] = None) -> str:
    L = []
    A = L.append
    A("=" * 78)
    A("SOLID ROCKET MOTOR - INTERNAL BALLISTICS SUMMARY")
    A(f"GRIBS v{PROGRAM_VERSION}   (schema {SCHEMA_VERSION})")
    A("=" * 78)
    A("")
    if migration_notes:
        A("[ configuration migration ]")
        for note in migration_notes:
            A(f"  WARNING: {note}")
        A("")
    A("[ configuration ]")
    A(f"  propellant      rho_p = {c.rho_p:.2f} kg/m3,  "
      f"r = {c.a_burn*1e3:.4g} mm/s x (p0/{c.p_ref/1e6:g} MPa)^{c.n_burn:g}")
    A(f"                  sigma_p = {c.sigma_p:g} 1/K, T_grain = {c.T_grain:g} K, "
      f"erosive alpha = {c.ero_alpha:g}")
    A(f"  grain           Ri0 = {c.R_i0*1e3:.3f} mm, Rp = {c.R_p*1e3:.3f} mm, "
      f"Lp0 = {c.L_p0*1e3:.3f} mm, end faces = {c.n_end:g}")
    A(f"                  web = {s['web_mm']:.3f} mm, propellant mass = "
      f"{s['propellant_mass_g']:.3f} g, Vg0 = {c.V_g0*1e6:.3f} cm3")
    A(f"  nozzle          Rt0 = {c.R_t0*1e3:.3f} mm, At0 = {math.pi*c.R_t0**2:.4e} m2, "
      f"Ae/At = {c.eps_nozzle:g}, Cd = {c.Cd:g}, eta_F = {c.eta_thrust:g}")
    A(f"  environment     pa = {c.p_a/1e3:.3f} kPa, p0(0) = {c.p0_init/1e3:.3f} kPa")
    A("")
    A("[ thermochemistry ]")
    A(f"  backend                          : {s['thermochemistry_backend']}")
    if thermo_md:
        A(f"  property range                   : "
          f"{s['property_range_min_Pa']/1e5:.6g} - {s['property_range_max_Pa']/1e5:.6g} bar "
          f"({s['property_range_min_Pa']:.6g} - {s['property_range_max_Pa']:.6g} Pa)")
        if thermo_md.get("cea_python_version"):
            A(f"  official CEA package version     : {thermo_md['cea_python_version']} "
              f"(library {thermo_md.get('cea_library_version')})")
            A(f"  official CEA module path         : {thermo_md.get('cea_python_module_path')}")
            A(f"  equilibrium formulation          : "
              f"{thermo_md.get('equilibrium_constraint')}")
        if thermo_md.get("compatibility_backend"):
            A(f"  NOTE                             : compatibility backend "
              f"({thermo_md.get('warning')})")
        A(f"  ideal mixture density            : "
          f"{thermo_md.get('ideal_mixture_density_kg_m3'):.6f} kg/m3")
        A(f"  packing fraction                 : {thermo_md.get('packing_fraction'):.6f}")
        A(f"  bulk propellant density          : "
          f"{thermo_md.get('bulk_propellant_density_kg_m3'):.6f} kg/m3")
        A(f"  cache key                        : {thermo_md.get('cache_key')}")
        A(f"  cache file                       : {thermo_md.get('cache_file')}")
        A(f"  cache reused                     : {thermo_md.get('cache_used')}")
    A(f"  property policy                  : {c.property_policy}")
    A(f"  values outside property range    : {s['property_extrapolation']} "
      f"({s['property_extrapolation_points']} sample points, "
      f"{100.0*s['property_extrapolation_fraction']:.3f} %)")
    A(f"  pressure range visited           : {s['pressure_min_Pa']:.6g} - "
      f"{s['pressure_max_Pa']:.6g} Pa")
    cea_status, cea_missing, cea_invalid = cea_diagnostics_status(cea_diag)
    if cea_status == "ok":
        A("  CEA theoretical comparison only (never used by the GRIBS model):")
        A(f"    cea_theoretical_cstar_m_s      : {cea_diag['cea_theoretical_cstar_m_s']:.3f}")
        A(f"    cea_theoretical_cf             : {cea_diag['cea_theoretical_cf']:.5f}")
        A(f"    cea_theoretical_isp_s          : {cea_diag['cea_theoretical_isp_s']:.4f}")
        A(f"    cea_theoretical_isp_vacuum_s   : {cea_diag['cea_theoretical_isp_vacuum_s']:.4f}")
        A(f"    evaluation                     : "
          f"{cea_diag['cea_theoretical_exit_definition']} at "
          f"p_c = {cea_diag['cea_theoretical_chamber_pressure_Pa']/1e6:.6f} MPa")
    elif cea_status != "not_requested":
        A("  CEA theoretical comparison only (never used by the GRIBS model):")
        A(f"    status                         : {cea_status} "
          "(no CEA comparison value is substituted)")
        if cea_status == "failed":
            A(f"    error                          : {cea_diag.get('cea_theoretical_error')}")
        else:
            if cea_missing:
                A(f"    missing fields                 : {', '.join(cea_missing)}")
            if cea_invalid:
                A(f"    invalid fields                 : {', '.join(cea_invalid)}")
            if isinstance(cea_diag, dict):
                for _key, _fmt in (("cea_theoretical_cstar_m_s", ".3f"),
                                   ("cea_theoretical_cf", ".5f"),
                                   ("cea_theoretical_isp_s", ".4f"),
                                   ("cea_theoretical_isp_vacuum_s", ".4f")):
                    A(f"    {_key:<31}: "
                      f"{_cea_diag_display(cea_diag, _key, _fmt, cea_invalid)}")
    if cea_status != "not_requested":
        A(f"    GRIBS effective c*             : {s['cstar_eff_m_s']:.3f} m/s")
        A(f"    GRIBS effective CF             : {s['CF_eff']:.5f}")
        A(f"    GRIBS effective Isp            : {s['Isp_s']:.4f} s")
    if cea_status == "ok":
        A("    note: both sets are reported side by side for context.  The CEA")
        A("          station values follow the conventions of the official package")
        A("          (throat c* is the theoretical characteristic velocity and the")
        A("          throat Cf equals v_e/c*), while cstar_eff_m_s, CF_eff and Isp_s")
        A("          are GRIBS results that include the pressure-thrust term of the")
        A("          unsteady nozzle model; they are not expected to be equal and are")
        A("          never substituted for one another.")
    A("")
    A("[ two-phase model (homogeneous equilibrium, model B) ]")
    tpm = s["two_phase_model"]
    A(f"  model                            : {tpm['mode']}"
      f"{'  (ACTIVE)' if tpm['active'] else '  (phase split inactive: Yg := 1)'}")
    A(f"  chamber EOS                      : p0*Vg = mt*Yg(p)*Rg(p)*T0(p) = mt*Psi(p)")
    A(f"  conserved chamber mass           : TOTAL product mass mt = m_gas + m_condensed")
    A(f"  max condensed fraction (run)     : "
      f"{s['maximum_condensed_mass_fraction_during_run']:.6f} [-]")
    A(f"  mean condensed fraction (burn)   : "
      f"{s['mean_condensed_mass_fraction_during_burn']:.6f} [-]")
    A(f"  condensed fraction at MEOP       : "
      f"{s['condensed_mass_fraction_at_max_pressure']:.6f} [-]")
    A(f"  Yc over the property table       : "
      f"{s['condensed_mass_fraction_table_min']:.6f} - "
      f"{s['condensed_mass_fraction_table_max']:.6f} [-]")
    A(f"  min (1 - p*Psi'/Psi) on table    : "
      f"{s['min_pressure_factor_1_minus_pPsi_over_Psi']:.6f} [-]  (must be > 0)")
    A(f"  nozzle entrainment assumption    : {s['nozzle_entrainment_assumption']}")
    A(f"  condensed volume assumption      : {s['condensed_volume_assumption']}")
    A("  NOTE: particle slip is NOT modelled; the complete-entrainment momentum")
    A("        term can overestimate the thrust for propellants with a large")
    A("        condensed fraction.  eta_F may additionally lump two-phase losses.")
    A("        These results are NOT experimentally validated - independent")
    A("        validation is required before any engineering use.")
    A("")
    A("[ initial balance ]")
    A(f"  initial nozzle regime            : {s['initial_regime']}")
    if s["choke_limit_pressure_Pa"]:
        A(f"  choking-limit chamber pressure   : {s['choke_limit_pressure_Pa']/1e3:.3f} kPa")
    A(f"  mdot_gen(0) = {s['initial_mdot_gen_kg_s']:.6e} kg/s")
    A(f"  mdot_out(0) = {s['initial_mdot_out_kg_s']:.6e} kg/s"
      f"   ->  chamber { 'fills (pressure rises)' if s['initial_mdot_gen_kg_s'] > s['initial_mdot_out_kg_s'] else 'empties (pressure falls)'}")
    A("")
    A("[ events ]")
    for name, times in s["event_times"].items():
        if times:
            A(f"  {name:<10s}: " + ", ".join(f"{v:.6f} s" for v in times[:6]))
    A(f"  burning phase terminated by      : {s['burn_stop_reason']}")
    A(f"  blowdown phase terminated by     : {s['blowdown_stop_reason']}")
    A("")
    A("[ performance ]")
    A(f"  web burn time                    : {s['burn_time_s']:.6f} s")
    A(f"  total computed duration          : {s['total_time_s']:.6f} s")
    A(f"  max chamber pressure (MEOP)      : {s['p_max_Pa']/1e6:.6f} MPa "
      f"at t = {s['t_p_max_s']:.6f} s")
    A(f"  min chamber pressure (burning)   : {s['p_min_burn_Pa']/1e6:.6f} MPa "
      f"at t = {s['t_p_min_s']:.6f} s")
    A(f"  mean chamber pressure (burning)  : {s['p_mean_burn_Pa']/1e6:.6f} MPa")
    A(f"  max thrust                       : {s['F_max_N']:.4f} N "
      f"at t = {s['t_F_max_s']:.6f} s")
    A(f"  mean thrust (burning)            : {s['F_mean_burn_N']:.4f} N")
    A(f"  impulse, burning phase           : {s['impulse_burn_Ns']:.4f} N.s")
    A(f"  impulse, blowdown                : {s['impulse_blowdown_Ns']:.4f} N.s")
    A(f"  total impulse                    : {s['impulse_total_Ns']:.4f} N.s")
    A(f"  effective specific impulse       : {s['Isp_s']:.2f} s")
    A(f"  effective c*                     : {s['cstar_eff_m_s']:.2f} m/s")
    A(f"  effective thrust coefficient CF  : {s['CF_eff']:.5f}")
    A(f"  Kn = Ab/At  (start/max/end)      : {s['Kn_initial']:.1f} / "
      f"{s['Kn_max']:.1f} / {s['Kn_final']:.1f}")
    A(f"  final throat radius              : {s['throat_radius_final_mm']:.5f} mm")
    A(f"  nozzle regimes visited           : {', '.join(s['regimes_visited'])}")
    A("")
    for line in nozzle_transition_summary_lines(
        s["nozzle_transition"],
        s["warnings"],
    ):
        A(line)
    A("")
    run_warnings = [
        warning
        for warning in s["warnings"]
        if warning.get("code") in {
            "W_POST_RUN_CONSISTENCY",
            "W_MASS_LEDGER_FAILED",
            "W_HIGH_CONDENSATE_MODEL_LIMITATION",
            "W_HIGH_CONDENSATE_G4B_UNVALIDATED",
            "W_TWO_PHASE_NOZZLE_APPROXIMATION",
            "W_THERMAL_EFFICIENCY_PARTIAL",
            "E_OUT_OF_VALIDATED_RANGE",
        }
    ]
    if run_warnings:
        A("[ verification warnings ]")
        for warning in run_warnings:
            A(
                f"  WARNING [{warning['code']}] "
                f"severity={warning['severity']}"
            )
            A(f"    {warning['message']}")
        A("")

    A("[ mass ledger ]")
    ledger = s["ledger"]
    mass = ledger["mass"]
    total = mass["total_system"]
    propellant = mass["propellant"]
    phase = mass["phase_discharge_diagnostic"]
    A(f"  basis                            : {mass['basis']}")
    A(f"  initial chamber inventory        : "
      f"{mass['initial_chamber_inventory_kg']:.9e} kg")
    A(f"  generated products, total        : "
      f"{mass['generated_mass_total_kg']:.9e} kg")
    A(f"    propellant-derived             : "
      f"{mass['generated_mass_propellant_kg']:.9e} kg")
    A(f"    igniter-derived, injected      : "
      f"{mass['igniter_mass_injected_kg']:.9e} kg")
    A(f"    igniter charge, commanded      : "
      f"{mass['igniter_mass_commanded_kg']:.9e} kg")
    A(f"  discharged products, ODE state   : "
      f"{mass['discharged_mass_total_kg']:.9e} kg")
    A(f"  final chamber inventory, EOS     : "
      f"{mass['final_chamber_inventory_kg']:.9e} kg")
    A(f"  initial solid propellant         : "
      f"{mass['initial_propellant_mass_kg']:.9e} kg")
    A(f"  final unburned propellant        : "
      f"{mass['final_unburned_propellant_kg']:.9e} kg")
    A(f"  burned propellant, geometry      : "
      f"{mass['burned_propellant_geometry_kg']:.9e} kg")
    A(f"  total-system signed residual     : "
      f"{total['signed_residual_kg']:+.9e} kg")
    A(f"  total-system relative residual   : "
      f"{total['relative_residual']:.3e} "
      f"(limit {total['tolerance']:.1e}, "
      f"{'PASS' if total['passed'] else 'FAIL'})")
    A(f"  propellant signed residual       : "
      f"{propellant['signed_residual_kg']:+.9e} kg")
    A(f"  propellant relative residual     : "
      f"{propellant['relative_residual']:.3e} "
      f"(limit {propellant['tolerance']:.1e}, "
      f"{'PASS' if propellant['passed'] else 'FAIL'})")
    A(f"  phase-discharge diagnostic       : {phase['method']} "
      f"(grid dependent)")
    A(f"    gas / condensed                : "
      f"{phase['gas_mass_kg']:.9e} / "
      f"{phase['condensed_mass_kg']:.9e} kg")
    A(f"    quadrature mismatch            : "
      f"{phase['relative_quadrature_mismatch']:.3e}")
    A(f"  mass-ledger result               : "
      f"{'PASS' if mass['passed'] else 'FAIL'}")
    A("")
    A("[ energy ledger ]")
    energy = ledger["energy"]
    A(f"  closure residual available       : "
      f"{energy['closure_residual_available']}")
    A(f"  reason                           : {energy['reason']}")
    A(f"  equilibrium constraint           : "
      f"{energy['model']['equilibrium_constraint']}")
    A(f"  temperature treatment            : "
      f"{energy['model']['temperature_treatment']}")
    A(f"  state relation                   : "
      f"{energy['model']['state_relation']}")
    A("  NOTE: the state relation is not a complete energy-conservation law.")
    A("")
    A("[ verification ]")
    A(f"  legacy propellant balance metric : "
      f"{s['propellant_mass_balance_error']:.3e} "
      f"(includes igniter in the legacy numerator)")
    A(f"  EOS-history invariant error      : "
      f"{s['total_mass_consistency_error_invariant']:.3e}")
    A(f"  invariant reference mass         : "
      f"{s['total_mass_consistency_reference_mass_kg']:.9e} kg")
    A(f"  invariant metric version         : "
      f"{s['total_mass_consistency_metric_version']}")
    A(f"  legacy shrinking-mass error      : "
      f"{s['total_mass_consistency_error']:.3e} (deprecated)")
    A(f"  property-range extrapolation used: {s['property_extrapolation']}"
      f"   (policy = {c.property_policy})")
    if s["property_extrapolation"] and c.property_policy == "extrapolate":
        A("  WARNING: results rely on extrapolated thermochemical properties.")
    # --- A-5.1-5: validation status and c* indicators ---
    val = s.get("validation", {})
    if val:
        A("")
        A("[ validation status ]")
        A(f"  semantic status                : {val.get('semantic_status', 'n/a')}")
        A(f"  numerical status               : {val.get('numerical_status', 'n/a')}")
        A(f"  model status                   : {val.get('model_status', 'n/a')}")
        A(f"  table pressure range [not validated] : {val.get('table_pressure_range_Pa', 'n/a')}")
        A(f"  table Yg range [not validated]       : {val.get('table_Yg_range', 'n/a')}")
        A(f"  mass basis semantics version   : {val.get('mass_basis_semantics_version', 'n/a')}")
    lmc = s.get("local_model_cstar_m_s")
    if lmc is not None and not (isinstance(lmc, float) and math.isnan(lmc)):
        A(f"  local model c* (burn median)   : {lmc:.2f} m/s")
    A("=" * 78)
    return "\n".join(L) + "\n"


# Fields of the CEA theoretical diagnostics that the text summary and the
# figure display.  Numeric fields must be finite real numbers; the exit
# definition must be a string.
_CEA_DIAG_DISPLAY_NUMERIC_FIELDS = (
    "cea_theoretical_cstar_m_s",
    "cea_theoretical_cf",
    "cea_theoretical_isp_s",
    "cea_theoretical_isp_vacuum_s",
    "cea_theoretical_chamber_pressure_Pa",
)
_CEA_DIAG_DISPLAY_TEXT_FIELDS = ("cea_theoretical_exit_definition",)


def _cea_diag_value_ok(key: str, value: Any) -> bool:
    if key in _CEA_DIAG_DISPLAY_TEXT_FIELDS:
        return isinstance(value, str)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return False
    try:
        return math.isfinite(float(value))
    except (OverflowError, ValueError, TypeError):
        return False


def cea_diagnostics_status(cea_diag: Optional[dict]) -> Tuple[str, List[str], List[str]]:
    """Classify CEA theoretical diagnostics for the output writers.

    Returns ``(status, missing_fields, invalid_fields)`` with status
      'not_requested' : None (diagnostics disabled or not applicable),
      'failed'        : the diagnostic computation raised (error recorded),
      'ok'            : every displayed field is present and valid, and no
                        other cea_theoretical_* number is non-finite,
      'incomplete'    : displayed fields are missing (none invalid),
      'invalid'       : a field is non-finite / of the wrong type.
    Never used by the GRIBS model or the run verdict.
    """
    if cea_diag is None:
        return "not_requested", [], []
    if not isinstance(cea_diag, dict):
        return "invalid", [], ["<not a mapping>"]
    if "cea_theoretical_error" in cea_diag:
        return "failed", [], []
    display = _CEA_DIAG_DISPLAY_NUMERIC_FIELDS + _CEA_DIAG_DISPLAY_TEXT_FIELDS
    missing = [k for k in display if k not in cea_diag]
    invalid = [k for k in display if k in cea_diag and not _cea_diag_value_ok(k, cea_diag[k])]
    for k, v in cea_diag.items():
        if (k not in display and isinstance(v, float) and not math.isfinite(v)):
            invalid.append(k)
    if invalid:
        return "invalid", missing, invalid
    if missing:
        return "incomplete", missing, invalid
    return "ok", [], []


def _cea_diag_display(cea_diag: dict, key: str, fmt: str, invalid: Sequence[str]) -> str:
    """Formatted value, or an explicit marker; never a substituted number."""
    if key not in cea_diag:
        return "unavailable (missing)"
    if key in invalid:
        return f"invalid ({cea_diag[key]!r})"
    return format(cea_diag[key], fmt)
