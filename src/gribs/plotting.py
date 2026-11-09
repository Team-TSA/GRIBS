"""Plotting (A-6 unit 7). Moved verbatim from section 10 of the entry: the one result figure (make_figure), its
colours and regime styles, and the key-results block with its automatic font fit.

Tests patch make_figure here ("gribs.plotting.make_figure"); the entry re-exports every name and reads
make_figure as gribs.plotting.make_figure, so a patch on the entry would not reach it. Names of other package
modules that tests patch are read as <module>.<name> (ADR decision 1). The matplotlib backend is set to Agg
before pyplot is imported, as the entry does, so the figure does not depend on which module imports pyplot
first.
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.gridspec import GridSpec  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

import gribs.thermo  # noqa: E402
from gribs.config import Config, web_thickness  # noqa: E402
from gribs.constants import PROGRAM_VERSION  # noqa: E402
from gribs.reporting import cea_diagnostics_status  # noqa: E402


C_PRESS = "#14507d"
C_THRUST = "#b3331f"
C_GEN = "#2e8b57"
C_OUT = "#7b4fa8"
C_ACC = "#d98c00"
C_GREY = "#5a5a5a"
REGIME_FACE = {"choked": "#dce9f5", "subsonic": "#fdecd2",
               "separated": "#ece2f7", "no-flow": "#ededed"}
REGIME_NAME = {"choked": "choked (sonic throat)", "subsonic": "subsonic nozzle",
               "separated": "separated (overexpanded)", "no-flow": "no outflow"}


def _style() -> None:
    plt.rcParams.update({
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "axes.facecolor": "#fdfdfe",
        "axes.edgecolor": "#54606b",
        "axes.linewidth": 0.9,
        "axes.grid": True,
        "axes.axisbelow": True,
        "grid.color": "#c5ccd3",
        "grid.linewidth": 0.55,
        "grid.alpha": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "font.family": "DejaVu Sans",
        "font.size": 9.0,
        "axes.titlesize": 10.0,
        "axes.titleweight": "bold",
        "axes.labelsize": 9.0,
        "legend.fontsize": 7.8,
        "legend.frameon": True,
        "legend.framealpha": 0.92,
        "legend.edgecolor": "#b9c2ca",
        "xtick.labelsize": 8.2,
        "ytick.labelsize": 8.2,
        "xtick.direction": "out",
        "ytick.direction": "out",
    })


def _regime_spans(t, regime):
    spans, i0 = [], 0
    for i in range(1, len(t)):
        if regime[i] != regime[i0]:
            spans.append((t[i0], t[i], regime[i0]))
            i0 = i
    spans.append((t[i0], t[-1], regime[i0]))
    return spans


def _shade(ax, spans, t_lo=None, t_hi=None, scale=1.0):
    for a, b, key in spans:
        if t_lo is not None:
            a, b = max(a, t_lo), min(b, t_hi)
        if b <= a:
            continue
        ax.axvspan(a * scale, b * scale, color=REGIME_FACE.get(key, "#eeeeee"),
                   alpha=0.85, lw=0, zorder=0)


#: Style of the key-results block and the geometry of its automatic fit.
#: The block is sized to its own axes: it grows up to KEY_RESULTS_MAX_FONTSIZE for
#: readability and shrinks (with a warning) if a configuration produces so many
#: lines that it would otherwise run over the neighbouring panels.
KEY_RESULTS_MAX_FONTSIZE = 9.6      # never grow beyond this (readability)
KEY_RESULTS_MIN_FONTSIZE = 5.2      # never shrink below this (a warning is printed)
KEY_RESULTS_LINESPACING = 1.28      # multiple of the font size
KEY_RESULTS_BOX_PAD = 0.55          # bbox pad, in font-size units
KEY_RESULTS_CHAR_WIDTH = 0.602      # DejaVu Sans Mono advance width, in em


def key_results_fontsize(fig, ax, n_lines: int, longest_line: int) -> float:
    """Point size whose text box just fits inside ``ax`` (clamped to sane limits).

    The box size is estimated analytically (line count, monospace advance width
    and the ``bbox`` padding) from the axes geometry, so the fit is decided
    before anything is drawn and the key-results block can never spill over the
    neighbouring panels - whatever the configuration produces (extra event
    lines, CEA diagnostics, long backend names, ...).  A negative return value
    means "does not fit even at the minimum size".
    """
    width_pt = fig.get_size_inches()[0] * ax.get_position().width * 72.0
    height_pt = fig.get_size_inches()[1] * ax.get_position().height * 72.0
    fs_h = height_pt / (max(1, n_lines) * KEY_RESULTS_LINESPACING + 2.0 * KEY_RESULTS_BOX_PAD + 0.8)
    fs_w = width_pt / (max(1, longest_line) * KEY_RESULTS_CHAR_WIDTH + 2.0 * KEY_RESULTS_BOX_PAD)
    return float(min(KEY_RESULTS_MAX_FONTSIZE, fs_h, fs_w))


def _draw_key_results(fig, ax, lines: Sequence[str]):
    """Draw the key-results block, fitted to (and contained in) its own axes."""
    text = "\n".join(lines)
    fontsize = key_results_fontsize(fig, ax, len(lines), max(len(l) for l in lines))
    if fontsize < KEY_RESULTS_MIN_FONTSIZE:
        print(f"WARNING: the key-results block needs a font size of "
              f"{fontsize:.2f} pt to fit its panel; clamping to "
              f"{KEY_RESULTS_MIN_FONTSIZE:.2f} pt (the block may overflow).")
        fontsize = KEY_RESULTS_MIN_FONTSIZE
    return ax.text(0.0, 1.0, text, transform=ax.transAxes, va="top", ha="left",
                   family="monospace", fontsize=fontsize,
                   linespacing=KEY_RESULTS_LINESPACING,
                   bbox=dict(boxstyle=f"round,pad={KEY_RESULTS_BOX_PAD}",
                             fc="#f4f7fa", ec="#a9bbcc", lw=1.0))


def make_figure(res, h, s, c: Config, path: Path, cea_diag: Optional[dict] = None) -> None:
    _style()
    t = h["t"]
    burn = h["phase"] == "burn"
    t_bo = s["burn_time_s"]
    spans = _regime_spans(t, h["regime"])
    has_blow = res["sol_blow"] is not None

    # The key-results box needs much more vertical room than the tallest plot
    # panel, so it gets its own full-height column on the right-hand side
    # (gs[:, 2]) instead of sharing a cell with a plot panel.  With the previous
    # 3x3 layout the text ran over the panel below it.
    fig = plt.figure(figsize=(16.6, 15.0))
    gs = GridSpec(4, 3, figure=fig, height_ratios=[1.45, 1.0, 1.0, 1.0],
                  hspace=0.40, wspace=0.28,
                  left=0.055, right=0.975, top=0.940, bottom=0.072)

    ax_p = fig.add_subplot(gs[0, 0:2])       # (a) chamber pressure (wide)
    ax_key = fig.add_subplot(gs[:, 2])       # key results (full-height column)
    ax_F = fig.add_subplot(gs[1, 0])         # (b) thrust
    ax_md = fig.add_subplot(gs[2, 0])        # (c) mass-flow balance
    ax_geo = fig.add_subplot(gs[3, 0])       # (d) grain geometry evolution
    ax_zm = fig.add_subplot(gs[1, 1])        # (e) ignition transient
    ax_nz = fig.add_subplot(gs[2, 1])        # (f) nozzle operating point
    ax_r = fig.add_subplot(gs[3, 1])         # (g) regression of the grain

    # ---------------- (a) chamber pressure ----------------
    _shade(ax_p, spans)
    ax_p.fill_between(t[burn], 0.0, h["p0"][burn] / 1e6, color=C_PRESS, alpha=0.13, lw=0)
    ax_p.plot(t[burn], h["p0"][burn] / 1e6, color=C_PRESS, lw=2.1, label="chamber pressure $p_0$")
    if has_blow:
        ax_p.plot(t[~burn], h["p0"][~burn] / 1e6, color=C_PRESS, lw=1.3, ls=":",
                  label="blowdown (post-burnout)")
    if np.isfinite(h["p_eq"]).any():
        mk = np.isfinite(h["p_eq"])
        ax_p.plot(t[mk], h["p_eq"][mk] / 1e6, color=C_GREY, lw=1.1, ls="--",
                  label="quasi-steady equilibrium")
    if s["choke_limit_pressure_Pa"]:
        ax_p.axhline(s["choke_limit_pressure_Pa"] / 1e6, color="#7a7a7a", lw=0.9, ls="-.")
        ax_p.text(0.995, s["choke_limit_pressure_Pa"] / 1e6, "  choking limit",
                  transform=ax_p.get_yaxis_transform(), ha="right", va="bottom",
                  fontsize=7.4, color="#5a5a5a")
    if s["property_extrapolation"]:
        ex = h["extrap"]
        ax_p.plot(t[ex], h["p0"][ex] / 1e6, color="#c0392b", lw=2.6, alpha=0.55,
                  solid_capstyle="butt", label="outside property-table range")
    ax_p.axvline(t_bo, color="#404040", lw=1.0, ls="--")
    ax_p.annotate("burnout", xy=(t_bo, ax_p.get_ylim()[1]),
                  xytext=(-4, -10), textcoords="offset points",
                  ha="right", va="top", fontsize=8, color="#404040", rotation=90)
    ax_p.plot([s["t_p_max_s"]], [s["p_max_Pa"] / 1e6], "o", ms=5.0,
              mfc="white", mec=C_PRESS, mew=1.6, zorder=5)
    ax_p.annotate(f"MEOP {s['p_max_Pa']/1e6:.3f} MPa @ {s['t_p_max_s']:.3f} s",
                  xy=(s["t_p_max_s"], s["p_max_Pa"] / 1e6),
                  xytext=(-14, 14), textcoords="offset points", ha="right",
                  fontsize=8.2, color=C_PRESS,
                  arrowprops=dict(arrowstyle="->", color=C_PRESS, lw=0.9))
    ax_p.set_xlim(0.0, t[-1])
    ax_p.set_ylim(0.0, max(h["p0"]) / 1e6 * 1.22)
    ax_p.set_xlabel("time  $t$  [s]")
    ax_p.set_ylabel("chamber pressure  $p_0$  [MPa]")
    ax_p.set_title("(a)  Chamber pressure history", loc="left")
    # condensed-phase mass fraction on a second axis (two-phase model B);
    # drawn only when a condensed phase is actually present, so the panel is
    # unchanged for condensed-free propellants
    yc_hist = np.asarray(h["Y_condensed"], dtype=float)
    if float(np.max(yc_hist)) > 0.0:
        ax_yc = ax_p.twinx()
        ax_yc.grid(False)
        ax_yc.plot(t, yc_hist, color="#8d6e63", lw=1.4, ls="-.",
                   label=r"condensed fraction $Y_c$")
        ax_yc.set_ylim(0.0, max(float(np.max(yc_hist)) * 1.35, 1e-6))
        ax_yc.set_ylabel(r"condensed mass fraction  $Y_c$  [-]", color="#8d6e63")
        ax_yc.tick_params(axis="y", colors="#8d6e63")
        ax_yc.spines["right"].set_visible(True)
        hl = ax_p.get_legend_handles_labels()
        h2 = ax_yc.get_legend_handles_labels()
        ax_p.legend(hl[0] + h2[0], hl[1] + h2[1], loc="lower right", ncol=1)
    else:
        ax_p.legend(loc="lower right", ncol=1)

    # ---------------- key results panel ----------------
    ax_key.axis("off")
    ev = s["event_times"]
    lines = [
        "KEY RESULTS",
        "-" * 42,
        f"{'web thickness':<24}{s['web_mm']:>11.3f} mm",
        f"{'propellant density':<24}{s['propellant_density_kg_m3']:>11.2f} kg/m3",
        f"{'propellant mass':<24}{s['propellant_mass_g']:>11.3f} g",
        f"{'web burn time':<24}{s['burn_time_s']:>11.4f} s",
        f"{'total time computed':<24}{s['total_time_s']:>11.4f} s",
        "",
        f"{'max chamber pressure':<24}{s['p_max_Pa']/1e6:>11.4f} MPa",
        f"{'min chamber pressure':<24}{s['p_min_burn_Pa']/1e6:>11.4f} MPa",
        f"{'mean chamber pressure':<24}{s['p_mean_burn_Pa']/1e6:>11.4f} MPa",
        f"{'max thrust':<24}{s['F_max_N']:>11.3f} N",
        f"{'mean thrust (burn)':<24}{s['F_mean_burn_N']:>11.3f} N",
        "",
        f"{'total impulse':<24}{s['impulse_total_Ns']:>11.3f} N s",
        f"{'  of which blowdown':<24}{s['impulse_blowdown_Ns']:>11.3f} N s",
        f"{'effective Isp':<24}{s['Isp_s']:>11.2f} s",
        f"{'effective c*':<24}{s['cstar_eff_m_s']:>11.1f} m/s",
        f"{'effective CF':<24}{s['CF_eff']:>11.4f} -",
        "",
        f"{'Kn  (start / max / end)':<24}",
        f"{'':<10}{s['Kn_initial']:>8.1f} /{s['Kn_max']:>8.1f} /{s['Kn_final']:>8.1f}",
        f"{'nozzle regimes':<24}{', '.join(s['regimes_visited']):>11s}",
        "",
        "THERMOCHEMISTRY",
        f"  backend            {s['thermochemistry_backend']}",
        f"  property range     {s['property_range_min_Pa']/1e5:.3g}"
        f"-{s['property_range_max_Pa']/1e5:.3g} bar",
        f"  two-phase model    {s['two_phase_model']['mode']}",
        f"{'max Yc (run)':<24}{s['maximum_condensed_mass_fraction_during_run']:>11.4f} -",
        f"{'mean Yc (burn)':<24}{s['mean_condensed_mass_fraction_during_burn']:>11.4f} -",
        "",
        "VERIFICATION",
        f"{'propellant mass balance':<24}{s['propellant_mass_balance_error']:>11.2e}",
        f"{'total mass consistency':<24}{s['total_mass_consistency_error']:>11.2e}",
        f"{'property extrapolation':<24}{str(s['property_extrapolation']):>11s}",
    ]
    cea_status, _cea_missing, cea_invalid = cea_diagnostics_status(cea_diag)
    if cea_status == "ok":
        lines += [
            "",
            "CEA theoretical (comparison only)",
            f"  {'c*':<18}{cea_diag['cea_theoretical_cstar_m_s']:>11.1f} m/s",
            f"  {'CF':<18}{cea_diag['cea_theoretical_cf']:>11.4f} -",
            f"  {'Isp':<18}{cea_diag['cea_theoretical_isp_s']:>11.2f} s",
            f"  {'GRIBS c* (effective)':<18}{s['cstar_eff_m_s']:>11.1f} m/s",
        ]
    elif cea_status != "not_requested":
        # Show the diagnostic state; never draw a substituted CEA value.
        def _fig_value(key: str, fmt: str, unit: str) -> str:
            if not isinstance(cea_diag, dict) or key not in cea_diag:
                return f"{'n/a':>11}"
            if key in cea_invalid:
                return f"{'invalid':>11}"
            return f"{format(cea_diag[key], fmt):>11} {unit}"
        lines += [
            "",
            f"CEA theoretical (comparison only): {cea_status.upper()}",
        ]
        if cea_status != "failed":
            lines += [
                f"  {'c*':<18}{_fig_value('cea_theoretical_cstar_m_s', '.1f', 'm/s')}",
                f"  {'CF':<18}{_fig_value('cea_theoretical_cf', '.4f', '-')}",
                f"  {'Isp':<18}{_fig_value('cea_theoretical_isp_s', '.2f', 's')}",
            ]
        lines.append(f"  {'GRIBS c* (effective)':<18}{s['cstar_eff_m_s']:>11.1f} m/s")
    if ev.get("unchoked"):
        lines.append(f"{'first unchoking at':<24}{ev['unchoked'][0]*1e3:>11.3f} ms")
    if ev.get("rechoked"):
        lines.append(f"{'re-choking at':<24}{ev['rechoked'][0]*1e3:>11.3f} ms")
    _draw_key_results(fig, ax_key, lines)

    # ---------------- (b) thrust ----------------
    _shade(ax_F, spans)
    ax_F.fill_between(t[burn], 0.0, h["F"][burn], color=C_THRUST, alpha=0.15, lw=0,
                      label=f"impulse = {s['impulse_burn_Ns']:.2f} N s")
    ax_F.plot(t[burn], h["F"][burn], color=C_THRUST, lw=2.0, label="thrust $F$")
    if has_blow:
        ax_F.plot(t[~burn], h["F"][~burn], color=C_THRUST, lw=1.2, ls=":",
                  label="blowdown")
    ax_F.axvline(t_bo, color="#404040", lw=1.0, ls="--")
    ax_F.plot([s["t_F_max_s"]], [s["F_max_N"]], "o", ms=5, mfc="white",
              mec=C_THRUST, mew=1.6, zorder=5)
    ax_F.annotate(f"{s['F_max_N']:.2f} N", xy=(s["t_F_max_s"], s["F_max_N"]),
                  xytext=(-12, 10), textcoords="offset points", ha="right",
                  fontsize=8.2, color=C_THRUST,
                  arrowprops=dict(arrowstyle="->", color=C_THRUST, lw=0.9))
    ax_F.set_xlim(0.0, t[-1])
    ax_F.set_ylim(0.0, max(h["F"]) * 1.22)
    ax_F.set_xlabel("time  $t$  [s]")
    ax_F.set_ylabel("thrust  $F$  [N]")
    ax_F.set_title("(b)  Thrust history", loc="left")
    ax_F.legend(loc="upper left")

    # ---------------- (c) mass flow balance ----------------
    _shade(ax_md, spans)
    pos = h["mdot_gen_total"] > 0
    ax_md.plot(t[pos], h["mdot_gen_total"][pos] * 1e3, color=C_GEN, lw=1.8,
               label=r"generated (total)  $\dot{m}_{gen}=\rho_p A_b r$")
    ax_md.plot(t, np.maximum(h["mdot_out_total"], 1e-12) * 1e3, color=C_OUT, lw=1.6,
               ls="--", label=r"nozzle (total)  $\dot{m}_{out}$")
    if np.any(np.asarray(h["mdot_out_condensed"]) > 0.0):
        ax_md.plot(t, np.maximum(h["mdot_out_condensed"], 1e-12) * 1e3,
                   color="#8d6e63", lw=1.2, ls=":",
                   label=r"nozzle condensed  $\dot{m}_{out,c}$")
    acc = h["mdot_gen_total"] > h["mdot_out_total"]
    ax_md.fill_between(t, np.maximum(h["mdot_out_total"], 1e-12) * 1e3,
                       np.maximum(h["mdot_gen_total"], 1e-12) * 1e3, where=acc,
                       color=C_ACC, alpha=0.25, lw=0, label="chamber filling")
    ax_md.axvline(t_bo, color="#404040", lw=1.0, ls="--")
    ax_md.set_yscale("log")
    ax_md.set_xlim(0.0, t[-1])
    ax_md.set_xlabel("time  $t$  [s]")
    ax_md.set_ylabel(r"mass flow (total)  [g/s]")
    ax_md.set_title("(c)  Mass-flow balance", loc="left")
    ax_md.legend(loc="lower right")

    # ---------------- (d) burning area / Kn / free volume ----------------
    ax_geo.plot(t[burn], h["Ab"][burn] * 1e4, color="#1b7f79", lw=1.9,
                label="burning area $A_b$")
    ax_geo.set_xlabel("time  $t$  [s]")
    ax_geo.set_ylabel("$A_b$  [cm$^2$]", color="#1b7f79")
    ax_geo.tick_params(axis="y", colors="#1b7f79")
    ax_geo.set_xlim(0.0, t[-1])
    ax_geo.axvline(t_bo, color="#404040", lw=1.0, ls="--")
    ax_g2 = ax_geo.twinx()
    ax_g2.grid(False)
    ax_g2.plot(t[burn], h["Kn"][burn], color="#8a5a00", lw=1.5, ls="--",
               label="$K_n=A_b/A_t$")
    ax_g2.plot(t[burn], h["Vg"][burn] * 1e6, color="#7a7a7a", lw=1.2, ls=":",
               label="free volume $V_g$ [cm$^3$]")
    ax_g2.set_ylabel("$K_n$ [-]   /   $V_g$ [cm$^3$]", color="#8a5a00")
    ax_g2.tick_params(axis="y", colors="#8a5a00")
    ax_g2.spines["right"].set_visible(True)
    hl = ax_geo.get_legend_handles_labels()
    h2 = ax_g2.get_legend_handles_labels()
    ax_geo.legend(hl[0] + h2[0], hl[1] + h2[1], loc="center right")
    ax_geo.set_title("(d)  Grain geometry evolution", loc="left")

    # ---------------- (e) start-up transient ----------------
    R0, T00, g00 = gribs.thermo.gas_props(c.p0_init, c)
    gam = np.sqrt(g00) * (2.0 / (g00 + 1.0)) ** ((g00 + 1.0) / (2.0 * (g00 - 1.0)))
    tau = c.V_g0 * gam / (res["At0"] * math.sqrt(R0 * T00))
    t_zoom = float(min(t[-1], max(25.0 * tau, 3.0e-3)))
    mz = t <= t_zoom + 1e-15
    _shade(ax_zm, spans, 0.0, t_zoom, scale=1e3)
    ax_zm.plot(t[mz] * 1e3, h["p0"][mz] / 1e3, color=C_PRESS, lw=2.0,
               label="$p_0$ [kPa]")
    ax_zm.axhline(c.p_a / 1e3, color=C_GREY, lw=0.9, ls=":", label="ambient $p_a$")
    if s["choke_limit_pressure_Pa"]:
        ax_zm.axhline(s["choke_limit_pressure_Pa"] / 1e3, color="#7a7a7a",
                      lw=0.9, ls="-.", label="choking limit")
    ax_z2 = ax_zm.twinx()
    ax_z2.grid(False)
    ax_z2.plot(t[mz] * 1e3, h["F"][mz], color=C_THRUST, lw=1.4, ls="--",
               label="$F$ [N]")
    ax_z2.set_ylabel("thrust  $F$  [N]", color=C_THRUST)
    ax_z2.tick_params(axis="y", colors=C_THRUST)
    ax_z2.spines["right"].set_visible(True)
    ax_zm.set_xlim(0.0, t_zoom * 1e3)
    ax_zm.set_xlabel("time  $t$  [ms]   (start-up zoom)")
    ax_zm.set_ylabel("$p_0$  [kPa]", color=C_PRESS)
    ax_zm.tick_params(axis="y", colors=C_PRESS)
    ax_zm.set_title(f"(e)  Ignition transient  ($\\tau_{{fill}}\\approx${tau*1e3:.2f} ms)",
                    loc="left")
    hl = ax_zm.get_legend_handles_labels()
    h2 = ax_z2.get_legend_handles_labels()
    ax_zm.legend(hl[0] + h2[0], hl[1] + h2[1], loc="best")

    # ---------------- (f) nozzle diagnostics ----------------
    _shade(ax_nz, spans)
    ax_nz.plot(t, h["Me"], color="#2d6a9f", lw=1.8, label="exit Mach $M_e$")
    ax_nz.axhline(1.0, color="#7a7a7a", lw=0.8, ls="-.")
    ax_nz.set_ylim(0.0, max(1.15, float(np.nanmax(h["Me"])) * 1.15))
    ax_nz.set_xlim(0.0, t[-1])
    ax_nz.set_xlabel("time  $t$  [s]")
    ax_nz.set_ylabel("$M_e$  [-]", color="#2d6a9f")
    ax_nz.tick_params(axis="y", colors="#2d6a9f")
    ax_n2 = ax_nz.twinx()
    ax_n2.grid(False)
    ax_n2.plot(t, h["pe"] / c.p_a, color="#9b59b6", lw=1.3, ls="--",
               label="$p_e/p_a$")
    ax_n2.plot(t, h["CF"], color="#c0392b", lw=1.3, ls=":", label="$C_F$")
    ax_n2.set_ylabel("$p_e/p_a$  ,  $C_F$  [-]", color="#7b4fa8")
    ax_n2.tick_params(axis="y", colors="#7b4fa8")
    ax_n2.spines["right"].set_visible(True)
    hl = ax_nz.get_legend_handles_labels()
    h2 = ax_n2.get_legend_handles_labels()
    ax_nz.legend(hl[0] + h2[0], hl[1] + h2[1], loc="center right")
    ax_nz.set_title("(f)  Nozzle operating point", loc="left")

    # ---------------- (g) burn rate / burn depth ----------------
    ax_r.plot(t[burn], h["x"][burn] * 1e3, color="#d2691e", lw=1.9,
              label="burn depth $x$")
    ax_r.axhline(web_thickness(c) * 1e3, color="#7a7a7a", lw=0.9, ls="-.",
                 label="web")
    ax_r.set_xlim(0.0, t[-1])
    ax_r.set_xlabel("time  $t$  [s]")
    ax_r.set_ylabel("$x$  [mm]", color="#d2691e")
    ax_r.tick_params(axis="y", colors="#d2691e")
    ax_r.axvline(t_bo, color="#404040", lw=1.0, ls="--")
    ax_r2 = ax_r.twinx()
    ax_r2.grid(False)
    ax_r2.plot(t[burn], h["r"][burn] * 1e3, color="#6b4f2a", lw=1.4, ls="--",
               label="burn rate $r$")
    ax_r2.set_ylabel("$r$  [mm/s]", color="#6b4f2a")
    ax_r2.tick_params(axis="y", colors="#6b4f2a")
    ax_r2.spines["right"].set_visible(True)
    hl = ax_r.get_legend_handles_labels()
    h2 = ax_r2.get_legend_handles_labels()
    ax_r.legend(hl[0] + h2[0], hl[1] + h2[1], loc="center right")
    ax_r.set_title("(g)  Regression of the grain", loc="left")

    # ---------------- title, regime legend, footer ----------------
    fig.suptitle("Solid Rocket Motor - Internal Ballistics\n"
                 f"bore grain ($R_i$={c.R_i0*1e3:.2f} mm, $R_p$={c.R_p*1e3:.2f} mm, "
                 f"$L_p$={c.L_p0*1e3:.1f} mm, {c.n_end:g} burning end face(s))   |   "
                 f"nozzle $R_t$={c.R_t0*1e3:.2f} mm, $A_e/A_t$={c.eps_nozzle:g}   |   "
                 f"$r=${c.a_burn*1e3:.3g} mm/s $\\times(p_0/{c.p_ref/1e6:g}\\,$MPa$)^{{{c.n_burn:g}}}$",
                 fontsize=12.5, fontweight="bold", y=0.985)

    used = [k for k in ("choked", "subsonic", "separated", "no-flow")
            if any(sp[2] == k for sp in spans)]
    handles = [Patch(facecolor=REGIME_FACE[k], edgecolor="#b8c2cc",
                     label=REGIME_NAME[k]) for k in used]
    handles.append(Line2D([0], [0], color="#404040", lw=1.0, ls="--", label="burnout"))
    fig.legend(handles=handles, loc="lower center", ncol=len(handles),
               bbox_to_anchor=(0.5, 0.012), fontsize=8.6, frameon=False)

    note = (f"GRIBS v{PROGRAM_VERSION}  |  thermochemistry: {s['thermochemistry_backend']}   |   "
            f"two-phase: {s['two_phase_model']['mode']}"
            f"{' (complete entrainment, no particle slip)' if s['two_phase_model']['active'] else ''}   |   "
            f"solver: {c.method}, rtol={c.rtol:g}   |   "
            f"unchoked policy: {c.unchoked_policy}   |   "
            f"property policy: {c.property_policy} "
            f"({s['property_range_min_Pa']/1e5:g}-{s['property_range_max_Pa']/1e5:g} bar)   |   "
            f"Cd={c.Cd:g}, eta_F={c.eta_thrust:g}, eta_T0={c.eta_T0:g}   |   "
            f"erosive burning: {'on' if c.ero_alpha > 0 else 'off'}   |   "
            f"throat erosion: {'on' if c.ero_throat_c > 0 else 'off'}")
    fig.text(0.5, 0.037, note, ha="center", va="center", fontsize=7.8,
             color="#4a5560")

    fig.savefig(path, dpi=190)
    plt.close(fig)
