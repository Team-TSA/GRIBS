#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
srm_internal_ballistics.py
================================================================================
HIGH-FIDELITY INTERNAL BALLISTICS OF A SOLID ROCKET MOTOR
  Grain    : cylindrical bore + N burning end face(s) (outer surface inhibited)
  Nozzle   : converging (Ae = At) or converging-diverging (Ae/At > 1)
  Chamber  : unsteady mass/state equation with pressure-dependent gas properties

--------------------------------------------------------------------------------
GOVERNING MODEL
--------------------------------------------------------------------------------
  Equation of state in the chamber        :  p0 * Vg = m_g * Theta(p0),
                                             Theta = R(p0) * T0(p0)
  Differentiating and substituting
      dVg/dt = Ab * r ,  dm_g/dt = mdot_gen - mdot_out
  gives the pressure equation used here (valid for any nozzle regime):

      dp0/dt = [ Theta*(mdot_gen - mdot_out) - p0*Ab*r ] / [ Vg*(1 - p0*Theta'/Theta) ]

  Geometry (exactly consistent: dVg/dx = Ab is verified numerically)
      Ri = Ri0 + x ,  Lp = Lp0 - N_end*x
      Ab = 2*pi*Ri*Lp + N_end*pi*(Rp^2 - Ri^2)
      Vg = Vg0 + Vp0 - pi*(Rp^2 - Ri^2)*Lp

  Burn rate (Saint-Robert + temperature sensitivity + optional erosive burning,
  Lenoir-Robert form, solved implicitly):
      r0 = a*(p0/p_ref)^n * exp(sigma_p*(T_grain - T_ref))
      r  = r0 + alpha*G^0.8/D_h^0.2 * exp(-beta*rho_p*r/G)

  Nozzle (isentropic, discharge coefficient Cd, thrust efficiency eta_F)
      choked      :  mdot = Cd*At*p0*sqrt(g/(R*T0))*(2/(g+1))^((g+1)/(2(g-1)))
                     exit Mach from the area-Mach relation (Ae/At)
      unchoked    :  pe = pa, subsonic exit Mach from the pressure ratio
                     (this branch is *continuous* with the choked branch)
      separated   :  Summerfield criterion pe < f_sep*pa (only for Ae/At > 1)
      thrust      :  F = eta_F*mdot*ve + (pe - pa)*Ae_effective

  Throat erosion :  dRt/dt = C_ero*(p0/p_ref)^m_ero  (Ae is held fixed)

--------------------------------------------------------------------------------
NUMERICS
--------------------------------------------------------------------------------
  State y = [p0, x, Rt, m_out, Impulse, m_gen] integrated with solve_ivp
  (LSODA/BDF/Radau, dense output, tight tolerances).  Quadratures for impulse
  and integrated masses are carried as ODE states, so they inherit the solver
  error control instead of relying on post-hoc trapezoidal sums.
  Root-finding events capture: burnout, choking loss/recovery, ambient pressure.
  A set of numerical self-tests is executed before every production run.

--------------------------------------------------------------------------------
OUTPUT
--------------------------------------------------------------------------------
  results/ballistics.png   single figure, 8 panels, publication-style
  results/time_history.csv full time history
  results/summary.txt      human-readable summary (English)
  results/summary.json     machine-readable summary + the exact input set

  Every physical quantity can be overridden, e.g.
      python srm_internal_ballistics.py --a_burn 1.0e-3 --n_burn 0.35 \
             --R_t0 1.8e-3 --eps_nozzle 4.0 --Cd 0.98 --eta_thrust 0.985
      python srm_internal_ballistics.py --list_params
      python srm_internal_ballistics.py --config my_motor.json
================================================================================
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import dataclass, fields, asdict
from pathlib import Path
from typing import Callable, Optional

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

G0 = 9.80665                      # standard gravity [m/s^2]
_TRAPZ = getattr(np, "trapezoid", np.trapz)


# ==============================================================================
# 1. CONFIGURATION  (every quantity is user-settable)
# ==============================================================================
@dataclass
class Config:
    # ---------------- propellant / combustion ----------------
    rho_p: float = 1769.91        # propellant density [kg/m^3]
    a_burn: float = 0.7e-3        # burn-rate coefficient [m/s] at p = p_ref
    n_burn: float = 0.5           # pressure exponent [-]
    p_ref: float = 1.0e6          # reference pressure of the burn law [Pa]
    sigma_p: float = 0.0          # burn-rate temperature sensitivity [1/K]
    T_grain: float = 294.0        # initial grain temperature [K]
    T_ref: float = 294.0          # reference grain temperature [K]
    ero_alpha: float = 0.0        # Lenoir-Robert erosive coefficient (0 = off)
    ero_beta: float = 53.0        # Lenoir-Robert damping constant [-]

    # ---------------- grain geometry ----------------
    R_i0: float = 3.00e-3         # initial bore radius [m]
    R_p: float = 1.20e-2          # grain outer radius [m]
    L_p0: float = 0.080           # initial grain length [m]
    V_g0: float = 6.79e-6         # initial free (gas) volume [m^3]
    n_end: float = 1.0            # number of burning end faces (0, 1 or 2)

    # ---------------- nozzle ----------------
    R_t0: float = 1.50e-3         # initial throat radius [m]
    eps_nozzle: float = 1.0       # area ratio Ae/At (1.0 = converging nozzle)
    Cd: float = 1.0               # nozzle discharge coefficient [-]
    eta_thrust: float = 1.0       # thrust (momentum) efficiency [-]
    use_separation: bool = True   # model overexpanded flow separation
    sep_ratio: float = 0.4        # Summerfield separation criterion pe/pa [-]
    ero_throat_c: float = 0.0     # throat erosion rate [m/s] at p_ref
    ero_throat_m: float = 1.0     # throat erosion pressure exponent [-]

    # ---------------- environment / initial state ----------------
    p_a: float = 101325.0         # ambient (back) pressure [Pa]
    p0_init: float = 0.2e6        # initial chamber pressure [Pa, absolute]

    # ---------------- igniter (optional gas injection) ----------------
    ign_mdot: float = 0.0         # igniter mass flow [kg/s]
    ign_time: float = 0.0         # igniter duration [s]

    # ---------------- gas property fits  X = A + B*ln(p0/1e5) ----------------
    R_a: float = 385.984          # R  intercept [J/(kg K)]
    R_b: float = -4.00590         # R  slope     [J/(kg K)]
    T_a: float = 2980.31          # T0 intercept [K]
    T_b: float = 115.009          # T0 slope     [K]
    g_a: float = 1.11488          # gamma intercept [-]
    g_b: float = 0.00525420       # gamma slope     [-]
    p_fit_min: float = 1.0e5      # lower validity limit of the fits [Pa]
    p_fit_max: float = 8.0e6      # upper validity limit of the fits [Pa]
    eta_T0: float = 1.0           # combustion/heat-loss efficiency on T0 [-]
    property_policy: str = "clamp"  # 'clamp' | 'extrapolate' (outside fit range)

    # ---------------- nozzle low-pressure policy ----------------
    unchoked_policy: str = "switch"  # 'switch' (subsonic branch) | 'stop'

    # ---------------- solver ----------------
    method: str = "LSODA"         # LSODA | BDF | Radau
    rtol: float = 1.0e-9
    atol_p: float = 1.0e-3        # [Pa]
    atol_x: float = 1.0e-13       # [m]
    t_max: float = 600.0          # integration horizon of the burning phase [s]
    max_step_burn: float = 0.02   # [s]
    blowdown: bool = True         # integrate the post-burnout blowdown
    blowdown_tmax: float = 5.0    # [s]
    max_step_blow: float = 5.0e-4  # [s]


PARAM_DOC = {
    "rho_p": "propellant density [kg/m^3]",
    "a_burn": "burn-rate coefficient [m/s at p_ref]",
    "n_burn": "burn-rate pressure exponent [-]",
    "p_ref": "reference pressure of the burn law [Pa]",
    "sigma_p": "burn-rate temperature sensitivity [1/K]",
    "T_grain": "initial grain temperature [K]",
    "T_ref": "reference grain temperature [K]",
    "ero_alpha": "Lenoir-Robert erosive-burning coefficient (0 disables it)",
    "ero_beta": "Lenoir-Robert damping constant [-]",
    "R_i0": "initial bore radius [m]",
    "R_p": "grain outer radius [m]",
    "L_p0": "initial grain length [m]",
    "V_g0": "initial free gas volume [m^3]",
    "n_end": "number of burning end faces (0, 1 or 2)",
    "R_t0": "initial throat radius [m]",
    "eps_nozzle": "nozzle area ratio Ae/At [-]",
    "Cd": "nozzle discharge coefficient [-]",
    "eta_thrust": "thrust (momentum) efficiency [-]",
    "use_separation": "model overexpanded flow separation (true/false)",
    "sep_ratio": "separation criterion pe/pa [-]",
    "ero_throat_c": "throat erosion rate at p_ref [m/s]",
    "ero_throat_m": "throat erosion pressure exponent [-]",
    "p_a": "ambient back pressure [Pa]",
    "p0_init": "initial chamber pressure [Pa abs]",
    "ign_mdot": "igniter mass flow [kg/s]",
    "ign_time": "igniter duration [s]",
    "R_a": "gas-constant fit intercept [J/(kg K)]",
    "R_b": "gas-constant fit slope [J/(kg K)]",
    "T_a": "flame-temperature fit intercept [K]",
    "T_b": "flame-temperature fit slope [K]",
    "g_a": "specific-heat-ratio fit intercept [-]",
    "g_b": "specific-heat-ratio fit slope [-]",
    "p_fit_min": "lower validity limit of the property fits [Pa]",
    "p_fit_max": "upper validity limit of the property fits [Pa]",
    "eta_T0": "combustion efficiency applied to T0 [-]",
    "property_policy": "'clamp' or 'extrapolate' outside the fit range",
    "unchoked_policy": "'switch' to subsonic branch, or 'stop'",
    "method": "ODE method: LSODA | BDF | Radau",
    "rtol": "relative tolerance",
    "atol_p": "absolute tolerance on pressure [Pa]",
    "atol_x": "absolute tolerance on burn depth [m]",
    "t_max": "integration horizon of the burning phase [s]",
    "max_step_burn": "maximum solver step, burning phase [s]",
    "blowdown": "integrate post-burnout blowdown (true/false)",
    "blowdown_tmax": "blowdown horizon [s]",
    "max_step_blow": "maximum solver step, blowdown [s]",
}


# ==============================================================================
# 2. DERIVED GEOMETRY HELPERS
# ==============================================================================
def web_thickness(c: Config) -> float:
    """Burn depth at which the grain is consumed (first geometric limit)."""
    radial = c.R_p - c.R_i0
    axial = c.L_p0 / c.n_end if c.n_end > 0.0 else math.inf
    return min(radial, axial)


def initial_propellant_volume(c: Config) -> float:
    return math.pi * (c.R_p ** 2 - c.R_i0 ** 2) * c.L_p0


def propellant_mass(c: Config) -> float:
    return c.rho_p * initial_propellant_volume(c)


def geometry(x: float, c: Config):
    """Return (Ri, Lp, Ab, Vg) for burn depth x [m] (clamped to the web)."""
    xw = web_thickness(c)
    xe = min(max(x, 0.0), xw)
    Ri = min(c.R_i0 + xe, c.R_p)
    Lp = max(c.L_p0 - c.n_end * xe, 0.0)
    dA = max(c.R_p ** 2 - Ri ** 2, 0.0)
    Ab = 2.0 * math.pi * Ri * Lp + c.n_end * math.pi * dA
    Vprop = math.pi * dA * Lp
    Vg = c.V_g0 + initial_propellant_volume(c) - Vprop
    return Ri, Lp, Ab, Vg


# ==============================================================================
# 3. GAS PROPERTIES
# ==============================================================================
def gas_props(p0: float, c: Config):
    """Return (R [J/(kg K)], T0 [K], gamma [-]) for absolute pressure p0 [Pa]."""
    pe = p0
    if c.property_policy == "clamp":
        pe = min(max(p0, c.p_fit_min), c.p_fit_max)
    pe = max(pe, 1.0)
    ell = math.log(pe / 1.0e5)
    R = c.R_a + c.R_b * ell
    T0 = (c.T_a + c.T_b * ell) * c.eta_T0
    g = c.g_a + c.g_b * ell
    if R <= 0.0 or T0 <= 0.0 or g <= 1.0:
        raise ValueError("Non-physical thermodynamic properties from the fit.")
    return R, T0, g


def theta_and_deriv(p0: float, c: Config):
    """Return Theta = R*T0 [J/kg] and dTheta/dp0 [J/(kg Pa)]."""
    R, T0, _ = gas_props(p0, c)
    theta = R * T0
    clamped = (c.property_policy == "clamp"
               and (p0 < c.p_fit_min or p0 > c.p_fit_max))
    dtheta = 0.0 if clamped else (c.R_b * T0 + R * c.eta_T0 * c.T_b) / p0
    return theta, dtheta


def is_extrapolated(p0: float, c: Config) -> bool:
    return bool(p0 < c.p_fit_min or p0 > c.p_fit_max)


# ==============================================================================
# 4. BURN RATE (Saint-Robert + temperature sensitivity + erosive burning)
# ==============================================================================
def base_burn_rate(p0: float, c: Config) -> float:
    return (c.a_burn * (max(p0, 0.0) / c.p_ref) ** c.n_burn
            * math.exp(c.sigma_p * (c.T_grain - c.T_ref)))


def burn_rate(p0: float, c: Config, Ab: float, Ri: float) -> float:
    """Total regression rate including the optional erosive contribution."""
    r0 = base_burn_rate(p0, c)
    if c.ero_alpha <= 0.0 or Ri <= 0.0 or Ab <= 0.0:
        return r0
    A_port = math.pi * Ri ** 2
    D_h = 2.0 * Ri
    r = r0
    for _ in range(60):                     # damped fixed-point iteration
        G = c.rho_p * r * Ab / A_port       # port mass flux at the aft end
        if G <= 1.0e-12:
            return r0
        r_new = r0 + c.ero_alpha * G ** 0.8 / D_h ** 0.2 \
            * math.exp(-c.ero_beta * c.rho_p * r / G)
        if abs(r_new - r) <= 1.0e-13 * max(r, 1.0e-12):
            return r_new
        r = 0.5 * (r + r_new)
    return r


# ==============================================================================
# 5. NOZZLE MODEL
# ==============================================================================
def area_mach(M: float, g: float) -> float:
    """A/At as a function of Mach number (isentropic)."""
    return (1.0 / M) * ((2.0 / (g + 1.0))
                        * (1.0 + 0.5 * (g - 1.0) * M * M)) ** ((g + 1.0) / (2.0 * (g - 1.0)))


def mach_from_area(eps: float, g: float, supersonic: bool) -> float:
    """Invert the area-Mach relation."""
    if eps <= 1.0 + 1.0e-14:
        return 1.0
    f = lambda M: area_mach(M, g) - eps
    if supersonic:
        return brentq(f, 1.0 + 1.0e-12, 200.0, xtol=1e-14, rtol=8.9e-16, maxiter=200)
    return brentq(f, 1.0e-7, 1.0, xtol=1e-16, rtol=8.9e-16, maxiter=200)


def critical_ratio(g: float) -> float:
    """pe/p0 at the throat when choked."""
    return (2.0 / (g + 1.0)) ** (g / (g - 1.0))


def choke_threshold_pressure(p0: float, At: float, Ae: float, c: Config) -> float:
    """Chamber pressure above which the nozzle is choked (evaluated at gamma(p0))."""
    _, _, g = gas_props(p0, c)
    eps = Ae / At
    if eps <= 1.0 + 1.0e-12:
        return c.p_a / critical_ratio(g)
    Me_sub = mach_from_area(eps, g, supersonic=False)
    return c.p_a * (1.0 + 0.5 * (g - 1.0) * Me_sub ** 2) ** (g / (g - 1.0))


def choke_margin(p0: float, At: float, Ae: float, c: Config) -> float:
    """> 0 when choked, < 0 when subsonic."""
    return p0 - choke_threshold_pressure(p0, At, Ae, c)


def choke_limit_pressure(c: Config, At: float, Ae: float) -> Optional[float]:
    """Solve p0 = p_threshold(p0) (the fixed point marking the choking limit)."""
    f = lambda p: choke_margin(p, At, Ae, c)
    lo, hi = c.p_a * (1.0 + 1e-9), c.p_a * 1.0e4
    try:
        if f(lo) * f(hi) > 0.0:
            return None
        return brentq(f, lo, hi, xtol=1e-8, rtol=1e-14)
    except (ValueError, RuntimeError):
        return None


def nozzle_state(p0: float, At: float, Ae: float, c: Config) -> dict:
    """Mass flow, thrust, exit conditions and flow regime."""
    if p0 <= c.p_a:
        return dict(mdot=0.0, F=0.0, Me=0.0, pe=p0, ve=0.0,
                    regime="no-flow", A_eff=Ae)

    R, T0, g = gas_props(p0, c)
    eps = Ae / At
    sq = math.sqrt(g / (R * T0))
    p_thr = choke_threshold_pressure(p0, At, Ae, c)

    if p0 >= p_thr:                                   # ---- choked ----
        mdot = c.Cd * At * p0 * sq \
            * (2.0 / (g + 1.0)) ** ((g + 1.0) / (2.0 * (g - 1.0)))
        Me = 1.0 if eps <= 1.0 + 1.0e-12 else mach_from_area(eps, g, True)
        pe = p0 * (1.0 + 0.5 * (g - 1.0) * Me * Me) ** (-g / (g - 1.0))
        Te = T0 / (1.0 + 0.5 * (g - 1.0) * Me * Me)
        ve = Me * math.sqrt(g * R * Te)
        A_eff, pe_eff, ve_eff, regime = Ae, pe, ve, "choked"

        if c.use_separation and eps > 1.0 and pe < c.sep_ratio * c.p_a:
            p_sep = c.sep_ratio * c.p_a
            arg = (p0 / p_sep) ** ((g - 1.0) / g) - 1.0
            if arg > 0.0:
                M_sep = math.sqrt(2.0 / (g - 1.0) * arg)
                eps_sep = area_mach(M_sep, g)
                if 1.0 < eps_sep < eps:
                    Te_s = T0 / (1.0 + 0.5 * (g - 1.0) * M_sep ** 2)
                    A_eff = eps_sep * At
                    pe_eff = p_sep
                    ve_eff = M_sep * math.sqrt(g * R * Te_s)
                    Me = M_sep
                    regime = "separated"
        F = c.eta_thrust * mdot * ve_eff + (pe_eff - c.p_a) * A_eff
        return dict(mdot=mdot, F=F, Me=Me, pe=pe_eff, ve=ve_eff,
                    regime=regime, A_eff=A_eff)

    # ---- subsonic (exit pressure equals ambient); continuous with the above ----
    Me = math.sqrt(2.0 / (g - 1.0) * ((p0 / c.p_a) ** ((g - 1.0) / g) - 1.0))
    mdot = c.Cd * Ae * p0 * sq * Me \
        * (1.0 + 0.5 * (g - 1.0) * Me * Me) ** (-(g + 1.0) / (2.0 * (g - 1.0)))
    Te = T0 / (1.0 + 0.5 * (g - 1.0) * Me * Me)
    ve = Me * math.sqrt(g * R * Te)
    F = c.eta_thrust * mdot * ve
    return dict(mdot=mdot, F=F, Me=Me, pe=c.p_a, ve=ve,
                regime="subsonic", A_eff=Ae)


# ==============================================================================
# 6. QUASI-STEADY EQUILIBRIUM PRESSURE (diagnostic only)
# ==============================================================================
def equilibrium_pressure(Ab: float, Ri: float, At: float, Ae: float,
                         c: Config) -> Optional[float]:
    """Lowest pressure satisfying mdot_gen(p) = mdot_out(p); None if not found."""
    if not np.isfinite(Ab) or Ab <= 0.0:
        return None

    def bal(p):
        try:
            v = c.rho_p * burn_rate(p, c, Ab, Ri) * Ab - nozzle_state(p, At, Ae, c)["mdot"]
            return float(v) if np.isfinite(v) else np.nan
        except (ValueError, OverflowError, FloatingPointError):
            return np.nan

    grid = np.geomspace(c.p_a * (1.0 + 1e-9), 1.0e11, 140)
    vals = np.array([bal(p) for p in grid])
    for p1, p2, f1, f2 in zip(grid[:-1], grid[1:], vals[:-1], vals[1:]):
        if not (np.isfinite(f1) and np.isfinite(f2)):
            continue
        if f1 == 0.0:
            return float(p1)
        if f1 * f2 < 0.0:
            return float(brentq(bal, p1, p2, xtol=1.0, rtol=1e-12))
    return None


# ==============================================================================
# 7. ODE SYSTEM
#    y = [p0, x, Rt, m_out, Impulse, m_gen]
# ==============================================================================
IP, IX, IRT, IMO, IIMP, IMG = range(6)


def igniter_mdot(t: float, c: Config) -> float:
    return c.ign_mdot if (c.ign_time > 0.0 and t <= c.ign_time) else 0.0


def make_rhs(c: Config, Ae: float, burning: bool) -> Callable:
    xw = web_thickness(c)
    p_floor = 500.0

    def rhs(t, y):
        p0 = max(float(y[IP]), p_floor)
        x = min(max(float(y[IX]), 0.0), xw)
        Rt = max(float(y[IRT]), 1.0e-9)
        At = math.pi * Rt * Rt

        Ri, Lp, Ab, Vg = geometry(x, c)
        if burning and Ab > 0.0:
            r = burn_rate(p0, c, Ab, Ri)
        else:
            r, Ab = 0.0, 0.0
        m_gen = c.rho_p * r * Ab + igniter_mdot(t, c)

        nz = nozzle_state(p0, At, Ae, c)
        theta, dtheta = theta_and_deriv(p0, c)
        denom = Vg * (1.0 - p0 * dtheta / theta)
        if denom <= 0.0:
            raise ValueError("Non-physical factor (1 - p0*Theta'/Theta) <= 0.")

        dp = (theta * (m_gen - nz["mdot"]) - p0 * Ab * r) / denom
        dRt = c.ero_throat_c * (p0 / c.p_ref) ** c.ero_throat_m
        return np.array([dp, r, dRt, nz["mdot"], nz["F"], m_gen])

    return rhs


def _event(fn, terminal=True, direction=0.0):
    fn.terminal = terminal
    fn.direction = direction
    return fn


def run_model(c: Config) -> dict:
    """Integrate the burning phase and (optionally) the blowdown phase."""
    validate(c)
    At0 = math.pi * c.R_t0 ** 2
    Ae = c.eps_nozzle * At0                      # geometric exit area, fixed
    xw = web_thickness(c)

    init_regime = nozzle_state(c.p0_init, At0, Ae, c)["regime"]
    if c.unchoked_policy == "stop" and init_regime != "choked":
        raise ValueError("Initial state is not choked: policy 'stop' has "
                         "nothing to integrate.")

    theta0, _ = theta_and_deriv(c.p0_init, c)
    y0 = np.array([c.p0_init, 0.0, c.R_t0, 0.0, 0.0, 0.0])
    m_gas0 = c.p0_init * c.V_g0 / theta0

    atol = np.array([c.atol_p, c.atol_x, 1e-14, 1e-12, 1e-9, 1e-12])

    ev_burn = _event(lambda t, y: y[IX] - xw, True, +1.0)
    ev_unch = _event(lambda t, y: choke_margin(max(y[IP], 1.0),
                                               math.pi * max(y[IRT], 1e-9) ** 2,
                                               Ae, c),
                     c.unchoked_policy == "stop", -1.0)
    ev_rech = _event(lambda t, y: choke_margin(max(y[IP], 1.0),
                                               math.pi * max(y[IRT], 1e-9) ** 2,
                                               Ae, c), False, +1.0)
    ev_amb = _event(lambda t, y: y[IP] - c.p_a * 1.0005, True, -1.0)
    names = ["burnout", "unchoked", "rechoked", "ambient"]

    sol1 = solve_ivp(make_rhs(c, Ae, True), (0.0, c.t_max), y0,
                     method=c.method, events=[ev_burn, ev_unch, ev_rech, ev_amb],
                     dense_output=True, rtol=c.rtol, atol=atol,
                     max_step=c.max_step_burn, first_step=1.0e-7)
    if not sol1.success:
        raise RuntimeError(f"Burning-phase integration failed: {sol1.message}")

    stop1 = "time-limit"
    for nm, te in zip(names, sol1.t_events):
        if nm in ("burnout", "ambient") and len(te):
            stop1 = nm
            break
        if nm == "unchoked" and c.unchoked_policy == "stop" and len(te):
            stop1 = "choking-loss"
            break

    sol2, stop2 = None, "not-run"
    if c.blowdown and stop1 == "burnout":
        y_bo = sol1.y[:, -1].copy()
        y_bo[IX] = xw
        ev_amb2 = _event(lambda t, y: y[IP] - c.p_a * 1.0005, True, -1.0)
        ev_unch2 = _event(lambda t, y: choke_margin(max(y[IP], 1.0),
                                                    math.pi * max(y[IRT], 1e-9) ** 2,
                                                    Ae, c),
                          c.unchoked_policy == "stop", -1.0)
        sol2 = solve_ivp(make_rhs(c, Ae, False),
                         (sol1.t[-1], sol1.t[-1] + c.blowdown_tmax), y_bo,
                         method=c.method, events=[ev_amb2, ev_unch2],
                         dense_output=True, rtol=c.rtol, atol=atol,
                         max_step=c.max_step_blow)
        if not sol2.success:
            raise RuntimeError(f"Blowdown integration failed: {sol2.message}")
        stop2 = "ambient" if len(sol2.t_events[0]) else "time-limit"
        if c.unchoked_policy == "stop" and len(sol2.t_events[1]):
            stop2 = "choking-loss"

    return dict(sol_burn=sol1, sol_blow=sol2, stop_burn=stop1, stop_blow=stop2,
                Ae=Ae, At0=At0, m_gas0=m_gas0, init_regime=init_regime,
                t_events={nm: [float(v) for v in te]
                          for nm, te in zip(names, sol1.t_events)})


def validate(c: Config) -> None:
    if c.R_i0 >= c.R_p:
        raise ValueError("R_i0 must be smaller than R_p.")
    for name in ("rho_p", "a_burn", "p_ref", "R_p", "L_p0", "V_g0", "R_t0",
                 "p_a", "p0_init", "T_grain", "T_ref"):
        if getattr(c, name) <= 0.0:
            raise ValueError(f"{name} must be positive.")
    if c.eps_nozzle < 1.0:
        raise ValueError("eps_nozzle (Ae/At) must be >= 1.")
    if c.n_end < 0.0 or c.n_end > 2.0:
        raise ValueError("n_end must lie between 0 and 2.")
    if c.property_policy not in ("clamp", "extrapolate"):
        raise ValueError("property_policy must be 'clamp' or 'extrapolate'.")
    if c.unchoked_policy not in ("switch", "stop"):
        raise ValueError("unchoked_policy must be 'switch' or 'stop'.")


# ==============================================================================
# 8. POST-PROCESSING
# ==============================================================================
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
    Ae = res["Ae"]
    xw = web_thickness(c)
    segs = [(res["sol_burn"], "burn", 4000, True)]
    if res["sol_blow"] is not None:
        segs.append((res["sol_blow"], "blowdown", 1600, True))

    rec = {k: [] for k in
           ("t", "phase", "p0", "x", "Rt", "At", "Kn", "Ri", "Lp", "Ab", "Vg",
            "r", "mdot_gen", "mdot_out", "F", "Me", "pe", "ve", "regime",
            "R", "T0", "gamma", "m_out", "impulse", "m_gen", "m_gas_eos",
            "m_gas_bal", "extrap", "cstar", "CF")}

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
            R, T0, g = gas_props(p0, c)
            nz = nozzle_state(p0, At, Ae, c)
            if phase == "burn":
                r = burn_rate(p0, c, Ab, Ri)
                mg = c.rho_p * r * Ab + igniter_mdot(t, c)
            else:
                r, Ab, mg = 0.0, 0.0, 0.0
            theta = R * T0
            m_eos = p0 * Vg / theta
            m_bal = res["m_gas0"] + float(ys[IMG, k]) - float(ys[IMO, k])

            rec["t"].append(float(t));          rec["phase"].append(phase)
            rec["p0"].append(p0);               rec["x"].append(x)
            rec["Rt"].append(Rt);               rec["At"].append(At)
            rec["Kn"].append(Ab / At);          rec["Ri"].append(Ri)
            rec["Lp"].append(Lp);               rec["Ab"].append(Ab)
            rec["Vg"].append(Vg);               rec["r"].append(r)
            rec["mdot_gen"].append(mg);         rec["mdot_out"].append(nz["mdot"])
            rec["F"].append(nz["F"]);           rec["Me"].append(nz["Me"])
            rec["pe"].append(nz["pe"]);         rec["ve"].append(nz["ve"])
            rec["regime"].append(nz["regime"]); rec["R"].append(R)
            rec["T0"].append(T0);               rec["gamma"].append(g)
            rec["m_out"].append(float(ys[IMO, k]))
            rec["impulse"].append(float(ys[IIMP, k]))
            rec["m_gen"].append(float(ys[IMG, k]))
            rec["m_gas_eos"].append(m_eos);     rec["m_gas_bal"].append(m_bal)
            rec["extrap"].append(is_extrapolated(p0, c))
            rec["cstar"].append(p0 * At / nz["mdot"] if nz["mdot"] > 0 else np.nan)
            rec["CF"].append(nz["F"] / (p0 * At) if p0 * At > 0 else np.nan)

    out = {k: (np.asarray(v) if k not in ("phase", "regime") else np.asarray(v, dtype=object))
           for k, v in rec.items()}

    # quasi-steady equilibrium pressure (diagnostic, burning phase only)
    t = out["t"]
    p_eq = np.full_like(t, np.nan)
    burn_idx = np.nonzero(out["phase"] == "burn")[0]
    if burn_idx.size:
        pick = burn_idx[np.unique(np.linspace(0, burn_idx.size - 1,
                                              min(200, burn_idx.size)).astype(int))]
        vals = np.array([equilibrium_pressure(out["Ab"][i], out["Ri"][i],
                                              out["At"][i], Ae, c) or np.nan
                         for i in pick])
        ok = np.isfinite(vals)
        if ok.sum() >= 2:
            p_eq[burn_idx] = np.interp(t[burn_idx], t[pick][ok], vals[ok],
                                       left=np.nan, right=np.nan)
    out["p_eq"] = p_eq
    return out


def summarize(res: dict, h: dict, c: Config) -> dict:
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

    pAt_int = float(_TRAPZ(h["p0"][burn] * h["At"][burn], t[burn]))
    cstar_eff = pAt_int / m_out_tot if m_out_tot > 0 else float("nan")
    CF_eff = I_burn / pAt_int if pAt_int > 0 else float("nan")
    p_mean = float(pAt_int / (At0 * (t_bo - t[0]))) if t_bo > t[0] else float("nan")

    denom = np.maximum(np.abs(h["m_gas_eos"]), 1e-15)
    mass_err = float(np.max(np.abs(h["m_gas_bal"] - h["m_gas_eos"]) / denom))

    p_star = choke_limit_pressure(c, At0, Ae)
    mg0 = c.rho_p * base_burn_rate(c.p0_init, c) * geometry(0.0, c)[2]
    mo0 = nozzle_state(c.p0_init, At0, Ae, c)["mdot"]

    s = dict(
        web_mm=web_thickness(c) * 1e3,
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
        Isp_s=I_tot / (m_prop * G0) if m_prop > 0 else float("nan"),
        cstar_eff_m_s=cstar_eff,
        CF_eff=CF_eff,
        Kn_initial=float(h["Kn"][0]), Kn_max=float(np.max(h["Kn"][burn])),
        Kn_final=float(h["Kn"][burn][-1]),
        throat_radius_final_mm=float(h["Rt"][-1]) * 1e3,
        propellant_mass_balance_error=abs(m_gen_tot / m_prop - 1.0) if m_prop > 0 else float("nan"),
        gas_mass_consistency_error=mass_err,
        property_extrapolation=bool(np.any(h["extrap"])),
        regimes_visited=sorted(set(h["regime"].tolist())),
        event_times=res["t_events"],
    )
    return s


# ==============================================================================
# 9. SELF-TESTS
# ==============================================================================
CEA_TABLE = {1: (385.980, 2981.28, 1.1151), 2: (383.223, 3059.37, 1.1185),
             5: (379.549, 3164.47, 1.1232), 10: (376.746, 3244.82, 1.1268),
             20: (373.967, 3325.29, 1.1305), 30: (372.342, 3372.15, 1.1327),
             50: (370.302, 3430.66, 1.1355), 80: (368.464, 3483.70, 1.1381)}


def self_tests(c: Config, verbose=True) -> bool:
    out, ok_all = [], True

    def chk(name, ok, detail=""):
        nonlocal ok_all
        ok_all &= bool(ok)
        out.append(f"  [{'PASS' if ok else 'FAIL'}] {name:<46s} {detail}")

    # T1 property fits versus the tabulated CEA data (only if defaults are used)
    ref = Config()
    if all(getattr(c, k) == getattr(ref, k)
           for k in ("R_a", "R_b", "T_a", "T_b", "g_a", "g_b", "eta_T0")):
        e = 0.0
        for pbar, (Rv, Tv, gv) in CEA_TABLE.items():
            R, T0, g = gas_props(pbar * 1e5, c)
            e = max(e, abs(R / Rv - 1), abs(T0 / Tv - 1), abs(g / gv - 1))
        chk("property fits vs CEA table", e < 5e-3, f"max rel. err = {e:.2e}")
    else:
        chk("property fits vs CEA table", True, "skipped (custom coefficients)")

    # T2 analytic dTheta/dp against a centred difference
    e = 0.0
    for p in np.geomspace(max(c.p_fit_min * 1.05, 1.1e5), c.p_fit_max * 0.95, 25):
        hstep = p * 1e-5
        fd = (theta_and_deriv(p + hstep, c)[0] - theta_and_deriv(p - hstep, c)[0]) / (2 * hstep)
        e = max(e, abs(fd / theta_and_deriv(p, c)[1] - 1.0))
    chk("analytic vs numerical dTheta/dp", e < 1e-7, f"max rel. err = {e:.2e}")

    # T3 geometric consistency dVg/dx = Ab
    xw = web_thickness(c)
    e = 0.0
    for x in np.linspace(0.02 * xw, 0.98 * xw, 25):
        hstep = xw * 1e-6
        fd = (geometry(x + hstep, c)[3] - geometry(x - hstep, c)[3]) / (2 * hstep)
        e = max(e, abs(fd / geometry(x, c)[2] - 1.0))
    chk("geometry consistency dVg/dx = Ab", e < 1e-7, f"max rel. err = {e:.2e}")

    # T4 area-Mach inversion round trip
    e = 0.0
    for eps in (1.5, 2.5, 6.0, 20.0):
        for g in (1.10, 1.15, 1.25):
            e = max(e, abs(area_mach(mach_from_area(eps, g, True), g) / eps - 1.0),
                    abs(area_mach(mach_from_area(eps, g, False), g) / eps - 1.0))
    chk("area-Mach inversion round trip", e < 1e-10, f"max rel. err = {e:.2e}")

    # T5 continuity of mdot and F across the choking boundary
    # Test only the nozzle geometry specified in the current configuration.
    worst_mdot = 0.0
    worst_thrust = 0.0

    eps = c.eps_nozzle

    cc = Config(**{
        **asdict(c),
        "eps_nozzle": eps,
        "use_separation": False
    })

    At0 = math.pi * cc.R_t0 ** 2
    Ae = eps * At0

    pstar = choke_limit_pressure(cc, At0, Ae)

    if pstar is None:
        chk("choked/subsonic continuity (mdot, F)", False,
            "choking-limit pressure could not be determined")
    else:
        lo = nozzle_state(pstar * (1.0 - 1e-9), At0, Ae, cc)
        hi = nozzle_state(pstar * (1.0 + 1e-9), At0, Ae, cc)
        worst_mdot = abs(lo["mdot"] / hi["mdot"] - 1.0)
        worst_thrust = abs((lo["F"] + 1e-12) / (hi["F"] + 1e-12) - 1.0)
        worst = max(worst_mdot, worst_thrust)
        chk("choked/subsonic continuity (mdot, F)", worst < 1e-6,
            f"Ae/At = {eps:.6g}, mdot jump = {worst_mdot:.2e}, "
            f"F jump = {worst_thrust:.2e}")

        # T6 closed-form check of the converging-nozzle choked thrust
        cc = Config(**{**asdict(c), "eps_nozzle": 1.0, "Cd": 1.0,
                       "eta_thrust": 1.0, "use_separation": False})
        At0 = math.pi * cc.R_t0 ** 2
        e = 0.0
        for p in np.geomspace(cc.p_a / 0.5, 50e6, 20):
            nz = nozzle_state(p, At0, At0, cc)
            _, _, g = gas_props(p, cc)
            F_ref = At0 * ((g + 1.0) * critical_ratio(g) * p - cc.p_a)
            e = max(e, abs(nz["F"] / F_ref - 1.0))
        chk("closed-form choked thrust identity", e < 1e-12,
            f"max rel. err = {e:.2e}")

        # T7 finite and physical right-hand side at t = 0
        try:
            dy = make_rhs(c, c.eps_nozzle * math.pi * c.R_t0 ** 2, True)(
                0.0, np.array([c.p0_init, 0.0, c.R_t0, 0.0, 0.0, 0.0]))
            chk("finite initial derivatives", np.all(np.isfinite(dy)),
                f"dp/dt = {dy[0]:+.3e} Pa/s, r = {dy[1]:.3e} m/s")
        except Exception as exc:                              # pragma: no cover
            chk("finite initial derivatives", False, str(exc))

    if verbose:
        print("Self-tests")
        print("\n".join(out))
        print(f"  -> {'all checks passed' if ok_all else 'FAILURES DETECTED'}\n")
    return ok_all


# ==============================================================================
# 10. PLOTTING (one single, carefully designed figure)
# ==============================================================================
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


def make_figure(res, h, s, c: Config, path: Path) -> None:
    _style()
    t = h["t"]
    burn = h["phase"] == "burn"
    t_bo = s["burn_time_s"]
    spans = _regime_spans(t, h["regime"])
    has_blow = res["sol_blow"] is not None

    fig = plt.figure(figsize=(16.6, 13.0))
    gs = GridSpec(3, 3, figure=fig, height_ratios=[1.55, 1.0, 1.0],
                  hspace=0.40, wspace=0.28,
                  left=0.058, right=0.978, top=0.915, bottom=0.105)

    ax_p = fig.add_subplot(gs[0, 0:2])
    ax_key = fig.add_subplot(gs[0, 2])
    ax_F = fig.add_subplot(gs[1, 0])
    ax_md = fig.add_subplot(gs[1, 1])
    ax_geo = fig.add_subplot(gs[1, 2])
    ax_zm = fig.add_subplot(gs[2, 0])
    ax_nz = fig.add_subplot(gs[2, 1])
    ax_r = fig.add_subplot(gs[2, 2])

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
                  solid_capstyle="butt", label="outside property-fit range")
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
    ax_p.legend(loc="lower right", ncol=1)

    # ---------------- key results panel ----------------
    ax_key.axis("off")
    ev = s["event_times"]
    lines = [
        "KEY RESULTS",
        "-" * 42,
        f"{'web thickness':<24}{s['web_mm']:>11.3f} mm",
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
        "VERIFICATION",
        f"{'propellant mass balance':<24}{s['propellant_mass_balance_error']:>11.2e}",
        f"{'gas mass consistency':<24}{s['gas_mass_consistency_error']:>11.2e}",
        f"{'property extrapolation':<24}{str(s['property_extrapolation']):>11s}",
    ]
    if ev.get("unchoked"):
        lines.append(f"{'first unchoking at':<24}{ev['unchoked'][0]*1e3:>11.3f} ms")
    if ev.get("rechoked"):
        lines.append(f"{'re-choking at':<24}{ev['rechoked'][0]*1e3:>11.3f} ms")
    ax_key.text(0.0, 1.0, "\n".join(lines), transform=ax_key.transAxes,
                va="top", ha="left", family="monospace", fontsize=8.1,
                linespacing=1.28,
                bbox=dict(boxstyle="round,pad=0.55", fc="#f4f7fa", ec="#a9bbcc", lw=1.0))

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
    pos = h["mdot_gen"] > 0
    ax_md.plot(t[pos], h["mdot_gen"][pos] * 1e3, color=C_GEN, lw=1.8,
               label=r"generated  $\dot{m}_{gen}=\rho_p A_b r$")
    ax_md.plot(t, np.maximum(h["mdot_out"], 1e-12) * 1e3, color=C_OUT, lw=1.6,
               ls="--", label=r"nozzle  $\dot{m}_{out}$")
    acc = h["mdot_gen"] > h["mdot_out"]
    ax_md.fill_between(t, np.maximum(h["mdot_out"], 1e-12) * 1e3,
                       np.maximum(h["mdot_gen"], 1e-12) * 1e3, where=acc,
                       color=C_ACC, alpha=0.25, lw=0, label="chamber filling")
    ax_md.axvline(t_bo, color="#404040", lw=1.0, ls="--")
    ax_md.set_yscale("log")
    ax_md.set_xlim(0.0, t[-1])
    ax_md.set_xlabel("time  $t$  [s]")
    ax_md.set_ylabel(r"mass flow  [g/s]")
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
    R0, T00, g00 = gas_props(c.p0_init, c)
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

    note = (f"solver: {c.method}, rtol={c.rtol:g}   |   "
            f"unchoked policy: {c.unchoked_policy}   |   "
            f"property policy: {c.property_policy} "
            f"({c.p_fit_min/1e5:g}-{c.p_fit_max/1e5:g} bar)   |   "
            f"Cd={c.Cd:g}, eta_F={c.eta_thrust:g}, eta_T0={c.eta_T0:g}   |   "
            f"erosive burning: {'on' if c.ero_alpha > 0 else 'off'}   |   "
            f"throat erosion: {'on' if c.ero_throat_c > 0 else 'off'}")
    fig.text(0.5, 0.037, note, ha="center", va="center", fontsize=7.8,
             color="#4a5560")

    fig.savefig(path, dpi=190)
    plt.close(fig)


# ==============================================================================
# 11. FILE OUTPUT
# ==============================================================================
CSV_COLUMNS = ["t", "phase", "p0", "x", "Ri", "Lp", "Ab", "Vg", "Kn", "Rt", "At",
               "r", "mdot_gen", "mdot_out", "m_gen", "m_out", "m_gas_eos",
               "m_gas_bal", "regime", "Me", "pe", "ve", "F", "impulse",
               "cstar", "CF", "R", "T0", "gamma", "extrap"]
CSV_UNITS = {"t": "s", "p0": "Pa", "x": "m", "Ri": "m", "Lp": "m", "Ab": "m2",
             "Vg": "m3", "Kn": "-", "Rt": "m", "At": "m2", "r": "m/s",
             "mdot_gen": "kg/s", "mdot_out": "kg/s", "m_gen": "kg", "m_out": "kg",
             "m_gas_eos": "kg", "m_gas_bal": "kg", "Me": "-", "pe": "Pa",
             "ve": "m/s", "F": "N", "impulse": "N.s", "cstar": "m/s", "CF": "-",
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


def summary_text(s: dict, c: Config) -> str:
    L = []
    A = L.append
    A("=" * 78)
    A("SOLID ROCKET MOTOR - INTERNAL BALLISTICS SUMMARY")
    A("=" * 78)
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
    A("[ verification ]")
    A(f"  propellant mass balance error    : {s['propellant_mass_balance_error']:.3e}")
    A(f"  gas-mass (EOS vs balance) error  : {s['gas_mass_consistency_error']:.3e}")
    A(f"  property-fit extrapolation used  : {s['property_extrapolation']}"
      f"   (policy = {c.property_policy})")
    if s["property_extrapolation"] and c.property_policy == "extrapolate":
        A("  WARNING: results rely on extrapolated CEA property fits.")
    A("=" * 78)
    return "\n".join(L) + "\n"


# ==============================================================================
# 12. CLI
# ==============================================================================
def _str2bool(v: str) -> bool:
    if str(v).lower() in ("1", "true", "yes", "on", "y"):
        return True
    if str(v).lower() in ("0", "false", "no", "off", "n"):
        return False
    raise argparse.ArgumentTypeError("expected a boolean value")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="High-fidelity internal ballistics of a solid rocket motor.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outdir", type=Path,
                   default=Path(__file__).resolve().parent / "results",
                   help="output directory")
    p.add_argument("--config", type=Path, default=None,
                   help="JSON file with parameter overrides")
    p.add_argument("--selftest", action="store_true",
                   help="run the numerical self-tests only")
    p.add_argument("--list_params", action="store_true",
                   help="list every adjustable parameter and exit")
    grp = p.add_argument_group("physical / numerical parameters "
                               "(any value can be overridden)")
    for f in fields(Config):
        tp = f.type if isinstance(f.type, str) else getattr(f.type, "__name__", "str")
        kind = {"float": float, "int": int, "bool": _str2bool}.get(tp, str)
        grp.add_argument(f"--{f.name}", type=kind, default=None,
                         metavar=tp.upper(),
                         help=f"{PARAM_DOC.get(f.name, '')} (default: {f.default})")
    return p


def build_config(args) -> Config:
    data = asdict(Config())
    if args.config:
        data.update(json.loads(Path(args.config).read_text(encoding="utf-8")))
    for f in fields(Config):
        v = getattr(args, f.name, None)
        if v is not None:
            data[f.name] = v
    return Config(**data)


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    if args.list_params:
        print(f"{'parameter':<20s}{'default':>14s}   description")
        print("-" * 86)
        for f in fields(Config):
            print(f"{f.name:<20s}{str(f.default):>14s}   {PARAM_DOC.get(f.name, '')}")
        return

    cfg = build_config(args)
    validate(cfg)

    if args.selftest:
        sys.exit(0 if self_tests(cfg) else 1)

    if not self_tests(cfg):
        print("Self-tests failed; aborting before the production run.")
        sys.exit(1)

    print("Running the internal-ballistics simulation ...")
    res = run_model(cfg)
    hist = sample(res, cfg)
    summ = summarize(res, hist, cfg)

    outdir: Path = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    png = outdir / "ballistics.png"
    csvp = outdir / "time_history.csv"
    txtp = outdir / "summary.txt"
    jsnp = outdir / "summary.json"

    write_csv(hist, csvp)
    make_figure(res, hist, summ, cfg, png)
    text = summary_text(summ, cfg)
    txtp.write_text(text, encoding="utf-8")
    jsnp.write_text(json.dumps({"input": asdict(cfg), "results": summ},
                               indent=2, ensure_ascii=False) + "\n",
                    encoding="utf-8")

    print()
    print(text)
    print("Output files")
    for pth in (png, csvp, txtp, jsnp):
        print(f"  {pth}")


if __name__ == "__main__":
    main()