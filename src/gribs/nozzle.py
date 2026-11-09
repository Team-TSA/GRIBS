"""Nozzle model (A-6 unit 4): isentropic reference flow, choking limit and the two-phase mass split. Moved verbatim
from section 5 of the entry. Tests patch nozzle_state, choke_limit_pressure and NozzleTwoPhaseModel.split here,
e.g. "gribs.nozzle.nozzle_state"; the thermochemistry is read as gribs.thermo.<name>.
"""
from __future__ import annotations

import math
from typing import Optional, Tuple

from scipy.optimize import brentq

import gribs.thermo
from gribs.config import Config
from gribs.constants import KNOWN_TWO_PHASE_MODES, TWO_PHASE_HOMOGENEOUS, TWO_PHASE_SINGLE_LEGACY
from gribs.errors import ConfigurationError, ThermochemistryError


class NozzleTwoPhaseModel:
    """Two-phase mass split of the nozzle flow (isolated for future upgrades).

    Implemented model: HOMOGENEOUS EQUILIBRIUM, COMPLETE ENTRAINMENT.
      * The isentropic formula is evaluated with the effective gas constant
        R_eff = Ru / M_CEA = Yg * R_gas, so it returns the TOTAL-mixture
        flow directly.  No additional /Yg factor is applied anywhere.
      * Both phases are assumed to leave with the same velocity, so the
        total flow is decomposed with the chamber phase split:
            mdot_total     = mdot_total_reference
            mdot_gas       = Yg_chamber * mdot_total
            mdot_condensed = Yc_chamber * mdot_total
      * This is a FIRST APPROXIMATION: it is not a rigorous two-phase choking
        model (no particle slip, no two-phase sound speed, no non-equilibrium
        nozzle chemistry, no particle-size effects).  For large condensed
        fractions it tends to OVERESTIMATE the delivered mass flow and thrust.
      * In 'single_phase_legacy' mode Yg := 1, so R_eff is the gas constant
        itself and the total flow equals the gas reference flow, i.e. exactly
        the v0.4.3-alpha behaviour (regression comparison only).

    A future rigorous model (e.g. frozen/non-equilibrium two-phase choking,
    particle lag) should replace ``split`` without touching the ODE, which only
    consumes (mdot_total, mdot_gas, mdot_condensed).
    """

    def __init__(self, mode: str):
        if mode not in KNOWN_TWO_PHASE_MODES:
            raise ConfigurationError(
                f"Unknown two-phase mode {mode!r}; expected one of "
                f"{list(KNOWN_TWO_PHASE_MODES)}.")
        self.mode = mode

    @property
    def entrainment_assumption(self) -> str:
        return ("complete entrainment: the isentropic formula evaluated with "
                "R_eff returns mdot_total directly (no additional /Yg); both "
                "phases share the exit velocity (no particle slip)"
                if self.mode == TWO_PHASE_HOMOGENEOUS else
                "single-phase legacy: the total flow equals the gas reference "
                "flow (Yg := 1)")

    def split(self, mdot_total_reference: float,
              Yg_chamber: float) -> Tuple[float, float, float]:
        """Return (mdot_total, mdot_gas, mdot_condensed) [kg/s].

        ``mdot_total_reference`` is the total-mixture flow produced by the
        isentropic formula evaluated with R_eff = Ru/M_CEA = Yg*R_gas.  It
        already carries the phase split through R_eff, so it is returned
        unchanged as ``mdot_total``: no additional /Yg factor is applied.
        """
        if self.mode == TWO_PHASE_SINGLE_LEGACY:
            return mdot_total_reference, mdot_total_reference, 0.0
        if not (0.0 < Yg_chamber <= 1.0):
            raise ThermochemistryError(
                f"Nozzle two-phase split received Yg_chamber = {Yg_chamber!r}; "
                "require 0 < Yg <= 1.")
        mdot_total = mdot_total_reference
        mdot_gas = Yg_chamber * mdot_total
        return mdot_total, mdot_gas, mdot_total - mdot_gas


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
    _, _, g = gribs.thermo.gas_props(p0, c)
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
    """Nozzle mass flow, thrust, exit conditions and flow regime.

    MASS BASIS (v0.5.0-alpha): the returned ``mdot_total`` is the TOTAL product
    mass flow (gas + condensed) and is the quantity that must be fed to the
    pressure ODE, which conserves the total product mass.  ``mdot_gas`` and
    ``mdot_condensed`` decompose it with the chamber equilibrium phase split.
    ``mdot`` is retained as a backward-compatible ALIAS of ``mdot_total`` (it is
    NOT the gas-only flow); prefer the explicit names in new code.

    The isentropic reference flow is evaluated with the effective gas
    constant R_eff(p0) = Ru / M_CEA(p0) = Yg(p0) * R_gas(p0), so it returns
    the TOTAL-mixture flow directly with no additional Yg factor (complete
    entrainment, see ``NozzleTwoPhaseModel``).  In the single-phase legacy
    mode Yg := 1, so R_eff is the gas constant itself and the total flow
    equals the gas reference flow.
    """
    split_model = NozzleTwoPhaseModel(c.two_phase_mode)
    Yg, Yc = gribs.thermo.phase_fractions(p0, c)

    if p0 <= c.p_a:
        return dict(mdot=0.0,                   # alias of mdot_total (documented)
                    mdot_total=0.0,
                    mdot_gas=0.0,
                    mdot_condensed=0.0,
                    mdot_total_reference=0.0,
                    gas_mass_fraction=Yg,
                    condensed_mass_fraction=Yc,
                    F=0.0, Me=0.0, pe=p0, ve=0.0,
                    regime="no-flow", A_eff=Ae)

    R, T0, g = gribs.thermo.gas_props(p0, c)
    eps = Ae / At
    sq = math.sqrt(g / (R * T0))
    p_thr = choke_threshold_pressure(p0, At, Ae, c)

    if p0 >= p_thr:                                   # ---- choked ----
        mdot_ref = c.Cd * At * p0 * sq \
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
        # Complete velocity equilibrium: both phases leave at ve_eff.  With
        # particle slip the momentum term would be mdot_gas*vg + mdot_c*vc and
        # the thrust would be lower for large condensed fractions; eta_F may
        # additionally lump two-phase losses (avoid double correction).
        mdot_total, mdot_gas, mdot_cond = split_model.split(mdot_ref, Yg)
        F = c.eta_thrust * mdot_total * ve_eff + (pe_eff - c.p_a) * A_eff
        return dict(mdot=mdot_total,
                    mdot_total=mdot_total,
                    mdot_gas=mdot_gas,
                    mdot_condensed=mdot_cond,
                    mdot_total_reference=mdot_ref,
                    gas_mass_fraction=Yg,
                    condensed_mass_fraction=Yc,
                    F=F, Me=Me, pe=pe_eff, ve=ve_eff,
                    regime=regime, A_eff=A_eff)

    # ---- subsonic (exit pressure equals ambient); continuous with the above ----
    Me = math.sqrt(2.0 / (g - 1.0) * ((p0 / c.p_a) ** ((g - 1.0) / g) - 1.0))
    mdot_ref = c.Cd * Ae * p0 * sq * Me \
        * (1.0 + 0.5 * (g - 1.0) * Me * Me) ** (-(g + 1.0) / (2.0 * (g - 1.0)))
    Te = T0 / (1.0 + 0.5 * (g - 1.0) * Me * Me)
    ve = Me * math.sqrt(g * R * Te)
    mdot_total, mdot_gas, mdot_cond = split_model.split(mdot_ref, Yg)
    F = c.eta_thrust * mdot_total * ve
    return dict(mdot=mdot_total,
                mdot_total=mdot_total,
                mdot_gas=mdot_gas,
                mdot_condensed=mdot_cond,
                mdot_total_reference=mdot_ref,
                gas_mass_fraction=Yg,
                condensed_mass_fraction=Yc,
                F=F, Me=Me, pe=c.p_a, ve=ve,
                regime="subsonic", A_eff=Ae)
