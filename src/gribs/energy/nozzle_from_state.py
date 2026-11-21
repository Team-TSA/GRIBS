"""The nozzle flow from a chamber state (A-6 stage 2, C8; contract 164 sections 3-C, 2b-4.4).

nozzle_state (legacy) reads T0, gamma and the phase split from the property TABLE at p0.  In energy mode the chamber state
has its own temperature, gas constant and phase split, so this function takes them as arguments and calls the same flow code
(gribs.nozzle._nozzle_flow, C2).  The legacy nozzle_state is unchanged: given the state of the table, this function returns
what nozzle_state returns.
"""
import math

import gribs.nozzle


def _choke_threshold_pressure(p_a: float, eps: float, gamma_s: float) -> float:
    """The chamber pressure above which the nozzle is choked, for the GIVEN gamma_s (the legacy choke_threshold_pressure
    reads gamma from the table at p0; the expression is the same)."""
    if eps <= 1.0 + 1.0e-12:
        return p_a / gribs.nozzle.critical_ratio(gamma_s)
    Me_sub = gribs.nozzle.mach_from_area(eps, gamma_s, supersonic=False)
    return p_a * (1.0 + 0.5 * (gamma_s - 1.0) * Me_sub ** 2) ** (gamma_s / (gamma_s - 1.0))


def nozzle_state_from_chamber(p_stag, T_stag, R_eff, gamma_s, Yg, At, Ae, c, Yc=None) -> dict:
    """Nozzle mass flow, thrust, exit conditions and flow regime for the chamber state (p_stag, T_stag, R_eff, gamma_s, Yg).

    p_stag is the chamber STAGNATION pressure (163 T16: the legacy code also uses p0 as a stagnation pressure) and T_stag the
    stagnation temperature; R_eff is the effective gas constant of the total mass (Ru / M_CEA = Yg * R_gas), so the mass flow
    comes out for the total mixture with no extra Yg factor, and it is proportional to p_stag / sqrt(R_eff * T_stag).  Cd scales
    the mass flow only: the stagnation enthalpy per unit mass that leaves is not changed by it (an adiabatic nozzle, held H-14).
    The keys of the result are those of nozzle_state.

    The closed forms below (choked flow, exit Mach number, exit pressure, exit velocity written with p_stag/rho_stag) close with
    gamma_s and the stagnation quantities; they hold for a homogeneous equilibrium flow (no velocity or temperature lag, shifting
    equilibrium), gamma_s constant along the isentrope, a path that does not cross a phase-transition plateau, the V0 convention
    (condensed volume 0) and a one-dimensional steady adiabatic flow.  R_eff does not have to be constant for them; it has to be
    constant for the forms that write the temperature ratio with the pressure ratio, T_e/T_0 = (p_e/p_0)^((gamma_s-1)/gamma_s),
    and for the exit velocity written as sqrt(2 gamma_s/(gamma_s-1) R_eff (T_0 - T_e)).  gamma_s is the quantity of RP-1311
    eq. 2.71, not c_p/c_v.  The exit velocity of the closed form equals sqrt(2 (h_0 - h_e)) only for a constant gamma_s; otherwise
    the latter is the right one."""
    if Yc is None:
        Yc = 1.0 - Yg
    split_model = gribs.nozzle.NozzleTwoPhaseModel(c.two_phase_mode)
    if p_stag <= c.p_a:
        return dict(mdot=0.0, mdot_total=0.0, mdot_gas=0.0, mdot_condensed=0.0, mdot_total_reference=0.0,
                    gas_mass_fraction=Yg, condensed_mass_fraction=Yc, F=0.0, Me=0.0, pe=p_stag, ve=0.0,
                    regime="no-flow", A_eff=Ae)
    eps = Ae / At
    sq = math.sqrt(gamma_s / (R_eff * T_stag))
    p_thr = _choke_threshold_pressure(c.p_a, eps, gamma_s)
    return gribs.nozzle._nozzle_flow(p_stag, At, Ae, c, split_model, Yg, Yc, R_eff, T_stag, gamma_s, eps, sq, p_thr)
