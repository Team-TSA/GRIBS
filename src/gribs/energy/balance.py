"""First-law balance of the chamber control volume, as pure functions (A-6 stage 2, C4; contract 164 sections 1b, 2b-4.3, 2b-4.4).

dE_cv/dt = Qdot_in + sum_in mdot_rel*H0_in - sum_out mdot_rel*H0_out - Wdot_pressure_boundary_out - Wdot_other_out
(plus the burn-surface term below).  E_cv is the total energy of the contents (internal, kinetic, potential).  mdot_rel
is the mass flow relative to each face (never negative).  H0 = h + v^2/2 + g z with the velocity in an inertial frame.
Wdot_pressure_boundary_out = the integral of p v_b.n over the moving boundary (positive when the volume grows; zero for
fixed walls).  Wdot_other_out is positive outward.  Chemical energy enters through the formation enthalpy in h; there is
no separate chemical heat term.  Heat through the burn surface is included by evaluating h_s at the far temperature of the
grain (T_grain_inf), not added separately.  Quasi-steady only: there is no storage term for ignition and tail-off.

Standard library and gribs.energy.config_model only.  Nothing imports this module yet.
"""
import math
from dataclasses import dataclass
from typing import Optional, Tuple

from gribs.energy.config_model import BurnSurface, SolidEnthalpy

#: The enthalpy basis h_s must share with the products: absolute, formation enthalpy included (contract 2b-4.3 (iii)).
ABSOLUTE_ENTHALPY_BASIS = "absolute_formation_included"
STANDARD_GRAVITY = 9.80665        # [m/s^2]
_BURN_SURFACE_FORMS = ("a", "b")
_IDENTITY_REL_TOL = 1.0e-9


def _finite(name: str, value: float) -> float:
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite, got {value!r}")
    return value


@dataclass(frozen=True)
class Flow:
    """One flow across a face: the mass flow relative to the face [kg/s, never negative] and its stagnation enthalpy [J/kg]."""
    mdot_rel: float
    H0: float

    def __post_init__(self):
        _finite("mdot_rel", self.mdot_rel)
        _finite("H0", self.H0)
        if self.mdot_rel < 0.0:
            raise ValueError(f"mdot_rel must not be negative, got {self.mdot_rel!r} (the sign is carried by in / out)")


@dataclass(frozen=True)
class BalanceTerms:
    """The terms of dE_cv/dt [W].  burn_surface_net is the NET energy rate into the control volume across the burn surface
    (the enthalpy inflow and the pressure work of the receding surface together, from burn_surface_term): when it is used,
    mdot*h_s must not also be listed in inflows and p_s*A_b*r must not also be part of Wdot_pressure_boundary_out."""
    Qdot_in: float
    inflows: Tuple[Flow, ...] = ()
    outflows: Tuple[Flow, ...] = ()
    burn_surface_net: float = 0.0
    Wdot_pressure_boundary_out: float = 0.0
    Wdot_other_out: float = 0.0


def first_law_rhs(t: BalanceTerms) -> float:
    """dE_cv/dt [W] from the terms."""
    for name in ("Qdot_in", "burn_surface_net", "Wdot_pressure_boundary_out", "Wdot_other_out"):
        _finite(name, getattr(t, name))
    return (t.Qdot_in + t.burn_surface_net
            + sum(f.mdot_rel * f.H0 for f in t.inflows)
            - sum(f.mdot_rel * f.H0 for f in t.outflows)
            - t.Wdot_pressure_boundary_out - t.Wdot_other_out)


def h_stag(h: float, v: float, z: float = 0.0, g: float = STANDARD_GRAVITY) -> float:
    """H0 = h + v^2/2 + g z [J/kg] (v in an inertial frame)."""
    return _finite("h", h) + 0.5 * _finite("v", v) ** 2 + g * _finite("z", z)


def h_s_at(T: float, solid: SolidEnthalpy) -> float:
    """Specific enthalpy of the solid at temperature T [J/kg]: reference + cp (T - T_ref) (constant cp)."""
    if not (math.isfinite(T) and T > 0.0):
        raise ValueError(f"T must be a positive finite temperature, got {T!r}")
    return solid.reference_J_kg + solid.cp_J_kgK * (T - solid.reference_temperature_K)


def h_s_for_burn_surface(surface: BurnSurface, *, products_enthalpy_basis: str) -> float:
    """h_s of the burn surface [J/kg], always evaluated at the far temperature of the grain (T_grain_inf, contract 2b-4.3):
    the quasi-steady balance of the preheated layer gives the net inflow m_dot h_s(T_grain_inf) whatever the heat path is.
    There is no argument for the surface temperature.  products_enthalpy_basis is the enthalpy convention the products'
    provider declares; it must be the absolute basis with the formation enthalpy included."""
    if products_enthalpy_basis != ABSOLUTE_ENTHALPY_BASIS:
        raise ValueError(
            f"h_s must share the absolute enthalpy basis of the products ({ABSOLUTE_ENTHALPY_BASIS!r}, formation enthalpy "
            f"included); the products declare {products_enthalpy_basis!r}")
    return h_s_at(surface.grain_far_temperature_K, surface.solid_enthalpy)


def burn_surface_term(mdot_burn: float, h_s: float, p_s: float, A_b: float, r: float, form: str = "a",
                      u_s: Optional[float] = None, rho_s: Optional[float] = None) -> float:
    """Net energy rate into the control volume across the burn surface [W]; exactly one form per call (never both).

    form "a": mdot*h_s - p_s*A_b*r  (the pressure work of the receding surface, p_s A_b r, is part of the term).
    form "b": mdot*u_s, valid only when mdot = rho_s*A_b*r and h_s = u_s + p_s/rho_s hold (the interface state is the
    same on both sides); both identities are checked (relative 1e-9) and a violation is an error."""
    if form not in _BURN_SURFACE_FORMS:
        raise ValueError(f"form must be one of {list(_BURN_SURFACE_FORMS)}, got {form!r}")
    for name, value in (("mdot_burn", mdot_burn), ("h_s", h_s), ("p_s", p_s), ("A_b", A_b), ("r", r)):
        _finite(name, value)
    if form == "a":
        return mdot_burn * h_s - p_s * A_b * r
    if u_s is None or rho_s is None:
        raise ValueError("form 'b' needs u_s and rho_s (the identities of the interface are checked)")
    _finite("u_s", u_s)
    if not (math.isfinite(rho_s) and rho_s > 0.0):
        raise ValueError(f"rho_s must be positive and finite, got {rho_s!r}")
    flux = rho_s * A_b * r
    if abs(mdot_burn - flux) > _IDENTITY_REL_TOL * max(abs(mdot_burn), abs(flux), 1.0e-300):
        raise ValueError(f"identity mdot = rho_s*A_b*r violated: mdot = {mdot_burn!r}, rho_s*A_b*r = {flux!r}")
    h_from_u = u_s + p_s / rho_s
    if abs(h_s - h_from_u) > _IDENTITY_REL_TOL * max(abs(h_s), abs(h_from_u), 1.0):
        raise ValueError(f"identity h_s = u_s + p_s/rho_s violated: h_s = {h_s!r}, u_s + p_s/rho_s = {h_from_u!r}")
    return mdot_burn * u_s


def h0_out_from_chamber(h_static: float, v_inlet: float, *, slip: bool = False) -> float:
    """H0_out: the stagnation enthalpy at the nozzle inlet plane, h + v_inlet^2/2 with the STATIC enthalpy of the chamber
    gas at that plane (contract 2b-4.4 (a): no outlet species is needed).  The kinetic energy of the inlet velocity is
    included, and a stagnation value must not be passed as h_static.  Condition (iii): no slip or hold-up of a phase,
    which is rejected here."""
    if slip:
        raise ValueError("H0_out from the chamber state needs no slip or hold-up of a phase (contract 2b-4.4 (iii)); "
                         "a flux rule for the inlet plane is outside this contract")
    return h_stag(h_static, v_inlet)
