"""A-6 stage 2, C4: the first-law balance, the burn-surface term and H0_out as pure functions (contract 164 sections 1b,
2b-4.3, 2b-4.4; PRE_REGISTRATION_C4).  Synthetic fixtures; the numbers are not data.  No assignment, no monkeypatch."""
import dataclasses
import inspect
import math

import pytest

from gribs.energy import balance
from gribs.energy.balance import (BalanceTerms, Flow, burn_surface_term, first_law_rhs, h0_out_from_chamber, h_stag)

# a burn surface that satisfies the interface identities: mdot = rho_s A_b r and h_s = u_s + p_s / rho_s
RHO_S, A_B, R_BURN, P_S, U_S = 1800.0, 0.01, 0.005, 7.0e6, -2.1e6
MDOT = RHO_S * A_B * R_BURN
H_S = U_S + P_S / RHO_S


@pytest.mark.fast
def test_f01_a_closed_adiabatic_moving_piston_gives_dU_equal_minus_integral_p_dV():
    m, cv, R, V0, T0 = 0.05, 900.0, 300.0, 1.0e-3, 1000.0           # ideal gas, constant cv; V grows from V0 to 2 V0
    dVdt = V0

    def dEdt(E, V):
        p = R * E / (cv * V)                                          # p = m R T / V with T = E / (m cv)
        return first_law_rhs(BalanceTerms(Qdot_in=0.0, Wdot_pressure_boundary_out=p * dVdt))

    E, V, n = m * cv * T0, V0, 2000
    E_start, h = E, 1.0 / n
    for _ in range(n):                                                # RK4 in time, t from 0 to 1
        k1 = dEdt(E, V)
        k2 = dEdt(E + 0.5 * h * k1, V + 0.5 * h * dVdt)
        k3 = dEdt(E + 0.5 * h * k2, V + 0.5 * h * dVdt)
        k4 = dEdt(E + h * k3, V + h * dVdt)
        E, V = E + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0, V + h * dVdt
    T_exact = T0 * (V0 / (2.0 * V0)) ** (R / cv)                      # isentropic: T V^(gamma-1) = constant
    assert E < E_start                                                # the expansion work is taken out
    assert E - E_start == pytest.approx(m * cv * (T_exact - T0), rel=1e-6)


@pytest.mark.fast
def test_f02_constant_volume_heating_adds_q_times_t():
    rhs = first_law_rhs(BalanceTerms(Qdot_in=500.0))
    assert rhs == 500.0 and 1.0e5 + rhs * 2.0 == 1.0e5 + 1000.0


@pytest.mark.fast
def test_f03_inflow_adds_outflow_subtracts_and_a_negative_flow_is_rejected():
    t = BalanceTerms(Qdot_in=0.0, inflows=(Flow(0.5, -1.0e6),), outflows=(Flow(0.2, -2.0e6),))
    assert first_law_rhs(t) == pytest.approx(0.5 * -1.0e6 - 0.2 * -2.0e6)
    for bad in (-0.1, math.nan, math.inf):
        with pytest.raises(ValueError):
            Flow(bad, 1.0e6)


@pytest.mark.fast
def test_f04_boundary_pressure_work_is_positive_out_and_zero_for_fixed_walls():
    assert first_law_rhs(BalanceTerms(Qdot_in=0.0, Wdot_pressure_boundary_out=100.0)) == -100.0   # expansion
    assert first_law_rhs(BalanceTerms(Qdot_in=0.0)) == 0.0                                       # fixed walls
    assert first_law_rhs(BalanceTerms(Qdot_in=0.0, Wdot_pressure_boundary_out=-40.0)) == 40.0    # compression


@pytest.mark.fast
def test_f05_there_is_no_separate_chemical_heat_term():
    names = {f.name for f in dataclasses.fields(BalanceTerms)}
    assert names == {"Qdot_in", "inflows", "outflows", "burn_surface_net", "Wdot_pressure_boundary_out", "Wdot_other_out"}
    assert not [n for n in names if "chem" in n.lower()]


@pytest.mark.fast
def test_f06_forms_a_and_b_agree_when_the_interface_identities_hold():
    a = burn_surface_term(MDOT, H_S, P_S, A_B, R_BURN, form="a")
    b = burn_surface_term(MDOT, H_S, P_S, A_B, R_BURN, form="b", u_s=U_S, rho_s=RHO_S)
    assert a == pytest.approx(b, rel=1e-12) and b == pytest.approx(MDOT * U_S, rel=1e-15)


@pytest.mark.fast
def test_f07_form_b_is_rejected_when_an_identity_is_broken():
    with pytest.raises(ValueError, match="mdot = rho_s"):
        burn_surface_term(1.1 * MDOT, H_S, P_S, A_B, R_BURN, form="b", u_s=U_S, rho_s=RHO_S)
    with pytest.raises(ValueError, match="h_s = u_s"):
        burn_surface_term(MDOT, H_S + 5.0e3, P_S, A_B, R_BURN, form="b", u_s=U_S, rho_s=RHO_S)
    with pytest.raises(ValueError, match="needs u_s and rho_s"):
        burn_surface_term(MDOT, H_S, P_S, A_B, R_BURN, form="b")


@pytest.mark.fast
def test_f08_exactly_one_form_per_call():
    with pytest.raises(ValueError, match="form"):
        burn_surface_term(MDOT, H_S, P_S, A_B, R_BURN, form="ab")
    assert isinstance(burn_surface_term(MDOT, H_S, P_S, A_B, R_BURN), float)       # one number, never a pair
    assert [f.name for f in dataclasses.fields(BalanceTerms)].count("burn_surface_net") == 1


@pytest.mark.fast
def test_f09_h0_out_is_the_stagnation_enthalpy_at_the_inlet_plane():
    h = -1.0e6
    assert h0_out_from_chamber(h, 100.0) == pytest.approx(h + 5000.0)
    assert h0_out_from_chamber(h, 100.0) != h and h0_out_from_chamber(h, 0.0) == h


@pytest.mark.fast
def test_f10_kinetic_energy_sensitivity_is_v_squared_over_two():
    assert h0_out_from_chamber(-1.0e6, 200.0) - h0_out_from_chamber(-1.0e6, 100.0) == pytest.approx(0.5 * (200.0 ** 2 - 100.0 ** 2))


@pytest.mark.fast
def test_f11_slip_is_rejected():
    with pytest.raises(ValueError, match="slip"):
        h0_out_from_chamber(-1.0e6, 100.0, slip=True)


@pytest.mark.fast
def test_f12_h_stag_units_include_g_times_z():
    assert h_stag(-1.0e6, 0.0, z=10.0) == pytest.approx(-1.0e6 + 98.0665)
    assert h_stag(-1.0e6, 3.0) == -1.0e6 + 4.5


@pytest.mark.fast
def test_f13_there_is_no_storage_term_for_ignition_and_tail_off():
    params = set(inspect.signature(burn_surface_term).parameters)
    assert params == {"mdot_burn", "h_s", "p_s", "A_b", "r", "form", "u_s", "rho_s"}
    assert not [n for n in params | {f.name for f in dataclasses.fields(BalanceTerms)} if "storage" in n.lower()]
    assert not hasattr(balance, "dE_solid_layer")


@pytest.mark.fast
def test_f14_heat_into_the_chamber_raises_dE_dt_and_a_loss_lowers_it():
    assert first_law_rhs(BalanceTerms(Qdot_in=250.0)) == 250.0
    assert first_law_rhs(BalanceTerms(Qdot_in=-250.0)) == -250.0
