"""A-6 stage 2, C5: state recovery across a phase-transition plateau, with a SYNTHETIC provider (contract 164 section
2b-FIX-1; PRE_REGISTRATION_C5_REV1).  Every number is an ASSUMED INPUT, not data (held items H-01, H-02, H-08): the
densities 3990 / 3000 kg/m3 and T_m = 2327 K are the A-marked values of fable/probes/S182/u182_arith.py; L = 1.0e6 J/kg
and the gas data are synthetic.  The providers are frozen dataclasses passed as arguments: no monkeypatch, no assignment.

Tolerances (approval, 7th record, B-4, provisional): the phase fraction must agree with the exact coupled solution within
5 x delta_ln_p (delta_ln_p = ln((V - Vc_solid)/(V - Vc_liquid)), the change of ln p along the plateau that melting
causes), and the plateau edges are met to 1e-6 of the plateau width.
"""
import dataclasses
import math

import pytest
from scipy.optimize import brentq

from gribs.constants import R_UNIVERSAL
from gribs.energy import recovery
from gribs.energy.config_model import CondensedData, CondensedPhaseSpec, PhaseTransitionSpec, RecoverySettings
from gribs.energy.exceptions import NoRootOutsidePlateauError, StateRecoveryError
from gribs.energy.provider import TPState
from gribs.energy.recovery import latent_heat, plateau_bounds, plateau_phase_fraction, recover_state

RU, M = R_UNIVERSAL, 1.0
B_E = {"X": M}
T_M = 2327.0
SOLID = CondensedPhaseSpec("OX(s)", "solid", 3990.0, 1300.0, -1.6e7, 298.15)
LIQUID = CondensedPhaseSpec("OX(L)", "liquid", 3000.0, 1300.0, -1.5e7, 298.15)
TRANSITION = PhaseTransitionSpec("OX(s)", "OX(L)", T_M)
CONDENSED = CondensedData("external_density", (SOLID, LIQUID), (TRANSITION,))
NO_TRANSITION = CondensedData("external_density", (SOLID, LIQUID), ())
RECOVERY = RecoverySettings(phase_fraction_tolerance=1.0e-6, root_tolerance=1.0e-10, plateau_edge_tolerance=1.0e-6)
SEARCH, STARTS = (300.0, 6000.0), (600.0, 1500.0, 3500.0, 5000.0)
YG, MG, CV_G, U_G0, A_P = 0.6, 28.0, 900.0, -3.0e6, 1.0e-3     # gas fraction, molar mass, cv, offset, weak p dependence
VC_S, VC_L = (1.0 - YG) * M / SOLID.density_kg_m3, (1.0 - YG) * M / LIQUID.density_kg_m3
V_FREE = M * (YG / MG) * RU * T_M / 3.0e6 + VC_S                # about 3 MPa at T_m with the condensed phase solid
DLNP = math.log((V_FREE - VC_S) / (V_FREE - VC_L))              # change of ln p along the plateau
PHI_TOL = 5.0 * DLNP                                            # the provisional tolerance of the phase fraction


def h_of(phase, T):
    return phase.enthalpy_reference_J_kg + phase.cp_J_kgK * (T - phase.reference_temperature_K)


@dataclasses.dataclass(frozen=True)
class Melting:
    """Ideal gas (constant cv, a weak a*p term in u) plus one condensed species that is solid below T_m and liquid above."""
    phase_at_Tm: str = "solid"
    converged: bool = True
    jump_declared: bool = True          # False: the same jump, but the transition is not declared to the recovery

    def tp(self, T, p, B_e):
        liquid = T > T_M or (T == T_M and self.phase_at_Tm == "liquid")
        phase = LIQUID if liquid else SOLID
        yc = 1.0 - YG
        u = YG * (U_G0 + CV_G * (T - 300.0) + A_P * p) + yc * h_of(phase, T)
        return TPState(u_J_per_kg=u, n_gas_per_kg=YG / MG, Yg=YG, condensed_mass_fraction={phase.species: yc},
                       converged=self.converged)


@dataclasses.dataclass(frozen=True)
class Gradual:
    """Condensation that grows continuously as T falls (an H2O-like species): no jump of the energy, no plateau."""
    def tp(self, T, p, B_e):
        yc = 0.3 * min(max((1500.0 - T) / 400.0, 0.0), 1.0)
        yg = 1.0 - yc
        u = yg * (U_G0 + CV_G * (T - 300.0)) + yc * (-1.0e7 + 1300.0 * (T - 300.0))
        return TPState(u_J_per_kg=u, n_gas_per_kg=yg / MG, Yg=yg, condensed_mass_fraction={"OX(s)": yc})


@dataclasses.dataclass(frozen=True)
class ThreeRoots:
    """u(T) that is not monotone: it falls between about 1289 K and 1711 K, so E = 1.0e6 J has three roots (1135, 1500 and
    1865 K); a start below the first root and a start between the second and third root meet different upward crossings."""
    def tp(self, T, p, B_e):
        u = 1.0e6 - 1000.0 * (T - 1500.0) + 300.0 * (T - 1500.0) ** 3 / 200.0 ** 2
        return TPState(u_J_per_kg=u, n_gas_per_kg=1.0 / MG, Yg=1.0)


def energy(provider, T, V=V_FREE):
    """M u(T, p(T)) of the provider at temperature T, at the volume V (what the recovery inverts)."""
    return recovery._Solver(0.0, V, B_E, provider, CONDENSED, RECOVERY).residual(T)


def bounds(provider=None):
    return plateau_bounds(V_FREE, B_E, TRANSITION, provider or Melting(), CONDENSED, RECOVERY)


def recover(E, provider=None, condensed=CONDENSED, starts=STARTS, **kw):
    return recover_state(E, V_FREE, B_E, provider or Melting(), condensed, RECOVERY, T_search_K=SEARCH, starts_K=starts, **kw)


def exact_phi(E):
    """The coupled exact solution on the plateau (phi and p from the volume and the energy equations), found by a solver
    that does not use the linear approximation."""
    def p_of(phi):
        return M * (YG / MG) * RU * T_M / (V_FREE - (1.0 - YG) * M * ((1.0 - phi) / SOLID.density_kg_m3 + phi / LIQUID.density_kg_m3))

    def residual(phi):
        u = YG * (U_G0 + CV_G * (T_M - 300.0) + A_P * p_of(phi)) + (1.0 - YG) * ((1.0 - phi) * h_of(SOLID, T_M) + phi * h_of(LIQUID, T_M))
        return M * u - E
    return brentq(residual, 0.0, 1.0, xtol=1.0e-15, rtol=1.0e-14)


def width_and_center():
    b = bounds()
    return b.E_plus - b.E_minus, 0.5 * (b.E_plus + b.E_minus)


@pytest.mark.fast
def test_r01_ta_known_latent_heat_gives_the_phase_fraction_of_the_exact_solution_within_5_dlnp():
    width, _ = width_and_center()
    b = bounds()
    for phi_true in (0.1, 0.35, 0.5, 0.8):
        E = b.E_minus + phi_true * width                              # an energy on the plateau, chosen by position
        got = recover(E)
        assert got.on_plateau and got.T_K == T_M
        assert abs(got.phi - exact_phi(E)) <= PHI_TOL, (phi_true, got.phi, exact_phi(E), PHI_TOL)
    assert latent_heat(TRANSITION, CONDENSED.phases) == pytest.approx(1.0e6, rel=1e-12)    # L from the two phases' data
    assert 2.0e-4 < DLNP < 3.0e-4                                      # the synthetic change of ln p along the plateau


@pytest.mark.fast
def test_r02_the_phase_fraction_is_the_linear_approximation_and_close_to_the_coupled_solution():
    b = bounds()
    width = b.E_plus - b.E_minus
    E = b.E_minus + 0.4 * width
    assert plateau_phase_fraction(E, b.E_minus, b.E_plus) == pytest.approx(0.4, rel=1e-12)
    assert recover(E).phi == pytest.approx(0.4, rel=1e-12)             # the recovery reports the linear approximation
    deviation = abs(plateau_phase_fraction(E, b.E_minus, b.E_plus) - exact_phi(E))
    assert deviation > 0.0 and deviation <= PHI_TOL                    # an approximation (a*p couples), within tolerance
    with pytest.raises(StateRecoveryError):
        plateau_phase_fraction(E, 1.0, 1.0)


@pytest.mark.fast
def test_r03_tb_heating_through_the_plateau_is_continuous():
    width, center = width_and_center()
    sweep = [recover(center + (k - 20) * width / 10.0) for k in range(41)]      # E from -2 width to +2 width
    for prev, cur in zip(sweep, sweep[1:]):
        assert cur.T_K >= prev.T_K - 1.0e-6
        assert cur.T_K - prev.T_K <= 1.2 * (width / 10.0) / 1060.0 + 1.0e-6     # no jump of T (cp of the mixture is 1060 J/K)
        assert abs(math.log(cur.p_Pa / prev.p_Pa)) <= math.log(cur.T_K / prev.T_K) + 5.0 * DLNP + 1.0e-6   # p ~ T / (V - Vc)
    phis = [s.phi for s in sweep if s.on_plateau]
    assert phis == sorted(phis) and 0.0 <= phis[0] < 0.2 and 0.8 < phis[-1] <= 1.0
    assert not sweep[0].on_plateau and not sweep[-1].on_plateau and sweep[0].T_K < T_M < sweep[-1].T_K


@pytest.mark.fast
def test_r04_tc_cooling_through_the_plateau_is_continuous():
    width, center = width_and_center()
    sweep = [recover(center + (20 - k) * width / 10.0) for k in range(41)]      # E from +2 width down to -2 width
    for prev, cur in zip(sweep, sweep[1:]):
        assert cur.T_K <= prev.T_K + 1.0e-6
        assert prev.T_K - cur.T_K <= 1.2 * (width / 10.0) / 1060.0 + 1.0e-6
    phis = [s.phi for s in sweep if s.on_plateau]
    assert phis == sorted(phis, reverse=True) and phis[0] > 0.8 and phis[-1] < 0.2


@pytest.mark.fast
def test_r05_td_no_root_outside_the_plateau_and_a_failed_equilibrium_are_errors():
    with pytest.raises(NoRootOutsidePlateauError, match="out of the searched range"):
        recover(-1.0e8)                                                # below the energy at the lowest searched T
    with pytest.raises(NoRootOutsidePlateauError, match="out of the searched range"):
        recover(1.0e9)                                                 # above the energy at the highest searched T
    with pytest.raises(StateRecoveryError, match="did not converge"):
        recover(center_energy(), provider=Melting(converged=False))


def center_energy():
    return width_and_center()[1]


@pytest.mark.fast
def test_r06_the_plateau_energy_is_not_rejected_as_no_root():
    got = recover(center_energy())
    assert got.on_plateau and got.T_K == T_M and 0.4 < got.phi < 0.6
    assert got.V_gas_m3 + got.V_condensed_m3 == pytest.approx(V_FREE, rel=1e-12)


@pytest.mark.fast
def test_r07_start_points_that_agree_are_accepted():
    got = recover(center_energy() - 0.0)
    assert got.n_starts_agreed == len(STARTS)
    off = recover(width_and_center()[1] - 3.0 * width_and_center()[0])
    assert not off.on_plateau and off.n_starts_agreed == len(STARTS) and off.phi is None


@pytest.mark.fast
def test_r08_start_points_that_disagree_are_an_error():
    E = 1.0e6                                                          # three roots: 1135, 1500 and 1865 K
    with pytest.raises(StateRecoveryError, match="disagree"):
        recover_state(E, V_FREE, B_E, ThreeRoots(), NO_TRANSITION, RECOVERY, T_search_K=SEARCH, starts_K=(600.0, 1700.0))
    one = recover_state(E, V_FREE, B_E, ThreeRoots(), NO_TRANSITION, RECOVERY, T_search_K=SEARCH, starts_K=(600.0,))
    assert one.n_starts_agreed == 1 and one.T_K == pytest.approx(1135.0, abs=1.0)   # a single start is compared with nothing


@pytest.mark.fast
def test_r09_the_phase_returned_at_exactly_tm_does_not_make_the_starts_disagree():
    E = center_energy()
    states = [recover(E, provider=Melting(phase_at_Tm=ph), starts=(T_M, 1500.0, 3500.0)) for ph in ("solid", "liquid")]
    for s in states:
        assert s.on_plateau and s.T_K == T_M and s.n_starts_agreed == 3
    a, b = states
    assert a.phi == pytest.approx(b.phi, abs=1.0e-9) and a.p_Pa == pytest.approx(b.p_Pa, rel=1e-9)


@pytest.mark.fast
def test_r10_a_transition_without_its_data_is_rejected_and_an_undeclared_jump_is_an_error():
    only_solid = CondensedData("external_density", (SOLID,), (TRANSITION,))
    with pytest.raises(StateRecoveryError, match="missing"):
        recover(center_energy(), condensed=only_solid)
    with pytest.raises(NoRootOutsidePlateauError, match="no phase transition is declared"):
        recover(center_energy(), condensed=NO_TRANSITION)             # the provider jumps; the recovery was told of no transition


@pytest.mark.fast
def test_r11_an_independent_latent_heat_must_match_the_two_phases():
    L = latent_heat(TRANSITION, CONDENSED.phases)
    assert recover(center_energy(), latent_heat_J_per_kg=L).on_plateau            # consistent: accepted
    with pytest.raises(StateRecoveryError, match="latent heat"):
        recover(center_energy(), latent_heat_J_per_kg=1.2 * L)


@pytest.mark.fast
def test_r12_a_mixture_condensation_is_continuous_and_makes_no_plateau():
    energies = [energy(Gradual(), T) for T in (900.0, 1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1500.0, 1600.0)]
    assert energies == sorted(energies)                                # continuous and increasing: no jump
    got = recover_state(energies[3], V_FREE, B_E, Gradual(), NO_TRANSITION, RECOVERY, T_search_K=SEARCH, starts_K=(600.0, 2500.0))
    assert not got.on_plateau and got.phi is None and got.T_K == pytest.approx(1200.0, rel=1e-6)
    assert got.condensed_mass_fraction["OX(s)"] > 0.0


@pytest.mark.fast
def test_r13_just_outside_an_edge_is_an_ordinary_root_and_just_inside_is_a_plateau_state():
    b = bounds()
    width = b.E_plus - b.E_minus
    below = recover(b.E_minus - 1.0e-3 * width)
    inside = recover(b.E_minus + 1.0e-3 * width)
    assert not below.on_plateau and below.T_K < T_M and below.phi is None
    assert inside.on_plateau and inside.T_K == T_M and inside.phi == pytest.approx(1.0e-3, abs=1.0e-9)
    above = recover(b.E_plus + 1.0e-3 * width)
    assert not above.on_plateau and above.T_K > T_M
    assert recover(b.E_plus - 1.0e-3 * width).phi == pytest.approx(1.0 - 1.0e-3, abs=1.0e-9)


@pytest.mark.fast
def test_r14_there_is_no_direct_uv_solve():
    import inspect
    assert not [n for n in dir(recovery) if "uv" in n.lower()]
    assert not [n for n in inspect.signature(recover_state).parameters if "uv" in n.lower()]
    assert "UV" in recovery.__doc__ and "no direct solve" in recovery.__doc__
