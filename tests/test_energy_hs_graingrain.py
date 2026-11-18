"""A-6 stage 2, C4: h_s of the burn surface is evaluated at the far temperature of the grain (contract 164 section 2b-4.3).

ASSUMED INPUTS, not data (held item H-06; sources: the A-marked values of fable/probes/S182/u182_arith.py, report 182
section 6.2): solid specific heat 1.0 / 1.3 / 1.5 kJ/(kg K), surface temperature 700 / 900 / 1000 K, far temperature
293 K, rho_s = 1800 kg/m3; the reference enthalpy -2e6 J/kg is an arbitrary synthetic value on the absolute basis.
"""
import inspect

import pytest

from gribs.energy.balance import ABSOLUTE_ENTHALPY_BASIS, h_s_at, h_s_for_burn_surface
from gribs.energy.config_model import BurnSurface, SolidEnthalpy

T_FAR, T_REF, H_REF = 293.0, 293.0, -2.0e6


def surface(cp=1300.0, t_far=T_FAR):
    return BurnSurface(form="a", grain_far_temperature_K=t_far, solid_density_kg_m3=1800.0,
                       solid_enthalpy=SolidEnthalpy(reference_J_kg=H_REF, reference_temperature_K=T_REF, cp_J_kgK=cp))


@pytest.mark.fast
def test_ha_the_default_evaluation_temperature_is_the_far_temperature_of_the_grain():
    assert h_s_for_burn_surface(surface(), products_enthalpy_basis=ABSOLUTE_ENTHALPY_BASIS) == H_REF
    assert h_s_for_burn_surface(surface(t_far=350.0), products_enthalpy_basis=ABSOLUTE_ENTHALPY_BASIS) == \
        pytest.approx(H_REF + 1300.0 * (350.0 - T_REF))


@pytest.mark.fast
def test_hb_evaluating_at_the_surface_temperature_overestimates_by_about_4e5_to_1e6_J_per_kg():
    diffs = []
    for cp in (1000.0, 1300.0, 1500.0):                                # assumed inputs (H-06)
        for t_s in (700.0, 900.0, 1000.0):
            solid = surface(cp).solid_enthalpy
            diffs.append(h_s_at(t_s, solid) - h_s_for_burn_surface(surface(cp), products_enthalpy_basis=ABSOLUTE_ENTHALPY_BASIS))
    assert min(diffs) == pytest.approx(1000.0 * (700.0 - 293.0)) and max(diffs) == pytest.approx(1500.0 * (1000.0 - 293.0))
    assert all(4.0e5 <= d <= 1.1e6 for d in diffs)                    # the order of the estimate of report 182


@pytest.mark.fast
def test_hc_the_enthalpy_basis_must_be_the_absolute_basis_of_the_products():
    with pytest.raises(ValueError, match="absolute"):
        h_s_for_burn_surface(surface(), products_enthalpy_basis="sensible_only")
    with pytest.raises(ValueError, match="absolute"):
        h_s_for_burn_surface(surface(), products_enthalpy_basis="")


@pytest.mark.fast
def test_hd_there_is_no_argument_for_the_surface_temperature():
    assert list(inspect.signature(h_s_for_burn_surface).parameters) == ["surface", "products_enthalpy_basis"]
    shift = (h_s_for_burn_surface(surface(t_far=393.0), products_enthalpy_basis=ABSOLUTE_ENTHALPY_BASIS)
             - h_s_for_burn_surface(surface(t_far=293.0), products_enthalpy_basis=ABSOLUTE_ENTHALPY_BASIS))
    assert shift == pytest.approx(1300.0 * 100.0)


@pytest.mark.fast
def test_he_the_a_b_form_difference_p_over_rho_s_is_two_orders_below_the_overestimate():
    p_over_rho_s = 7.0e6 / 1800.0                                     # J/kg, the difference between the forms (not normalised)
    overestimate_mid = 1300.0 * (900.0 - 293.0)                       # the middle of the assumed inputs
    assert p_over_rho_s / overestimate_mid < 1.0 / 100.0
