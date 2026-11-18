"""A-6 stage 2, C1: the mass-basis relations of contract 164 section 1 (F-21), as pure functions.

The numbers are the F-21 example of the contract (golden g5 at 2 MPa: R_eff = 309.31 J/(kg-total K),
Ru/M_gas = 437.73 J/(kg-gas K), hence Yg = 0.7066); the tolerance of 1e-3 (relative) is the number of
significant digits those two values are quoted with.  Not a validation of any propellant.
"""
import pytest

from gribs.constants import R_UNIVERSAL
from gribs.energy import relations

R_EFF, R_GAS, YG = 309.31, 437.73, 0.7066


@pytest.mark.fast
def test_n01_r_eff_is_yg_times_r_gas():
    assert relations.r_gas(R_EFF, YG) == pytest.approx(R_GAS, rel=1e-3)
    assert YG * relations.r_gas(R_EFF, YG) == pytest.approx(R_EFF, rel=1e-12)
    with pytest.raises(ValueError):
        relations.r_gas(R_EFF, 0.0)


@pytest.mark.fast
def test_n02_m_gas_is_yg_times_m_cea():
    m_cea = R_UNIVERSAL / R_EFF                     # R_eff = Ru / M_CEA
    mg = relations.m_gas(m_cea, YG)
    assert mg == pytest.approx(YG * m_cea, rel=1e-12)
    assert R_UNIVERSAL / mg == pytest.approx(R_GAS, rel=1e-3)   # R_gas = Ru / M_gas


@pytest.mark.fast
def test_n03_psi_is_r_eff_times_t():
    assert relations.psi(R_EFF, 3000.0) == pytest.approx(R_EFF * 3000.0, rel=1e-15)
