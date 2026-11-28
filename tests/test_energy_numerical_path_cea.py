"""NumericalPath on the golden g2 propellant (NH4ClO4 86 wt%, C4H6 14 wt%) with the official CEA package (without it the module is skipped).

The values are those of report 400 (run r646, CEA 3.3.4, bundled thermo.lib): T_ladder_K lower end 300 K, p_floor_ratio = 1e-5.  The
chamber is p_c = 1.5 MPa, HP equilibrium of the reactants at 298.15 K.  Tolerance 1e-6 relative (report 400 had rel_diff 0.0).
These are the values for a golden input, not a claim about a real propellant.  A condensed species that the energy model does not
declare is recorded by NumericalPath.last.undeclared_condensed (the first supersonic grid point of the Expansion cache with it).
"""
import dataclasses
import json
from pathlib import Path

import pytest

import gribs.thermo

cea = pytest.importorskip("cea")

from gribs.energy.cea_tp import CEATPProvider  # noqa: E402
from gribs.energy.numerical_path import DEFAULT_SETTINGS  # noqa: E402
from test_energy_numerical_path import NOT_PLATEAU, numeric_path  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
G2 = ROOT / "cases" / "golden" / "g2_single_phase_cea.json"
ATOMIC = {"N": 14.0067, "H": 1.00794, "Cl": 35.453, "O": 15.9994, "C": 12.0107}
FORMULAS = {"NH4CLO4(I)": {"N": 1, "H": 4, "Cl": 1, "O": 4}, "C4H6,butadiene": {"C": 4, "H": 6}}
P_C, AT = 1.5e6, 1.0e-3
SETTINGS = dataclasses.replace(DEFAULT_SETTINGS, T_ladder_K=(300.0, 400.0, 500.0, 800.0, 1000.0, 1200.0, 1500.0, 1800.0), p_floor_ratio=1.0e-5)
# report 400 section 5
S_C, T_C, P_T = 10756.02701074797, 2909.009889448683, 855540.5481174608
P_E_2000, T_E_2000, T_E_1000 = 23.198974256955125, 479.72565244376875, 516.7138537227204
T_FIRST_NONE = 506.0883413320999        # the lowest T of a grid point without a condensed species (eps = 1200)
TOL = 1.0e-6


def _rel(a, b):
    return abs(a - b) / abs(b)


@pytest.fixture(scope="module")
def chamber():
    names, wt = ["NH4CLO4(I)", "C4H6,butadiene"], [86.0, 14.0]       # the reactants of G2 (checked below)
    doc = json.loads(G2.read_text(encoding="utf-8"))
    assert [(r["name"], r["wt_percent"]) for r in doc["propellant"]["reactants"]] == list(zip(names, wt))
    elements = {e: 0.0 for e in ATOMIC}
    for name, w in zip(names, wt):
        mass = sum(ATOMIC[e] * n for e, n in FORMULAS[name].items())
        for e, n in FORMULAS[name].items():
            elements[e] += w / 100.0 * ATOMIC[e] * n / mass
    import numpy as np
    module = gribs.thermo.import_official_cea_module()
    gribs.thermo.pin_bundled_cea_database(module, False)
    reactants = module.Mixture(names)
    products = module.Mixture(names, products_from_reactants=True)
    solver = module.EqSolver(products, reactants=reactants, transport=False, ions=False, trace=1.0e-10, insert=[])
    solution = module.EqSolution(solver)
    weights = np.array(wt, dtype=float)
    h_c = float(reactants.calc_property(module.ENTHALPY, weights, np.array([298.15] * 2)))
    solver.solve(solution, module.HP, h_c / float(module.R), P_C / 1.0e5, weights)
    assert bool(solution.converged)
    return CEATPProvider(names, wt, elements), {e: 2.0 * f for e, f in elements.items()}, float(solution.T)


def _run(chamber, eps):
    prov, B_e, T_c = chamber
    c, path = numeric_path(prov, SETTINGS, p_a=1.0e3, Cd=1.0, eta_thrust=1.0)
    r = path._numerical(NOT_PLATEAU, P_C, T_c, 1.0, AT, eps * AT, c, B_e, (T_c, P_C, eps), None, "test")
    assert r.numeric
    return r, path._exp.state_at(r.p_e).T_K, path


def test_c01_g2_chamber_and_throat_pressure_match_report_400(chamber):
    prov, B_e, T_c = chamber
    st = prov.tp(T_c, P_C, B_e)
    assert _rel(T_c, T_C) <= TOL and _rel(st.entropy_J_per_kgK, S_C) <= TOL
    r, _, _ = _run(chamber, 2000.0)
    assert _rel(r.p_t, P_T) <= TOL


def test_c02_g2_eps_2000_exit_and_the_undeclared_condensed_species(chamber):
    r, T_e, _ = _run(chamber, 2000.0)
    assert _rel(r.p_e, P_E_2000) <= TOL and _rel(T_e, T_E_2000) <= TOL
    d = r.undeclared_condensed
    assert d is not None and "C(gr)" in d["species"] and d["T_K"] < T_FIRST_NONE


def test_c03_g2_eps_1000_exit_temperature_and_no_undeclared_condensed_species(chamber):
    r, T_e, _ = _run(chamber, 1000.0)
    assert _rel(T_e, T_E_1000) <= TOL
    assert r.undeclared_condensed is None
