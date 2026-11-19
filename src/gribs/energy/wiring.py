"""The wiring of an energy-mode configuration to its equilibrium provider (A-6 stage 2, C11; contract 164 sections 2b-FIX-2, 3-A).

build_energy_inputs(c) turns the configuration into what the integrator core needs: the explicit formulas
(energy_model.propellant_formulas, required) and the reactants of the run give the element mass fractions (energy/elements.py);
the t = 0 element checks of the initial fill run on them (energy/initial.py: the basis covers every element, the weights are
unique and non-negative, the fill is the products of the main propellant); the provider is the one of the configured
thermochemistry backend:

* the official CEA package (backend cea_python): CEATPProvider, with the TP cache of the run's cache directory when the cache is
  enabled; gamma_s is EqSolution.gamma_s of the TP equilibrium, the definition of the legacy gamma_s (an equilibrium isentropic
  exponent), and a value that is not finite or not above 1 stops the run;
* the analytic backend reference_constant: ReferenceConstantProvider, an ideal gas with the configured gamma and molecular weight.
  The legacy backend has no enthalpy, so the convention here is stated in the provenance: the enthalpy at the reference
  temperature is that of the solid at its far temperature (h_s), i.e. the adiabatic flame temperature is the reference
  temperature, so the steady state of energy mode can be compared with the legacy model (algebraic limit, contract 2b-4.8).

The provenance returned with the inputs records every choice the run made (provider, gamma_s definition, formulas, element
fractions, atomic-mass table status, enthalpy basis, the source of rho_p, the recovery range).
"""
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Mapping, Tuple

import gribs.constants
from gribs.energy.balance import ABSOLUTE_ENTHALPY_BASIS, h_s_for_burn_surface
from gribs.energy.elements import element_mass_fractions, reactant_element_rows
from gribs.energy.initial import ElementBasis, check_initial_fill
from gribs.energy.provider import EquilibriumProvider, TPState
from gribs.energy.exceptions import StateRecoveryError
from gribs.errors import ConfigurationError, ThermochemistryError

GAMMA_DEFINITION_CEA = ("EqSolution.gamma_s of the CEA TP equilibrium: the equilibrium isentropic exponent of the mixture, the "
                        "definition of the legacy gamma_s; must be finite and > 1")
GAMMA_DEFINITION_REFERENCE = "the constant gamma of the reference_constant backend (reference_gamma); must be > 1"


_ELEMENT_SUM_TOL = 1.0e-9


@dataclass(frozen=True)
class EnergyRunInputs:
    """What an energy-mode run needs besides the configuration: the provider, the element mass fractions of the propellant and
    of the initial fill (contract U1), and the provenance of the choices that built them.  The searched temperature range of
    the state recovery and its start temperatures are settings (energy_model.recovery); gamma_s comes from the provider's
    state (a state without a finite gamma_s > 1 stops the run)."""
    provider: EquilibriumProvider
    element_fractions: Mapping[str, float]
    provenance: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self):
        fractions = dict(self.element_fractions)
        if not fractions or any(not (math.isfinite(v) and v >= 0.0) for v in fractions.values()):
            raise ValueError("element_fractions must be a non-empty table of finite, non-negative fractions")
        if abs(math.fsum(fractions.values()) - 1.0) > _ELEMENT_SUM_TOL:
            raise ValueError(f"element_fractions must add up to 1, got {math.fsum(fractions.values())!r}")


class ReferenceConstantProvider:
    """An ideal gas of constant gamma and molecular weight, one gas phase and no condensed phase: u(T) = h_ref - R T_ref +
    cv (T - T_ref) with R = Ru / M, cv = R / (gamma - 1); h(T_ref) = h_ref (the enthalpy convention is the caller's)."""
    enthalpy_convention = ABSOLUTE_ENTHALPY_BASIS

    def __init__(self, gamma: float, molecular_weight_kg_kmol: float, reference_temperature_K: float,
                 reference_enthalpy_J_kg: float):
        for name, value in (("gamma", gamma), ("molecular_weight_kg_kmol", molecular_weight_kg_kmol),
                            ("reference_temperature_K", reference_temperature_K),
                            ("reference_enthalpy_J_kg", reference_enthalpy_J_kg)):
            if not math.isfinite(value):
                raise ConfigurationError(f"reference_constant provider: {name} must be finite, got {value!r}")
        if not gamma > 1.0:
            raise ConfigurationError(f"reference_constant provider: gamma must be > 1, got {gamma!r}")
        if not (molecular_weight_kg_kmol > 0.0 and reference_temperature_K > 0.0):
            raise ConfigurationError("reference_constant provider: the molecular weight and the temperature must be positive")
        self.gamma, self.M, self.T_ref, self.h_ref = gamma, molecular_weight_kg_kmol, reference_temperature_K, reference_enthalpy_J_kg
        self.R = gribs.constants.R_UNIVERSAL / molecular_weight_kg_kmol
        self.cv = self.R / (gamma - 1.0)

    def tp(self, T: float, p: float, B_e) -> TPState:
        u = self.h_ref - self.R * self.T_ref + self.cv * (T - self.T_ref)
        return TPState(u_J_per_kg=u, n_gas_per_kg=1.0 / self.M, Yg=1.0, condensed_mass_fraction={}, converged=True,
                       gamma_s=self.gamma)


def _reactants(c) -> Tuple[Tuple[str, ...], Dict[str, float]]:
    specs = getattr(c, "reactants_spec_list", None)
    if not specs:
        raise ConfigurationError("energy mode needs the reactants of the propellant (propellant.reactants)")
    names = tuple(str(r.name).strip() for r in specs)
    return names, {str(r.name).strip(): float(r.wt_percent) for r in specs}


def build_energy_inputs(c):
    """The EnergyRunInputs of an energy-mode configuration (see the module docstring).  Imported lazily by the integrator."""
    e = c.energy
    if e.propellant_formulas is None:
        raise ConfigurationError("energy mode needs energy_model.propellant_formulas (reactant name -> {element: atoms}); "
                                 "formulas are never inferred from the CEA species names")
    names, wt = _reactants(c)
    fractions = element_mass_fractions(wt, e.propellant_formulas)
    basis = ElementBasis.from_mapping(reactant_element_rows(e.propellant_formulas))
    weights = check_initial_fill(fractions, fractions, basis, element_tolerance=e.initial_fill.element_tolerance)
    for name in names:
        if abs(weights[name] - wt[name] / 100.0) > max(e.initial_fill.element_tolerance, 1.0e-9):
            raise ConfigurationError(f"the weights of the reactant basis ({weights}) do not reproduce the configured weight "
                                     f"percent ({wt}); check energy_model.propellant_formulas")
    provenance = {"element_mass_fractions": dict(fractions), "propellant_formulas": {k: dict(v) for k, v in
                  e.propellant_formulas.items()}, "reactants_wt_percent": dict(wt),
                  "atomic_mass_table": "energy/elements.py ATOMIC_MASS_KG_PER_KMOL (IUPAC standard atomic weights, conventional "
                                       "digits; NOT audited against the atomic masses of CEA, held item H-17)",
                  "enthalpy_convention": ABSOLUTE_ENTHALPY_BASIS,
                  "solid_density_source": "Config.rho_p (the one density of the solid, B-15; there is no energy_model key)",
                  "temperature_search_K": list(e.recovery_search.temperature_search_K),
                  "start_temperatures_K": list(e.recovery_search.start_temperatures_K), "backend": c.thermo_backend}
    if c.thermo_backend == gribs.constants.BACKEND_REFERENCE_CONSTANT:
        h_ref = h_s_for_burn_surface(e.burn_surface, products_enthalpy_basis=ABSOLUTE_ENTHALPY_BASIS)
        provider = ReferenceConstantProvider(c.reference_gamma, c.reference_molecular_weight_kg_kmol,
                                             c.reference_temperature_K, h_ref)
        provenance.update(provider="ReferenceConstantProvider", gamma_s_definition=GAMMA_DEFINITION_REFERENCE,
                          reference_enthalpy_J_kg=h_ref,
                          reference_enthalpy_convention="h(T_ref) = h_s(T_grain_inf): the adiabatic flame temperature is "
                                                        "T_ref (an assumption of the algebraic limit; legacy has no enthalpy)")
    elif c.thermo_backend == gribs.constants.BACKEND_CEA_PYTHON:
        from gribs.energy.cea_tp import CEATPProvider, TPCache
        cache = TPCache(Path(c.cea_cache_dir).expanduser()) if (c.cea_cache_enabled and c.cea_cache_dir) else None
        try:
            provider = CEATPProvider(names, [wt[n] for n in names], fractions,
                                     products_from_reactants=c.cea_py_products_from_reactants,
                                     product_species=c.cea_py_product_species, omit_species=c.cea_py_omit_species,
                                     insert_species=c.cea_py_insert_species, ions=c.cea_py_ions, trace=c.cea_py_trace,
                                     smooth_truncation=c.cea_py_smooth_truncation,
                                     truncation_width=c.cea_py_truncation_width,
                                     element_tolerance=e.initial_fill.element_tolerance, cache=cache)
        except StateRecoveryError as exc:           # the provider could not be constructed: a backend initialisation failure
            raise ThermochemistryError(str(exc)) from exc
        provenance.update(provider="CEATPProvider", gamma_s_definition=GAMMA_DEFINITION_CEA, cea_cache=cache is not None,
                          cea_identity=provider.cache_identity())
    else:
        raise ConfigurationError(f"energy mode has no provider for the thermochemistry backend {c.thermo_backend!r}")
    return EnergyRunInputs(provider=provider, element_fractions=fractions, provenance=provenance)
