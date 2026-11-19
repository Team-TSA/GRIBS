"""The CEA TP equilibrium provider and a TP cache in its own namespace (A-6 stage 2, C7; contract 164 sections 3-A, 5; F-21, F-24).

The provider solves the TP equilibrium of a FIXED reactant composition (the main propellant; stage 2 has a single source and
a fill of the products, contract assumption U1) with the official NASA CEA package, imported and its database pinned through
gribs.thermo exactly as the legacy backend does.  The element inventory B_e passed to tp() must have the element fractions
that the caller gave (element_mass_fractions); a composition that changes is outside stage 2.  The caller's fractions are
NOT audited against the species of CEA here.

Units: the energy-family values of EqSolution are kJ-based and are converted to J in one place (_KJ).  Everything is per kg of
TOTAL products (gas plus condensed): n_gas_per_kg = EqSolution.n = 1/M_CEA, so M_CEA * n = 1 is checked; R_eff = Ru n,
R_gas = R_eff / Yg and M_gas = Yg M_CEA are derived with gribs.energy.relations.

The TP cache holds states at exact (T, p) nodes in files thermo_cea_tp_<hash>.json.  Its key holds the schema, the CEA version
and database digest, the reactants and weights, the element order and fractions, the product-species options and the enthalpy
convention.  The p-only caches of the legacy backend (thermo_cea_python_*.json) are neither read nor written.
"""
from __future__ import annotations

import hashlib
import json
import math
import threading
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Mapping, Optional, Sequence

import gribs.thermo
from gribs.energy import relations
from gribs.energy.balance import ABSOLUTE_ENTHALPY_BASIS
from gribs.energy.exceptions import StateRecoveryError
from gribs.errors import ThermochemistryCacheError
from gribs.energy.provider import TPState

TP_CACHE_SCHEMA = "gribs-energy-tp-cache-2"        # 2: the state carries gamma_s (C11); schema-1 files are not read
TP_CACHE_PREFIX = "thermo_cea_tp_"
_KJ = 1000.0                      # EqSolution energy, enthalpy, cp and cv are kJ-based; J = 1000 * kJ
_INVARIANT_TOL = 1.0e-10


@dataclass(frozen=True)
class CEATPDetail:
    """One validated TP equilibrium, in SI units and per kg of total products."""
    T_K: float
    p_Pa: float
    u_J_per_kg: float
    h_J_per_kg: float
    n_gas_per_kg: float             # kmol of gas per kg of total = 1 / M_CEA
    M_CEA: float                    # kg-total per kmol-gas
    Yg: float
    condensed_mass_fraction: Dict[str, float]
    cea_version: str
    gamma_s: float                  # EqSolution.gamma_s: equilibrium isentropic exponent (the legacy gamma_s), > 1


def _canonical(obj) -> str:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _tp_cache_write_failed(phase: str, step: str, path, original: BaseException, replaced: bool,
                           cleanup: Optional[Mapping]) -> ThermochemistryCacheError:
    """The error of a failed TP cache write (code E_THERMO_CACHE_WRITE_FAILED, exit code 3 on the command line, like the legacy
    chamber-property cache).  The wording is the energy mode's own and says "TP cache"; the legacy wording and code are unchanged."""
    c = dict(cleanup or {})
    if phase == "create_directory":
        state = "the TP cache directory could not be created; no TP cache file was written"
    elif replaced:
        state = ("the TP cache file name was changed by this run, so any previous TP cache file is gone; the new file's "
                 "identity or durability could not be confirmed")
    else:
        state = "this run did not replace the TP cache file (an existing one, if any, was not modified by this run)"
    return ThermochemistryCacheError(
        f"TP cache could not be written ({phase}, {step}): {path}\n"
        f"Original error: {type(original).__name__}: {original}\n{state}; temp file: {c.get('status', 'unknown')}.",
        code="E_THERMO_CACHE_WRITE_FAILED",
        context={"phase": phase, "step": step, "path": str(path), "error": f"{type(original).__name__}: {original}",
                 "errno": getattr(original, "errno", None), "replaced": bool(replaced), "temp_cleanup": c, "cache": "tp"})


class TPCache:
    """Exact-node cache of TP states (see the module docstring); files are written atomically through gribs.thermo."""

    def __init__(self, directory):
        self.directory = Path(directory)
        self.hits = 0
        self.stores = 0

    @staticmethod
    def key(identity: Mapping, T: float, p: float) -> str:
        return hashlib.sha256(_canonical({"schema": TP_CACHE_SCHEMA, "identity": identity, "T": T, "p": p}).encode()).hexdigest()

    def path(self, identity: Mapping, T: float, p: float) -> Path:
        return self.directory / f"{TP_CACHE_PREFIX}{self.key(identity, T, p)[:24]}.json"

    def load(self, identity: Mapping, T: float, p: float) -> Optional[TPState]:
        path = self.path(identity, T, p)
        if not path.is_file():
            return None
        try:
            doc = json.loads(path.read_text(encoding="utf-8"))
            state = doc["state"]
            ok = (doc["schema"] == TP_CACHE_SCHEMA and doc["key"] == self.key(identity, T, p)
                  and doc["state_sha256"] == hashlib.sha256(_canonical(state).encode()).hexdigest())
            restored = TPState(u_J_per_kg=state["u_J_per_kg"], n_gas_per_kg=state["n_gas_per_kg"], Yg=state["Yg"],
                               condensed_mass_fraction=dict(state["condensed_mass_fraction"]), converged=True,
                               gamma_s=state["gamma_s"])
        except (OSError, ValueError, KeyError, TypeError) as exc:
            raise StateRecoveryError(f"the TP cache file {path} is unreadable or malformed: {type(exc).__name__}: {exc}") from exc
        if not ok:
            raise StateRecoveryError(f"the TP cache file {path} does not match its key or its checksum")
        self.hits += 1
        return restored

    def store(self, identity: Mapping, T: float, p: float, state: TPState) -> None:
        body = {"u_J_per_kg": state.u_J_per_kg, "n_gas_per_kg": state.n_gas_per_kg, "Yg": state.Yg,
                "condensed_mass_fraction": dict(state.condensed_mass_fraction), "gamma_s": state.gamma_s}
        doc = {"schema": TP_CACHE_SCHEMA, "key": self.key(identity, T, p), "identity": identity, "T": T, "p": p,
               "state": body, "state_sha256": hashlib.sha256(_canonical(body).encode()).hexdigest()}
        try:                                                       # as the legacy cache: the directory first
            self.directory.mkdir(parents=True, exist_ok=True)
        except OSError as exc:
            raise _tp_cache_write_failed("create_directory", "create_directory", self.directory, exc, False,
                                         {"status": "not_created", "created": False}) from exc
        try:
            gribs.thermo._cache_atomic_write_text(self.path(identity, T, p), _canonical(doc) + "\n")
        except gribs.thermo._CacheWriteError as exc:
            raise _tp_cache_write_failed(exc.phase, exc.step, exc.path, exc.original, exc.replaced, exc.cleanup) from exc
        self.stores += 1


class CEATPProvider:
    #: The enthalpy basis of CEA's values: absolute, formation enthalpy included (what h_s must share, contract 2b-4.3).
    enthalpy_convention = ABSOLUTE_ENTHALPY_BASIS

    def __init__(self, reactant_names: Sequence[str], reactant_wt_percent: Sequence[float],
                 element_mass_fractions: Mapping[str, float], *, products_from_reactants: bool = True,
                 product_species: Sequence[str] = (), omit_species: Sequence[str] = (), insert_species: Sequence[str] = (),
                 ions: bool = False, trace: float = 1.0e-10, smooth_truncation: bool = False,
                 truncation_width: float = -1.0, element_tolerance: float = 1.0e-9, cache: Optional[TPCache] = None):
        names = tuple(str(n).strip() for n in reactant_names)
        weights = tuple(float(w) for w in reactant_wt_percent)
        if not names or any(not n for n in names) or len(set(names)) != len(names) or len(weights) != len(names):
            raise ValueError("reactant_names must be unique, non-empty and as many as reactant_wt_percent")
        if any(not (math.isfinite(w) and w > 0.0) for w in weights) or abs(sum(weights) - 100.0) > 1.0e-8:
            raise ValueError("reactant_wt_percent must be positive, finite and sum to 100")
        if not element_mass_fractions or any(not (math.isfinite(v) and v >= 0.0) for v in element_mass_fractions.values()):
            raise ValueError("element_mass_fractions must be non-negative and finite")
        total = math.fsum(element_mass_fractions.values())
        if not total > 0.0:
            raise ValueError("element_mass_fractions has no mass")
        self._elements = tuple(element_mass_fractions)
        self._element_fractions = tuple(element_mass_fractions[e] / total for e in self._elements)
        self._element_tolerance = float(element_tolerance)
        self._names, self._weights, self._cache = names, weights, cache
        self._options = {"products_from_reactants": bool(products_from_reactants),
                         "product_species": [str(s) for s in product_species], "omit_species": [str(s) for s in omit_species],
                         "insert_species": [str(s) for s in insert_species], "ions": bool(ions), "trace": float(trace),
                         "smooth_truncation": bool(smooth_truncation), "truncation_width": float(truncation_width)}
        cea = gribs.thermo.import_official_cea_module()
        self._cea = cea
        self.database = gribs.thermo.pin_bundled_cea_database(cea, False)
        if abs(float(cea.R) - 8314.51) > 1.0e-9:
            raise StateRecoveryError(f"cea.R = {cea.R!r} differs from the Ru = 8314.51 J/(kmol K) of the contract")
        try:
            reactants = cea.Mixture(list(names))
            if products_from_reactants:
                products = cea.Mixture(list(names), products_from_reactants=True, omit=list(omit_species), ions=bool(ions))
            else:
                products = cea.Mixture(list(product_species), omit=list(omit_species), ions=bool(ions))
            kwargs = {"reactants": reactants, "transport": False, "ions": bool(ions), "trace": float(trace),
                      "insert": list(insert_species)}
            if smooth_truncation:
                kwargs.update(smooth_truncation=True, truncation_width=float(truncation_width))
            self._solver = cea.EqSolver(products, **kwargs)
            self._solution = cea.EqSolution(self._solver)
        except Exception as exc:
            raise StateRecoveryError(f"CEA could not construct the TP equilibrium solver: {type(exc).__name__}: {exc}") from exc
        species = tuple(products.species_names)
        n_gas = int(self._solver.num_gas)
        self._gas_names, self._condensed_names = species[:n_gas], species[n_gas:]
        self._lock = threading.Lock()

    def cache_identity(self) -> Dict[str, object]:
        """What the key of the TP cache depends on (besides T and p)."""
        return {"cea_version": str(self._cea.__version__), "database_thermo_lib_sha256": self.database["thermo_lib_sha256"],
                "reactant_names": list(self._names), "reactant_wt_percent": list(self._weights),
                "element_order": list(self._elements), "element_mass_fractions": list(self._element_fractions),
                "options": dict(self._options), "enthalpy_convention": self.enthalpy_convention}

    def _check_composition(self, B_e: Mapping[str, float]) -> float:
        mass = math.fsum(B_e.values())
        if not (math.isfinite(mass) and mass > 0.0):
            raise StateRecoveryError(f"the element inventory must have a positive mass, got {mass!r} kg")
        given = {e: B_e.get(e, 0.0) / mass for e in set(B_e) | set(self._elements)}
        expected = dict(zip(self._elements, self._element_fractions))
        worst = max(abs(given.get(e, 0.0) - expected.get(e, 0.0)) for e in given)
        if worst > self._element_tolerance:
            raise StateRecoveryError(
                f"the element fractions of B_e differ from those of the provider's composition by {worst:.3g} (> "
                f"{self._element_tolerance:g}); a composition that changes is outside stage 2 (single source, fill of the products)")
        return mass

    def tp_detail(self, T: float, p: float, B_e: Mapping[str, float]) -> CEATPDetail:
        if not (math.isfinite(T) and T > 0.0 and math.isfinite(p) and p > 0.0):
            raise StateRecoveryError(f"T and p must be positive and finite, got T = {T!r}, p = {p!r}")
        self._check_composition(B_e)
        with self._lock:
            try:
                import numpy as np
                self._solver.solve(self._solution, self._cea.TP, float(T), float(p) / 1.0e5, np.asarray(self._weights, dtype=float))
            except Exception as exc:
                raise StateRecoveryError(f"CEA TP solve raised {type(exc).__name__}: {exc}") from exc
            sol = self._solution
            if not bool(sol.converged):
                raise StateRecoveryError(f"the CEA TP equilibrium did not converge at T = {T:.9g} K, p = {p:.9g} Pa "
                                         f"(last_error = {int(sol.last_error)})")
            mass_fractions = {str(k): float(v) for k, v in dict(sol.mass_fractions).items()}
            n, M = float(sol.n), float(sol.M)
            energy, enthalpy = float(sol.energy) * _KJ, float(sol.enthalpy) * _KJ
            gamma = float(sol.gamma_s)
        if not math.isfinite(gamma) or not gamma > 1.0:
            raise StateRecoveryError(f"CEA returned an isentropic exponent gamma_s = {gamma!r} (must be finite and > 1) at "
                                     f"T = {T:.9g} K, p = {p:.9g} Pa")
        values = [n, M, energy, enthalpy, *mass_fractions.values()]
        if not all(math.isfinite(v) for v in values) or not (n > 0.0 and M > 0.0):
            raise StateRecoveryError(f"CEA returned a non-finite or non-positive value at T = {T:.9g} K, p = {p:.9g} Pa")
        if abs(n * M - 1.0) > _INVARIANT_TOL:
            raise StateRecoveryError(f"CEA M and n do not satisfy M_CEA = 1/n_gas: M = {M!r}, n = {n!r}")
        if abs(math.fsum(mass_fractions.values()) - 1.0) > _INVARIANT_TOL:
            raise StateRecoveryError(f"the CEA mass fractions do not sum to 1: {math.fsum(mass_fractions.values())!r}")
        condensed = {k: mass_fractions[k] for k in self._condensed_names if mass_fractions.get(k, 0.0) > 0.0}
        yg = math.fsum(mass_fractions[k] for k in self._gas_names)
        return CEATPDetail(T_K=float(T), p_Pa=float(p), u_J_per_kg=energy, h_J_per_kg=enthalpy, n_gas_per_kg=n, M_CEA=M,
                           Yg=yg, condensed_mass_fraction=condensed, cea_version=str(self._cea.__version__), gamma_s=gamma)

    def tp(self, T: float, p: float, B_e: Mapping[str, float]) -> TPState:
        if self._cache is not None:
            self._check_composition(B_e)
            hit = self._cache.load(self.cache_identity(), float(T), float(p))
            if hit is not None:
                return hit
        d = self.tp_detail(T, p, B_e)
        state = TPState(u_J_per_kg=d.u_J_per_kg, n_gas_per_kg=d.n_gas_per_kg, Yg=d.Yg,
                        condensed_mass_fraction=dict(d.condensed_mass_fraction), converged=True, gamma_s=d.gamma_s)
        if self._cache is not None:
            self._cache.store(self.cache_identity(), float(T), float(p), state)
        return state

    def derived(self, T: float, p: float, B_e: Mapping[str, float]) -> Dict[str, float]:
        """n_gas [kmol], m_total and m_gas [kg], Yg, M_CEA and M_gas [kg/kmol], R_eff and R_gas [J/(kg K)] (F-21: two
        different masses and two different gas constants, never one name for both)."""
        d = self.tp_detail(T, p, B_e)
        mass = math.fsum(B_e.values())
        r_eff = float(self._cea.R) * d.n_gas_per_kg
        return {"n_gas_kmol": mass * d.n_gas_per_kg, "m_total_kg": mass, "m_gas_kg": mass * d.Yg, "Yg": d.Yg,
                "M_CEA": d.M_CEA, "M_gas": relations.m_gas(d.M_CEA, d.Yg), "R_eff": r_eff, "R_gas": relations.r_gas(r_eff, d.Yg)}
