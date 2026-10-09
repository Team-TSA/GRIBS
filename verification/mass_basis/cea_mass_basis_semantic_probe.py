"""Verify the mass basis of NASA CEA 3.3.4 EqSolution.M.

This probe reproduces RP-1311 Example 14, a TP equilibrium problem with
condensed species. It independently reconstructs the gas-phase molecular
weight and checks the following identities:

    EqSolution.M = 1 / n_gas
    M_CEA = M_gas / Yg
    Ru / M_CEA = Yg * Ru / M_gas
    p / (rho * T) = Ru / M_CEA

The probe also demonstrates that the former GRIBS expression

    Psi_old = Yg * (Ru / M_CEA) * T

applies Yg twice, because Ru / M_CEA is already on a total-mixture mass basis.
"""

from __future__ import annotations

import sys

import cea
import numpy as np


RU_J_PER_KMOL_K = float(cea.R)
ALL_POINT_TOLERANCE = 1.0e-8
CONDENSED_IDENTITY_TOLERANCE = 1.0e-10
OLD_RATIO_TOLERANCE = 1.0e-12


def relative_error(a: float, b: float) -> float:
    scale = max(abs(a), abs(b), 1.0e-300)
    return abs(a - b) / scale


def main() -> int:
    reactant_names = ["H2(L)", "O2(L)"]
    temperatures_K = np.array(
        [1000.0, 500.0, 350.0, 305.0, 304.3, 304.2, 304.0, 300.0],
        dtype=float,
    )
    pressure_bar = float(cea.units.atm_to_bar(0.05))
    pressure_Pa = pressure_bar * 1.0e5

    fuel_moles = np.array([100.0, 0.0])
    oxidant_moles = np.array([0.0, 60.0])

    reactants = cea.Mixture(reactant_names)
    products = cea.Mixture(
        reactant_names,
        products_from_reactants=True,
    )
    solver = cea.EqSolver(products, reactants=reactants)
    solution = cea.EqSolution(solver)

    weights = (
        reactants.moles_to_weights(fuel_moles)
        + reactants.moles_to_weights(oxidant_moles)
    )

    species_names = list(products.species_names)
    num_gas = int(solver.num_gas)
    num_condensed = int(solver.num_condensed)

    gas_names = species_names[:num_gas]
    condensed_names = species_names[
        num_gas:num_gas + num_condensed
    ]

    print("GRIBS CEA mass-basis semantic probe")
    print("===================================")
    print(f"cea_library_version={cea.lib_version()}")
    print(f"Ru_J_per_kmol_K={RU_J_PER_KMOL_K:.15g}")
    print("reference_case=NASA_CEA_RP-1311_Example_14")
    print("constraint=TP")
    print(f"reactants={reactant_names!r}")
    print(f"pressure_bar={pressure_bar:.15g}")
    print(f"pressure_Pa={pressure_Pa:.15g}")
    print(f"num_products={len(species_names)}")
    print(f"num_gas={num_gas}")
    print(f"num_condensed={num_condensed}")
    print(
        "species_ordering="
        "products.species_names[:num_gas] are gas species"
    )
    print()

    rows: list[dict[str, float]] = []

    print(
        "T_K,Yg,Yc,n_kmol_per_kg_total,M_CEA_kg_per_kmol,"
        "M_gas_kg_per_kmol,R_eff_J_per_kgK,p_over_rho_T,"
        "M_identity_relerr,R_identity_relerr,density_EOS_relerr,"
        "old_Psi_over_correct_Psi"
    )

    for temperature_K in temperatures_K:
        solver.solve(
            solution,
            cea.TP,
            float(temperature_K),
            pressure_bar,
            weights,
        )

        if not solution.converged:
            print(f"ERROR: solution did not converge at T={temperature_K}")
            return 1

        nj = np.asarray(solution.nj, dtype=float)
        mass_fractions = {
            name: float(value)
            for name, value in solution.mass_fractions.items()
        }

        n_reported = float(solution.n)
        M_cea = float(solution.M)
        density = float(solution.density)

        n_gas_from_nj = float(np.sum(nj[:num_gas]))
        Yg = float(
            sum(mass_fractions.get(name, 0.0) for name in gas_names)
        )
        Yc = float(
            sum(
                mass_fractions.get(name, 0.0)
                for name in condensed_names
            )
        )

        M_gas = Yg / n_gas_from_nj
        M_from_Mgas_Yg = M_gas / Yg

        R_eff = RU_J_PER_KMOL_K / M_cea
        R_eff_from_gas = Yg * RU_J_PER_KMOL_K / M_gas
        R_from_density = (
            pressure_Pa / (density * float(temperature_K))
        )

        M_identity_error = relative_error(
            M_cea,
            M_from_Mgas_Yg,
        )
        R_identity_error = relative_error(
            R_eff,
            R_eff_from_gas,
        )
        density_EOS_error = relative_error(
            R_eff,
            R_from_density,
        )
        n_identity_error = relative_error(
            n_reported,
            n_gas_from_nj,
        )

        psi_correct = R_eff * float(temperature_K)
        psi_old = Yg * R_eff * float(temperature_K)
        old_ratio = psi_old / psi_correct

        rows.append(
            {
                "T_K": float(temperature_K),
                "Yg": Yg,
                "Yc": Yc,
                "n_identity_error": n_identity_error,
                "M_identity_error": M_identity_error,
                "R_identity_error": R_identity_error,
                "density_EOS_error": density_EOS_error,
                "old_ratio": old_ratio,
            }
        )

        print(
            f"{temperature_K:.6f},"
            f"{Yg:.15g},"
            f"{Yc:.15g},"
            f"{n_reported:.15g},"
            f"{M_cea:.15g},"
            f"{M_gas:.15g},"
            f"{R_eff:.15g},"
            f"{R_from_density:.15g},"
            f"{M_identity_error:.15g},"
            f"{R_identity_error:.15g},"
            f"{density_EOS_error:.15g},"
            f"{old_ratio:.15g}"
        )

    condensed_rows = [
        row for row in rows
        if row["Yc"] > 1.0e-8
    ]

    all_point_maxima = {
        "n_vs_sum_gas_nj_relerr": max(
            row["n_identity_error"] for row in rows
        ),
        "M_vs_Mgas_over_Yg_relerr": max(
            row["M_identity_error"] for row in rows
        ),
        "Reff_vs_Yg_Rgas_relerr": max(
            row["R_identity_error"] for row in rows
        ),
        "Reff_vs_density_EOS_relerr": max(
            row["density_EOS_error"] for row in rows
        ),
    }

    all_point_pass = all(
        value <= ALL_POINT_TOLERANCE
        for value in all_point_maxima.values()
    )

    condensed_identity_pass = bool(condensed_rows) and all(
        row["M_identity_error"]
        <= CONDENSED_IDENTITY_TOLERANCE
        and row["R_identity_error"]
        <= CONDENSED_IDENTITY_TOLERANCE
        and row["density_EOS_error"]
        <= CONDENSED_IDENTITY_TOLERANCE
        and abs(row["old_ratio"] - row["Yg"])
        <= OLD_RATIO_TOLERANCE
        for row in condensed_rows
    )

    print()
    print("Acceptance summary")
    print("------------------")
    print(f"all_point_tolerance={ALL_POINT_TOLERANCE:.15g}")
    print(
        "condensed_identity_tolerance="
        f"{CONDENSED_IDENTITY_TOLERANCE:.15g}"
    )
    print(f"old_ratio_tolerance={OLD_RATIO_TOLERANCE:.15g}")

    for name, value in all_point_maxima.items():
        print(f"max_{name}={value:.15g}")

    print(f"converged_rows={len(rows)}")
    print(f"condensed_rows={len(condensed_rows)}")
    print(f"all_point_pass={all_point_pass}")
    print(f"condensed_identity_pass={condensed_identity_pass}")

    if condensed_rows:
        strongest = min(condensed_rows, key=lambda row: row["Yg"])
        print(
            "strongest_condensation_T_K="
            f"{strongest['T_K']:.15g}"
        )
        print(
            "strongest_condensation_Yg="
            f"{strongest['Yg']:.15g}"
        )
        print(
            "strongest_condensation_Yc="
            f"{strongest['Yc']:.15g}"
        )
        print(
            "strongest_old_Psi_over_correct_Psi="
            f"{strongest['old_ratio']:.15g}"
        )

    passed = all_point_pass and condensed_identity_pass

    if passed:
        print(
            "RESULT=PASS: EqSolution.M=M_gas/Yg; "
            "Ru/M=Yg*Ru/M_gas; the former GRIBS EOS applies Yg twice."
        )
        return 0

    print(
        "RESULT=FAIL: semantic identities were not established; "
        "do not advance G1."
    )
    return 1


if __name__ == "__main__":
    sys.exit(main())
