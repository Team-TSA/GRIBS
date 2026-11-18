# Model assumptions and limitations

GRIBS is a lumped-parameter numerical model. Users must evaluate whether each assumption is acceptable for their application.

## Principal assumptions

1. The chamber is spatially uniform. Pressure waves, acoustic modes, and local gradients are not solved.
2. Chamber properties come from CEA HP equilibrium (or the constant-property `reference_constant` backend) tabulated over a configured pressure range and interpolated as pressure-only functions. The gas phase follows an ideal-gas state relation; non-ideal effects at high pressure are not modeled or bounded.
3. Regression follows a quasi-steady pressure law. Ignition delay, extinction, transient combustion response, and pressure-history effects are not explicitly modeled.
4. Bore and enabled end surfaces regress uniformly with a common burn depth. The outer grain surface is inhibited.
5. The erosive-burning correction does not resolve the axial flow distribution. It uses one mass flux, the total product mass flux at the aft end of the bore, G = rho_p·r·Ab/(π·Ri²) (it includes the products of the burning end faces and the condensed phase, not the igniter flow), and applies the resulting burn rate r = r0 + alpha·G^0.8/(2Ri)^0.2·exp(−beta·rho_p·r/G) uniformly to the whole burning surface, end faces included. It is therefore not a port-average model. The form is that of Lenoir–Robert with the bore diameter 2Ri as the length scale; alpha is in m^2.8·kg^−0.8·s^−0.2 (SI) and beta is dimensionless. Because G is proportional to r, the equation is solved for r at each evaluation (Newton's method); if it does not converge, the run stops with `E_EROSIVE_ROOT_NOT_CONVERGED` (exit code 4) instead of using an unconverged value.
6. Nozzle flow is quasi-one-dimensional and nominally isentropic. Losses are represented by aggregate coefficients.
7. The separation treatment is simplified and does not predict hysteresis, side loads, or unsteady asymmetric flow.
8. Blowdown temperature remains tied to the pressure-dependent chamber-property table rather than an independent energy equation.
9. With throat erosion enabled, exit area remains fixed and full nozzle-contour evolution is not modeled. If the eroding throat area becomes larger than the fixed exit area (At > Ae, compared on the radius as Rt > Rt0·√ε), the nozzle is computed with the exit as its smallest section, which is outside the model domain. The run then records the first such time in `event_times.throat_exceeds_exit_s` (an empty list when it does not occur) and reports the warning `W_THROAT_EXCEEDS_EXIT_AREA`; it continues, and the warning does not change the run verdict or the exit code.
10. Numerical self-tests check selected identities and continuity properties, not experimental validity.
11. Two-phase model B (homogeneous equilibrium): the conserved chamber mass is the total product mass (gas + condensed); condensed products move with the gas at the same velocity and temperature, and the equation of state acts on the gas phase only.
12. The condensed-phase volume is neglected against the chamber free volume.
13. The nozzle `c*` uses a constant-gamma isentropic relation. Against the CEA equilibrium reference it differs by up to about 2.3% in the tested condensed AP/Al compositions; the comparison gate G4a has not passed (see the README, "Validation status").

14. The composition is not checked for plausibility. GRIBS does not decide whether a composition is combustible, ignites or sustains combustion. The burn rate is the configured law, so the shapes of the pressure and thrust histories are set mainly by the grain geometry and the burn law. With the geometry, burn law and nozzle fixed, the normalized pressure histories of different compositions nearly overlap: the composition changes mainly their level, and their shape only slightly.
15. Condensed-phase scope is classified by the minimum gas mass fraction `Yg_min` of the thermochemistry table (code constants, not linked to validation evidence): `model_status` is classified with 0.99 / 0.65, and 0.80 only selects an additional warning. For `0.65 < Yg_min <= 0.99` the run reports `W_HIGH_CONDENSATE_MODEL_LIMITATION` (unmodelled two-phase physics). For `0.65 < Yg_min <= 0.80` it also reports `W_HIGH_CONDENSATE_G4B_UNVALIDATED`: reproduction of reported or experimental data (G4b) has not been done for this band. For `Yg_min <= 0.65` the run stops (`E_OUT_OF_VALIDATED_RANGE`).

## Property and input limitations

The chamber-property tables have a configured pressure range. With `thermochemistry.outside_range_policy = "clamp"` (used by the supplied configurations), properties outside that range are held at the nearest table limit and the run is not stopped. A run that crosses the table limit reports `W_PROPERTY_TABLE_RANGE_CROSSED` in the summary warnings. A run that spends part of its history outside the table reports `W_PROPERTY_TABLE_EXTRAPOLATED` with the fraction of sampled points outside. Golden case g2 (table 0.1-8 MPa, MEOP about 44.5 MPa) spends about 84 % of its sampled history above the table limit: it is a clamped-table regression reference, not a physically representative case, and its frozen values are intentionally unchanged. Clamping improves numerical robustness outside that interval but suppresses pressure dependence there. A run that leaves the range reports the warning `W_PROPERTY_TABLE_EXTRAPOLATED` (number and fraction of samples, policy). Extrapolation remains smooth but may be physically unreliable. Combustion coefficients, efficiencies, flow coefficients, and erosion parameters require independent evidence for the intended formulation and operating range.

## Run state

The run summary reports six run-state keys: `solver_completed`, `mass_consistency_passed`, `nominal_termination_reached`, `model_domain_passed`, `warnings_present` and `results_usable`. `results_usable` is true only when the first four are true. They are reported only for a run that completed; a run that stops with an error before the summary is made (configuration, thermochemistry, model-applicability or solver error) has no summary and no keys. They do not change `run_passed` or the exit code.

- A run that ends other than by burnout (and, with blowdown enabled, by reaching ambient pressure) reports `W_NON_NOMINAL_TERMINATION` and is not usable.
- At > Ae (`W_THROAT_EXCEEDS_EXIT_AREA`, item 9) and `E_OUT_OF_VALIDATED_RANGE` are model-domain violations; the run is not usable.
- Leaving the property-table pressure range (`W_PROPERTY_TABLE_EXTRAPOLATED`, see below) and subsonic nozzle flow during burning (`W_SUBSONIC_BURNING`; burn-phase c* indicators lose their meaning there) are warnings only and do not make the run unusable. Golden g1 and g2 leave the table range and are reported usable with this warning.
- `results_usable` is a statement about these checks only. It is not evidence that the physical model is valid (see "Validation status").

## Interpretation limits

- MEOP is the largest simulated chamber pressure, not an allowable structural pressure.
- Effective performance metrics are model outputs, not certification values.
- Mass consistency is an internal numerical diagnostic, not evidence that the physical model is correct.
- Default inputs are examples and must not be treated as design recommendations.

## Validation status

No claim of experimental validation is made by this public repository. Physical validity has not been established. Reproduction of reported or experimental data (gate G4b) has not been done, because no usable public data was found: the reported AP/Al cases do not give complete inputs, and no experimental data set or independent two-phase model was available. The reproducibility and CEA-consistency evidence and its limits are described in the README, "Validation status". Users should document solver convergence, uncertainty, independent comparisons, and experimental validation before relying on results.
