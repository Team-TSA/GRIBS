# Model assumptions and limitations

GRIBS is a lumped-parameter numerical model. Users must evaluate whether each assumption is acceptable for their application.

## Principal assumptions

1. The chamber is spatially uniform. Pressure waves, acoustic modes, and local gradients are not solved.
2. Chamber properties come from CEA HP equilibrium (or the constant-property `reference_constant` backend) tabulated over a configured pressure range and interpolated as pressure-only functions. The gas phase follows an ideal-gas state relation; non-ideal effects at high pressure are not modeled or bounded.
3. Regression follows a quasi-steady pressure law. Ignition delay, extinction, transient combustion response, and pressure-history effects are not explicitly modeled.
4. Bore and enabled end surfaces regress uniformly with a common burn depth. The outer grain surface is inhibited.
5. The erosive-burning correction is a simplified port-average model and does not resolve axial flow distribution.
6. Nozzle flow is quasi-one-dimensional and nominally isentropic. Losses are represented by aggregate coefficients.
7. The separation treatment is simplified and does not predict hysteresis, side loads, or unsteady asymmetric flow.
8. Blowdown temperature remains tied to the pressure-dependent chamber-property table rather than an independent energy equation.
9. With throat erosion enabled, exit area remains fixed and full nozzle-contour evolution is not modeled.
10. Numerical self-tests check selected identities and continuity properties, not experimental validity.
11. Two-phase model B (homogeneous equilibrium): the conserved chamber mass is the total product mass (gas + condensed); condensed products move with the gas at the same velocity and temperature, and the equation of state acts on the gas phase only.
12. The condensed-phase volume is neglected against the chamber free volume.
13. The nozzle `c*` uses a constant-gamma isentropic relation. Against the CEA equilibrium reference it differs by up to about 2.3% in the tested condensed AP/Al compositions; the comparison gate G4a has not passed (see the README, "Validation status").

14. The composition is not checked for plausibility. GRIBS does not decide whether a composition is combustible, ignites or sustains combustion. The burn rate is the configured law, so the shapes of the pressure and thrust histories are set mainly by the grain geometry and the burn law. With the geometry, burn law and nozzle fixed, the normalized pressure histories of different compositions nearly overlap: the composition changes mainly their level, and their shape only slightly.
15. Condensed-phase scope is classified by the minimum gas mass fraction `Yg_min` of the thermochemistry table (code constants, not linked to validation evidence): `model_status` is classified with 0.99 / 0.65, and 0.80 only selects an additional warning. For `0.65 < Yg_min <= 0.99` the run reports `W_HIGH_CONDENSATE_MODEL_LIMITATION` (unmodelled two-phase physics). For `0.65 < Yg_min <= 0.80` it also reports `W_HIGH_CONDENSATE_G4B_UNVALIDATED`: reproduction of reported or experimental data (G4b) has not been done for this band. For `Yg_min <= 0.65` the run stops (`E_OUT_OF_VALIDATED_RANGE`).

## Property and input limitations

The chamber-property tables have a configured pressure range. With `thermochemistry.outside_range_policy = "clamp"` (used by the supplied configurations), properties outside that range are held at the nearest table limit and the run is not stopped. Clamping improves numerical robustness outside that interval but suppresses pressure dependence there. Extrapolation remains smooth but may be physically unreliable. Combustion coefficients, efficiencies, flow coefficients, and erosion parameters require independent evidence for the intended formulation and operating range.

## Interpretation limits

- MEOP is the largest simulated chamber pressure, not an allowable structural pressure.
- Effective performance metrics are model outputs, not certification values.
- Mass consistency is an internal numerical diagnostic, not evidence that the physical model is correct.
- Default inputs are examples and must not be treated as design recommendations.

## Validation status

No claim of experimental validation is made by this public repository. Physical validity has not been established; the reproducibility and CEA-consistency evidence and its limits are described in the README, "Validation status". Users should document solver convergence, uncertainty, independent comparisons, and experimental validation before relying on results.
