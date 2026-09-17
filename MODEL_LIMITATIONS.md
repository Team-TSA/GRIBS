# Model assumptions and limitations

GRIBS is a lumped-parameter numerical model. Users must evaluate whether each assumption is acceptable for their application.

## Principal assumptions

1. The chamber is spatially uniform. Pressure waves, acoustic modes, and local gradients are not solved.
2. Chamber gas follows an ideal-gas state relation with fitted `R`, `T0`, and `gamma` treated as pressure-only functions.
3. Regression follows a quasi-steady pressure law. Ignition delay, extinction, transient combustion response, and pressure-history effects are not explicitly modeled.
4. Bore and enabled end surfaces regress uniformly with a common burn depth. The outer grain surface is inhibited.
5. The erosive-burning correction is a simplified port-average model and does not resolve axial flow distribution.
6. Nozzle flow is quasi-one-dimensional and nominally isentropic. Losses are represented by aggregate coefficients.
7. The separation treatment is simplified and does not predict hysteresis, side loads, or unsteady asymmetric flow.
8. Blowdown temperature remains tied to the pressure-property fit rather than an independent energy equation.
9. With throat erosion enabled, exit area remains fixed and full nozzle-contour evolution is not modeled.
10. Numerical self-tests check selected identities and continuity properties, not experimental validity.

## Property and input limitations

The default property fits have a configured pressure range. Clamping improves numerical robustness outside that interval but suppresses pressure dependence there. Extrapolation remains smooth but may be physically unreliable. Combustion coefficients, efficiencies, flow coefficients, and erosion parameters require independent evidence for the intended formulation and operating range.

## Interpretation limits

- MEOP is the largest simulated chamber pressure, not an allowable structural pressure.
- Effective performance metrics are model outputs, not certification values.
- Mass consistency is an internal numerical diagnostic, not evidence that the physical model is correct.
- Default inputs are examples and must not be treated as design recommendations.

## Validation status

No claim of experimental validation is made by this public repository. Users should document solver convergence, uncertainty, independent comparisons, and experimental validation before relying on results.
