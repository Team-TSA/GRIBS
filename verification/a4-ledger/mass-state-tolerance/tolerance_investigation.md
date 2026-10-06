# A-4 Mass-State Tolerance Investigation

- LSODA accepts vector rtol: `True`
- m_gen relative term / atol: `6.2031184769999996e+01`
- m_igniter relative term / atol: `2.0000000000000004e+00`

## Interpretation

The installed SciPy LSODA interface accepts a state-wise relative tolerance array.
The generated-mass state is controlled primarily by the relative tolerance term. Tightening only its absolute tolerance is unlikely to materially reduce the observed residual.
The igniter-mass state may still be sensitive to its absolute tolerance.
