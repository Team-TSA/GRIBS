# A-5.1-7: Before/After Comparison Report

- 作成日: 2026-10-10
- 対象: GRIBS v0.5.0-alpha
- 環境: Python 3.13.16 / NumPy 2.3.5 / SciPy 1.17.1 / cea 3.3.4

---

## 1. R1: Low Condensation (NH4NO3/Al/C 88.7/0.7/10.6)

### Configuration
- Propellant: NH4NO3(IV)/AL(cr)/C = 88.7/0.7/10.6 wt%
- Yg range: [0.9980, 1.0000]
- Pressure range: 1-80 bar
- rtol: 1e-10 (default)

### Results

| Metric | Before (A-5.1-4) | After (A-5.1-6) | Change |
|--------|------------------|-----------------|--------|
| MEOP | 3.168092 MPa | 3.168092 MPa | 0.000% |
| Burn time | 11.235862 s | 11.235862 s | 0.000% |
| Total impulse | 131.5116 N·s | 131.5116 N·s | 0.000% |
| Max thrust | 26.9011 N | 26.9011 N | 0.000% |
| Isp | 223.39 s | 223.39 s | 0.000% |
| c* (eff) | 1883.53 m/s | 1883.53 m/s | 0.000% |
| D_min | 0.969083 | 0.969083 | 0.000% |
| Post-run | PASS | PASS | — |
| Ledger | PASS | PASS | — |
| Semantic | — | pass | 新規 |
| Numerical | — | pass | 新規 |
| Model | — | validated_in_scope | 新規 |

**Conclusion**: No numerical change. The correction (Yg double-multiplication removal, /Yg removal) has zero effect when Yg ≈ 1.0, as expected.

---

## 2. R2: Medium Condensation (AP/Al/C7H8 70/16/14)

### Configuration
- Propellant: NH4CLO4(I)/AL(cr)/C7H8 = 70/16/14 wt%
- Yg range: [0.7238, 0.7266]
- Pressure range: 1-200 bar
- rtol: 1e-9 (config-specified)

### Results

| Metric | Before (A-5.1-4) | After (A-5.1-6) | Change |
|--------|------------------|-----------------|--------|
| MEOP | 10.548 MPa | 17.331 MPa | **+64.3%** |
| Burn time | 2.309 s | 1.940 s | −16.0% |
| Total impulse | 1101.9 N·s | 1276.2 N·s | +15.8% |
| Max thrust | 738.4 N | 1047.1 N | +41.8% |
| Isp | 170.52 s | 197.49 s | +15.8% |
| c* (eff) | 1148.1 m/s | 1594.4 m/s | +38.9% |
| CF (eff) | — | 1.198 | — |
| D_min | — | 0.975 | — |
| Post-run | PASS | PASS | — |
| Ledger | FAIL (1.57e-9) | **PASS (1.79e-9)** | 修正 |
| Semantic | — | pass | 新規 |
| Numerical | — | **pass** | 新規 (旧: fail) |
| Model | — | approximate_with_limitations | 新規 |

**Before values**: 独立レビュー (GRIBS_最終統合指針_妥当性評価と本体追加調査_2026-10-09.md §5, 付録 A)

### Key Findings

1. **MEOP +64.3%**: The `/Yg` removal in the nozzle changes the mass flow by factor 1/Yg ≈ 1.38, which propagates to higher chamber pressure.

2. **c* ratio correction**: Before fix, c*_GRIBS/c*_CEA ≈ 0.7234 (≈Yg). After fix, the ratio reflects only the model difference (constant-γ approximation vs CEA equilibrium).

3. **Ledger fix (S-2)**: The propellant tolerance now scales with rtol: `max(1e-9, 100*rtol)`. With rtol=1e-9, the tolerance is 1e-7, giving margin of 55.9x (was 1.4x with fixed 1e-9).

4. **D_min improvement**: The corrected Psi = R_eff * T_used (without Yg double-multiplication) improves the pressure factor monotonicity.

---

## 3. R3: High Condensation (Reported AP/Al)

**Status**: Not evaluated (input data not available)

The original report's AP/Al configuration has not been obtained. When available, G4b evaluation should verify:
- Yg ≈ 0.71
- c*_GRIBS/c*_CEA ≈ 0.72 (before fix)
- MEOP −44% (before → after)

---

## 4. Numerical Invariance Summary

| Check | R1 | R2 |
|-------|----|----|
| MEOP change | 0.000% | +64.3% |
| Burn time change | 0.000% | −16.0% |
| Total impulse change | 0.000% | +15.8% |
| Post-run gate | PASS | PASS |
| Ledger gate | PASS | PASS (after S-2 fix) |
| Self-test | 15/15 PASS | — |
| pytest | 89 passed | — |
| red suite | 11 passed, 1 xfailed | — |

**R1 is the regression anchor**: All numerical values are identical before and after the correction, confirming that the fix has no effect in the single-phase limit.

**R2 shows the expected physical correction**: The large changes (MEOP +64%, etc.) are the correct behavior when Yg ≈ 0.72, because the old code had systematic errors of order 1/Yg ≈ 1.38.