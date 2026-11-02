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

---

## 5. R4: CEA-backend evidence candidates — low (anchor), medium and high condensation (appended 2026-11-02)

### 5.0 Reading notes (read first)

**§1〜§4 の status 欄について:**
- §1〜§4 の status 欄の値（とくに §1 の Model `validated_in_scope`）は、2026-10-10 の時点でコードが出力した token である。**検証の承認や合格を示すものではない**
- 現在の承認状態（2026-11-02）:
  - A-5.1 の G4a は、c\* の許容差 2e-3 で**未合格**
  - G4b は未評価（入力が無い）
  - G1 は測定済みで**未承認**
  - したがって、どのケースも「検証済み（validated）」ではない
- **model_status の決まり方**（上から順に判定し、最初に当てはまったものを返す）
  - Yg_min は property table の気相質量分率の最小値。「2相が有効」は `two_phase_model.mode = homogeneous_equilibrium` で、かつ backend が2相に対応していること

| 順 | 条件 | cleanup P1（`dc073bf`）より前の commit | `dc073bf` 以降の commit。§5 の測定に使った本体 `d13a9120` も同じ規則 |
|---|---|---|---|
| 1 | `two_phase_model.mode = single_phase_legacy` | `validated_in_scope`（Yg_min に関係なく） | `regression_comparison_only`（Yg_min に関係なく） |
| 2 | 2相が無効（backend が2相に非対応など） | `validated_in_scope` | `in_scope_not_validated` |
| 2a | 2相が有効で、Yg の範囲が無い・数値でない・非有限・0 < Yg_min ≤ Yg_max ≤ 1 を満たさない | `1810e67` より前: 既定値 1.0 として扱われた。`1810e67` 以降: `E_PHASE_FRACTION_VALIDATION_INVALID` で止まり、status は付かない（巨大な整数など、値を文字列にできない場合の報告は `3f21597`・`e124ab4` で修正） | `E_PHASE_FRACTION_VALIDATION_INVALID` で止まり、status は付かない |
| 3 | 2相が有効で Yg_min ≤ 0.65 | `out_of_scope`。`1810e67`（A-9）以降は、積分の前に `E_OUT_OF_VALIDATED_RANGE` で止まる | `out_of_scope`。積分の前に `E_OUT_OF_VALIDATED_RANGE` で止まる |
| 4 | 2相が有効で 0.65 < Yg_min ≤ 0.99 | `approximate_with_limitations` | `approximate_with_limitations` |
| 5 | 2相が有効で Yg_min > 0.99 | `validated_in_scope` | `in_scope_not_validated` |

  - single_phase_legacy は最優先なので、Yg の値や2相が有効かどうかには左右されない
  - 閾値（0.99 / 0.65）と順 3・4 の結果は、両者で同じである。違いは token の名前（順 1・2・5）と、範囲が使えない場合・範囲外の場合に止めるかどうか（順 2a・3）だけ
  - どちらの token も、検証に合格したことは意味しない
  - 下の §5 の status は、`dc073bf` 以降の token で書いた。§5 の8つの run はすべて `homogeneous_equilibrium`、backend は `cea_python`（2相が有効）で、Yg の範囲は有効だったので、順 3〜5 だけが当てはまる

**本節の値の区別:**
- **[出力]** = GRIBS の本体が summary.json / time_history.csv に出力した値を、そのまま書いたもの
- **[算出]** = 出力値から、この節のために計算したもの（変化率、比、修正前の本体が出力しない指標の再構成）

### 5.1 Conditions
- **Before:** 本体 `35b599c`（2026-10-06）。質量基準の修正 `36927b9`（Psi）/ `16b9d9f`（ノズル）の直前の commit。git の履歴から読み取りだけで取り出した
- **After:** 本体 `d13a9120`
  - 測定した時点（2026-11-02）では未 commit の working tree で、A-9〜A-11、cleanup P1〜P5、CLI v5 を含んでいた。**commit 済みのどの版とも違う**
  - その一部（S0 と修正 A＋B、P1）は、のちに `1810e67`〜`dc073bf` として commit された
  - **それらの commit の本体で、§5 の値は再測定していない**
- 同じ設定ファイルを、両方の本体がそのまま受理した。変えたのは出力先と cache の場所だけ
- 環境: Python 3.14 / cea 3.3.4
- 調査の記録（repo の外）: `GRIBS-cleanup-2026-11-01/runs/r160_a51_before_after`、`runs/r167_a51_before_after_columns`

| ケース | 組成（wt%） | Yg の範囲（table）[出力] | 役割 |
|---|---|---|---|
| g3 | golden の g3 | [0.99796, 1.0] | 低凝縮の基準（anchor） |
| Al5 | AP/Al/HTPB 70/5/25 | [0.90602, 0.90668] | 中凝縮（AL2O3 は固相のまま） |
| Al7 | AP/Al/HTPB 70/7/23 | [0.86918, 0.86977] | 中凝縮（AL2O3 は燃焼室で液相、スロートで固相） |
| g4 | golden の g4 | [0.78450, 0.78797] | 名前は medium だが、10/07 の定義では高凝縮 |

- Al5 / Al7 は、g4 から作った**証拠候補の設定で、golden ではない**

### 5.2 Output values（両方の本体が出力する項目）

値は [出力]、括弧内の変化率は [算出]。

| 指標 | g3 前 → 後 | Al5 前 → 後 | Al7 前 → 後 | g4 前 → 後 |
|---|---|---|---|---|
| MEOP | 6.9316 → 6.9316 MPa（−2.4e-10） | 3.9017 → 4.5694 MPa（+17.11%） | 3.9347 → 4.9268 MPa（+25.21%） | 3.9998 → 5.9132 MPa（+47.84%） |
| Burn time | 1.10045 → 1.10045 s | 1.15056 → 1.08338 s（−5.84%） | 1.14685 → 1.05277 s（−8.20%） | 1.14086 → 0.98358 s（−13.79%） |
| Total impulse | 137.780 → 137.780 N·s | 88.927 → 93.119（+4.71%） | 91.374 → 97.546（+6.75%） | 97.154 → 108.845（+12.03%） |
| Max thrust | 238.95 → 238.95 N | 142.88 → 158.55（+10.97%） | 147.25 → 170.73（+15.95%） | 157.68 → 204.04（+29.40%） |
| Isp | 234.57 → 234.57 s | 173.52 → 181.70（+4.71%） | 174.36 → 186.13（+6.75%） | 174.92 → 195.97（+12.03%） |
| c\* (eff, burn interval) | 1901.87 → 1901.87 m/s | 1331.28 → 1468.56（+10.31%） | 1308.77 → 1504.91（+14.99%） | 1246.16 → 1587.70（+27.41%） |
| CF (eff) | 1.2058 → 1.2058 | 1.2737 → 1.2088（−5.09%） | 1.3019 → 1.2083（−7.19%） | 1.3718 → 1.2056（−12.11%） |
| D_min（min(1 − pΨ'/Ψ)） | 0.967670 → 0.969083 | 0.995165 → 0.995076 | 0.991416 → 0.991029 | 0.980910 → 0.980709 |
| CEA の理論 c\*（各 run の平均燃焼圧） | 1909.34 → 1909.34 | 1470.81 → 1470.80 | 1544.42 → 1544.66 | 1588.75 → 1590.59 |
| CEA の比較圧力（同上） | 3.662 → 3.662 MPa | 2.139 → 2.506 | 2.157 → 2.702 | 2.188 → 3.234 |

### 5.3 Calculated values

**(a) c\*_eff / CEA の理論 c\*** [算出。5.2 の2つの出力値の比]

| g3 | Al5 | Al7 | g4 |
|---|---|---|---|
| 0.9961 → 0.9961 | **0.9051（≈ Yg）→ 0.9985** | **0.8474 → 0.9743** | **0.7844（≈ Yg）→ 0.9982** |

- 分母の CEA の値は、それぞれの run の平均燃焼圧で計算されている。その圧力は前と後で違う（5.2 の最後の行）。したがって、同じ圧力どうしの比較ではない

**(b) 修正前の本体が出力しない指標の再構成** [算出]
- 定義は、修正後の本体の `summarize()` と同じ式を、各本体の `time_history.csv` [出力] に当てはめたもの
- 修正後の側も同じ方法で算出している

| 指標 | 定義 | g3 | Al5 | Al7 | g4 |
|---|---|---|---|---|---|
| 初期の全質量 | `m_total_eos` の1行目 | 9.39986e-6 → 同じ | 1.58894e-5 → 1.43990e-5 kg（−9.38%） | 1.60157e-5 → 1.39214e-5（−13.08%） | 1.66361e-5 → 1.30680e-5（−21.45%） |
| 初期の気相質量 | `m_gas_equilibrium` の1行目 | 同じ | 1.43990e-5 → 1.30485e-5（−9.38%） | 1.39214e-5 → 1.21010e-5（−13.08%） | 1.30680e-5 → 1.02652e-5（−21.45%） |
| 局所 c\* | 燃焼中の `cstar` 列の中央値 | 1899.19 → 1899.16 | 1335.90 → 1474.16（+10.35%） | 1313.32 → 1510.83（+15.04%） | 1249.98 → 1592.65（+27.41%） |
| 全区間 c\* | 燃焼 ＋ blowdown の ∫p0·At dt / m_out | 1907.77 → 1907.77 | 1335.93 → 1474.03（+10.34%） | 1313.31 → 1510.62（+15.02%） | 1250.44 → 1594.05（+27.48%） |
| 範囲外の点の割合 | `extrap` 列の平均 | 0 → 0 | 0 → 0 | 0 → 0 | 0 → 0 |
| model_status | §5.0 の表の `dc073bf` 以降の列（順 3〜5。全 run が homogeneous_equilibrium・cea_python のため）を、table の Yg_min [出力] に当てはめたもの | in_scope_not_validated（前後とも） | approximate_with_limitations（前後とも） | 同左 | 同左 |

**再構成の検証:**
- 修正後の側で算出した値は、修正後の本体の出力値（`local_model_cstar_m_s`、`cstar_burn_interval_m_s`、`cstar_full_interval_m_s`、`property_extrapolation_fraction`、`validation.model_status`）と、差 ≤4.1e-10 で一致した（status は一致）
- 修正前の側で算出した燃焼区間 c\* は、修正前の出力値（`cstar_eff_m_s`）と ≤7.9e-10 で一致した
- 初期質量の前後の比 [算出] は 1/Yg に当たる（Al5: 1.1035 ≈ 1/0.906、g4: 1.2730 ≈ 1/0.785）

### 5.4 Interpretation
- g3 は Yg ≈ 1 なので、変化はほぼ 0（anchor）
- Al5 と g4 では、修正前の c\* 比が Yg に一致する。質量基準の系統誤差（1/Yg）と整合する
- Al7 では、修正前の比（0.847）が Yg（0.870）より約 2.6% 低い
  - 修正後に残る約 −2.6% は、CEA との比較で見られる定γの式の差（約 −2.2%）と同じ程度である
  - Al7 は、AL2O3 が燃焼室とスロートの間で液相 → 固相になる条件にある
  - **この相変化と偏差の因果は確定していない**

### 5.5 Unverified items and limitations
1. **本体の差は、質量基準の修正だけではない**
   - `35b599c` → `d13a9120` の間には、A-5.1-5〜A-8（cache・validation・台帳の許容差 S-2・eta_T0 の範囲 A-7・golden）と、未 commit の A-9〜A-11、P1〜P5、CLI v5 が含まれる
   - 表の差のすべてが、質量基準の修正によるとは限らない
2. **5.3(b) の値は再構成である。** 修正前の本体の出力ではない
   - とくに修正前の model_status は、修正後の規則を当てはめたもの。修正前の本体には status の出力が無い
3. 10/07 の「clamp の割合」は、`extrap` 列（範囲外の点の割合）で代用した。10/07 の定義と同じかは**未確認**
4. c\* 比の分母の圧力が、前後で違う（5.3(a)）
5. Al5 / Al7 は証拠候補で、golden ではない。G4a は 2e-3 で**未合格**のまま。この節は、どのケースの合格も示さない
6. R3（報告元の AP/Al、G4b）は、まだ入力が無い
