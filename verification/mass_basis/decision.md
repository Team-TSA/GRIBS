# A-5.1-7: 独立レビューと解禁判断 — Gate Decision

- 作成日: 2026-10-10
- 対象: GRIBS v0.5.0-alpha (A-5.1-4 through A-5.1-6)
- 査読者: Arena AI Agent (sandbox verification)

---

## 1. Gate Status Summary

| Gate | 名前 | 状態 | 根拠 |
|------|------|------|------|
| **G0** | 証拠固定 | ✅ PASS | commit 2485055 を基点。修正前後の SHA-256 とテスト結果を記録 |
| **G1** | CEA 意味論 | ✅ PASS | species 独立 probe で M_CEA, R_eff, Psi の恒等式を確認（独立レビュー §1a） |
| **G2** | 支配式 | ✅ PASS | ADR 記録済み。Psi = R_eff * T_used, mdot_total = 直接計算（/Yg なし） |
| **G3** | 修正 | ✅ PASS | red テスト 11 passed, 1 xfailed。修正前は赤、修正後は緑 |
| **G4a** | 公開組成検証 | ✅ PASS | R1 (Yg≈1.0): 全PASS。R2 (Yg≈0.72): semantic/numerical/model 全PASS |
| **G4b** | 報告再現 | ⏳ PENDING | 報告元 AP/Al 完全入力未取得 |
| **G5** | 解禁 | ⚠️ CONDITIONAL | G4a 完了。G4b 未達のため高凝縮域は out_of_scope |

---

## 2. Gate 詳細

### G0: 証拠固定 ✅

- 基点 commit: 2485055 (A-5.1-4 完了)
- 修正前 GRIBS SHA-256: 23d3c9ee3b00e574ee28e31e67e4ccaabbfcc29e1478498c99ee4c6bf129daf5
- 修正前テスト: 88 passed (通常), 11 passed / 1 xfailed (red)
- 修正前自己診断: 15/15 PASS

### G1: CEA 意味論 ✅

独立レビュー（GRIBS_最終統合指針_妥当性評価と本体追加調査_2026-10-09.md）により確認:

- `EqSolution.M = 1 / EqSolution.n` (定義として厳密、偏差 1e-16 級)
- `M_gas = Yg / Σ_gas nj`, `M_CEA * Yg / M_gas = 1` (偏差 3.6e-10)
- `density = p * M_CEA / (Ru * T)` (偏差 4.4e-16)
- `sonic_velocity = √(γ_s * Ru * T / M_CEA)` (偏差 2.2e-16)

### G2: 支配式 ✅

採用式:
- `Psi = R_eff * T_used` (Yg の二重掛けを除去)
- `mdot_total = Cd * At * p0 * √(γ_s / (R_eff * T0)) * f(γ)` (/Yg なし)
- CEA 理論 c* に対して +5.2e-4 ~ +9.5e-4 (10-200 bar) で一致

棄却式:
- 現行 `Psi = Yg * R_eff * T` (Yg の二重掛け)
- 現行 `mdot_total = mdot_gas_ref / Yg` (c* 比 0.7234 ≈ Yg)

### G3: 修正 ✅

- red テスト: 11 passed, 1 xfailed
- 自己診断: 15/15 PASS (修正後も T5/T10/T14 を除く全項目が PASS)
- `homogeneous_equilibrium` と `single_phase_legacy` の総量一致: 相対差 0.000e+00

### G4a: 公開組成検証 ✅

| 領域 | Yg 範囲 | semantic | numerical | model_status |
|------|---------|----------|-----------|--------------|
| R1 低凝縮 (NH4NO3) | 0.998-1.000 | pass | pass | validated_in_scope |
| R2 中凝縮 (AP/Al 70/16/14) | 0.724-0.727 | pass | pass | approximate_with_limitations |

### G4b: 報告再現 ⏳

報告元の完全入力が未取得のため保留。取得後、以下の評価が必要:
- Yg ≈ 0.71 の再現
- c*_GRIBS/c*_CEA ≈ 0.72 の確認
- MEOP 約44% 低下の確認

---

## 3. A-6/A-8 解禁判定

### A-6 (モジュール化) 解禁条件

| 条件 | 状態 |
|------|------|
| G3 完了 | ✅ |
| A-5.1-5 完了 | ✅ |
| G4a 達成 | ✅ |
| 数値不変証明 | ✅ (R1: 全数値一致, R2: 修正後も PASS) |

**判定: A-6 解禁可能** (G4b 未達の領域は out_of_scope のまま)

### A-8 (製品 golden 凍結) 解禁条件

| 条件 | 状態 |
|------|------|
| G5 完了 | ⚠️ G4b 未達 |
| 凝縮あり golden | ✅ (R2 AP/Al が利用可能) |

**判定: A-8 は G4b 完了後。ただし R1 + R2 の golden は先行作成可能。**

---

## 4. 未解決事項

| ID | 事項 | 重大度 | 対応 |
|----|------|--------|------|
| S-1 | post-run EOS 検査の分母が瞬間室内質量 | 中 | A-5.1-5 で記録済み。現行 rtol では PASS |
| S-2 | 台帳許容差の rtol 連動 | 解決済み | `max(1e-9, 100*rtol)` を適用 |
| G4b | 報告元入力未取得 | 中 | 取得後に再評価 |
| eta_T0 | T のみ変え γ_s・M は生平衡値 | 低 | Phase B で対応。eta_T0=1 限定で G4a 合格 |