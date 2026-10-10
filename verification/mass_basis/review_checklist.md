# A-5.1-7: Review Checklist

- 作成日: 2026-10-10
- 対象: GRIBS v0.5.0-alpha (A-5.1-4 through A-5.1-6)

---

## 1. Physical Correctness

- [x] `M_CEA = 1 / EqSolution.n` (total mixture mass basis)
- [x] `R_eff = Ru / M_CEA = Yg * R_gas` (恒等式確認)
- [x] `Psi = R_eff * T_used` (Yg の二重掛けなし)
- [x] `mdot_total` = 直接計算 (/Yg なし)
- [x] `density = p * M_CEA / (Ru * T)` (全質量密度)
- [x] `sonic_velocity = √(γ_s * R_eff * T)` (CEA と 2.2e-16 で一致)
- [x] Yg=1 単相極限で修正前後が数値一致

## 2. Cache Integrity

- [x] cache_schema = 4 (A-5.1-4 で設定)
- [x] schema 3 を明示拒否
- [x] `two_phase_mode` を cache key に追加
- [x] `mass_basis_semantics_version` を cache key に追加
- [x] 旧 cache は新しい key と異なるため自動的に無視される

## 3. Configuration Migration

- [x] v0.4.3 → v0.5.0 移行時に `W_MIGRATION_MODEL_CHANGE` を warnings に出力
- [x] `molecular_weight: total_MW` を拒否
- [x] 移行 notes が summary.json に保存される

## 4. Output Completeness

- [x] `summary.json` に `validation` ブロック追加
- [x] `semantic_status`, `numerical_status`, `model_status` を出力
- [x] `validated_pressure_range_Pa`, `validated_Yg_range` を出力
- [x] `mass_basis_semantics_version` を出力
- [x] `two_phase_capabilities` を出力 (機械可読)
- [x] c* 指標を再編 (local/burn/full/legacy)
- [x] summary テキストに validation status セクション追加

## 5. Metadata Corrections

- [x] `specific_gas_constant_definition` を訂正 (「GAS-PHASE」→「total mixture」)
- [x] `M_definition`, `R_definition`, `psi_definition`, `density_definition` を追加
- [x] `gamma_s_definition` を追加
- [x] `two_phase_capabilities` を追加 (modelled vs not modelled)

## 6. Numerical Gates

- [x] post-run EOS mass check: PASS (R1, R2)
- [x] mass ledger: PASS (R1, R2 with S-2 fix)
- [x] self-test: 15/15 PASS
- [x] pytest: 89 passed (通常), 11 passed / 1 xfailed (red)
- [x] 台帳許容差を rtol に連動 (`max(1e-9, 100*rtol)`)

## 7. Test Coverage

- [x] `test_corrected_mass_basis_metadata_contract` (A-5.1-5)
- [x] `test_migration_model_change_warning_code_in_source` (A-5.1-5)
- [x] モック Backend に `cache_schema`, `mass_basis_semantics_version` を追加
- [x] 台帳テストの許容差を新しい式に合わせて更新

## 8. Documentation

- [x] `verification/mass_basis/decision.md` (本ファイル)
- [x] `verification/mass_basis/before_after.md`
- [x] `verification/a5-1-6/regime_comparison.json`
- [x] `verification/a5-1-6/before_after_comparison.md`

---

## 9. Known Limitations

| Item | Status | Action |
|------|--------|--------|
| G4b: 報告元入力未取得 | PENDING | 取得後に再評価 |
| eta_T0 ≠ 1 の場合の γ_s・M 整合 | 記録済み | Phase B で対応 |
| 高凝縮域の実験検証 | 未実施 | Phase B |
| particle slip, finite-rate condensation | 未実装 | Phase B |
| post-run 検査の分母 (瞬間室内質量) | 記録済み | 将来修正 |

---

## 10. Final Verdict

**A-5.1-4 through A-5.1-6: COMPLETE**

- G0-G4a: All PASS
- G4b: PENDING (input not available)
- A-6: 解禁可能 (G4b 未達領域は out_of_scope)
- A-8: G4b 完了後に解禁

The CEA mass-basis correction is verified for:
- Single-phase limit (Yg ≈ 1.0): 数値不変
- Medium condensation (Yg ≈ 0.72): 物理的に正しい修正
- High condensation (Yg < 0.65): out_of_scope (G4b pending)