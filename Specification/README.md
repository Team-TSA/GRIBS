# GRIBS Phase-A starter kit

「GRIBS改良計画_評価と詳細手順.md」の **Phase A（v0.6.0-alpha）** に着手するための
最小スケルトンです。書き直しではなく「現行エンジンを包んで、公開 API と結果オブジェクトを
先に凍結する」ための足場として作ってあります。

```
gribs_phaseA_kit/
├── pyproject.toml                  # パッケージ定義・pytest マーカー
├── cases/                          # 基準ケース（決定論的な設定の例）
│   ├── ref_l1_converging.json      #   ε=1（現行版でも完走する基準）
│   ├── ref_l2_cd_nozzle.json       #   ε=4（現行版では T4 で停止 = F-01 の再現）
│   └── ref_l1_wide_window.json     #   物性窓を 0.1–20 MPa に広げた版（クランプ回避）
├── schema/gribs_config.schema.json # 設定スキーマ（WP A-7 の雛形、2020-12）
├── src/gribs/
│   ├── __init__.py                 # 公開 API: run_simulation / load_config / Result
│   ├── api.py                      # run_simulation の実装（現行エンジンをラップ）
│   └── result.py                   # Result / Warning データクラス（凍結対象）
└── tests/
    └── test_phase_a_acceptance.py  # 受入テスト（台帳・Result 契約・スキーマ・C-D 完走）
```

## 使い方

```bash
# 1) 依存（CEA 公式パッケージ。無くても reference/単相試験は動く）
python3 -m pip install cea pytest jsonschema

# 2) 現行エンジンの場所を教える（既定は ../gribs_baseline/gribs.py）
export GRIBS_LEGACY_PATH=/home/user/gribs_baseline/gribs.py

# 3) API を使ってみる
python3 - <<'PY'
import sys; sys.path.insert(0, "src")
from gribs import run_simulation
r = run_simulation("/home/user/gribs_baseline/gribs_config_eps1.json")
print(r.summary_line())
print("ledger:", {k: f"{v:.6e}" for k, v in r.ledger.items()})
for w in r.warnings:
    print(f"  [{w.severity}] {w.code}: {w.message} {w.context}")
PY

# 4) 受入テスト
python3 -m pytest -q tests
#   期待: 2 passed, 1 skipped, 1 xfailed
#   xfailed = WP A-2 の対象（現行版は Ae/At > 1 で自己試験 T4 が失敗して abort する）
```

## ここで実装した Phase A の成果物（対応表）

| ファイル | 対応 WP | 内容 |
|---|---|---|
| `src/gribs/result.py` | A-6 | 結果オブジェクトの契約（履歴・要約・警告・入力・CEA 情報・収束情報・台帳・モデルと仮定・版・ハッシュ） |
| `src/gribs/api.py` | A-4, A-6, A-9 | `run_simulation(config) -> Result`、グローバル状態の保存/復帰（並行実行可能性）、項目別質量台帳、警告コード（`W_PROPERTY_RANGE`, `W_ENTRAINMENT`, `W_NOZZLE_TRANSITION`） |
| `schema/gribs_config.schema.json` | A-7 | 設定スキーマ（単位・範囲・列挙・相互依存のコメント付き） |
| `tests/test_phase_a_acceptance.py` | A-2, A-4, A-6, A-7 | 受入テスト。F-01（ε>1 が走らない）を xfail として固定し、修正したら strict 失敗で気づけるようにしてある |

## このキットが意図的に「未完成」な点（= 残りの作業）

1. `api.py` は現行の単一ファイルエンジンを import している。**パッケージ分割（A-6）完了時に
   `_load_legacy()` を消して `gribs.simulation.engine` に置き換える**だけで、公開 API は変えない。
2. 歴史列のキー名は現行流儀（`sample()` は単位なし、CSV は単位付き）を両対応にしている。
   A-9 で片方に統一すること（片方に寄せたらこの互換コードを消す）。
3. `schema_version` は `0.6.0-alpha` を期待している。A-0/A-7 完了まではスキーマテストが
   skip される（現行設定は `0.4.4-alpha` のまま）。
4. GUI 用の `cli.py` は未作成（現行 `main()` を `gribs/cli.py` へ移すのが A-6 の一部）。

## 拡張の指針

- **必ず `run_simulation` 経由でしか計算しない**（GUI から内部へ触らせない）。lint で強制してよい。
- 数値を変える変更（物理モデル・補間・後処理）は、まず `tests/golden/` に基準を凍結してから行う。
- 警告は必ずコード付きで追加する（表示は消費側の責務）。
