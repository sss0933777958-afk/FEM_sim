# 規則 #4：每次執行都要先讀規則（meta-rule）

**使用者拍板（2026-06-26）**：在 `magnetic_sim/ANSYS/main/` 下動手做**任何**工作之前，**先把 `rules/` 資料夾內的全部規則檔讀一遍**，再開始動作。

## 🔒 規則
- 只要是在 `main/` 樹下工作——跑 sim、寫/改分析腳本、搬/建/刪檔、畫圖、清理、改幾何、出 IGES、跑 COMSOL…——**第一步先讀 `main/rules/` 全部 `.md`**（至少掃過標題與本檔下方清單，命中的規則讀全文）。
- 不靠記憶猜規則內容；規則可能已更新。
- 全域規則（repo-root `…/FEM_sim/.claude/rules/`）仍適用；本資料夾是 main/ 專屬補充。

## 目前 `main/rules/` 規則清單
| 規則 | 何時必讀 |
|---|---|
| `db-folder-retention.md`（#1） | 清/整理 `ANSYS_data/<model>/db/` |
| `matlab-output-layout.md`（#2） | 寫/搬 MATLAB `.mat`（一律放 `matlab/<model>/<activity>/data/`） |
| `results-pdf-only.md`（#3） | 動到 `Hall_sensor_base_fix_dir/results/` |
| `read-rules-first.md`（#4，本檔） | 每次開工 |
| `figure-style.md`（繪圖風格） | 畫任何圖之前（**先問使用者要哪個風格選項**） |
| `figure-output.md`（繪圖輸出） | 畫任何圖（**輸出實檔→覆蓋迭代到定案**） |

## 觸發片語
- 任何 `main/` 下的工作起手式（等同「開工前 checklist」）。

相關：`README.md`（本資料夾索引）、`main/CLAUDE.md` 最上方鐵則、memory `feedback_matlab_local_data_layout`。
