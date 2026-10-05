# FEM Simulation Workspace

`FEM_sim/` 是 FEM 模擬容器，session 一律開在這個根目錄。**本檔是全 repo 唯一的 CLAUDE.md。**
活躍工作全在 `magnetic_sim/ANSYS/main/`（下稱 `main/`），求解器有 APDL 與 Maxwell 兩條線。
`magnetic_sim/ANSYS/backup/` 是歸檔，活躍程式不得依賴它。

## 🔒 鐵則

1. **不擅自更動檔案架構**：未經使用者明確指示，不得移動 / 改名 / 刪除 / 新建資料夾 / 重組目錄 / 搬移檔案，一律**先問**。新增功能組資料夾也算。改既有檔的**內文**不受此限。
2. **改動同步 README**：改某夾內容 → 更新該夾 `README.md`；新增 / 改名 / 移動夾（須先問）→ 更新上層索引 README 與本檔「資料夾地圖」。範圍只限受影響的那幾份。
3. 所有新產物寫進 `main/` 下對應子目錄，不得寫到其他設計目錄或外部路徑。

## 資料夾地圖

`main/` 底下（每夾有自己的 `README.md`，動該夾前先讀）：

| 資料夾 | 放什麼 |
|---|---|
| `CAD_model/<model>/` | SolidWorks 原檔 + STEP（幾何 **source of truth**） |
| `model_check/<model>/` | 交付檢查用的 mm STEP |
| `apdl/<model>/{geom,mesh,sim,postproc}/` | APDL deck |
| `ANSYS_data/<model>/{data,db}/` | APDL 輸出：`.dat` 場、`.db` 模型；每個 model 有 `RESULTS_MAP.md` |
| `matlab/Flux/{APDL/Calibration_using_FEM_modeling,Maxwell}/` | 磁場校正（電流版 / 電壓版） |
| `matlab/Force/{APDL/Calibration_using_FEM_modeling,Maxwell}/` | 力模型校正 |
| `figures/` | 論文圖與其繪圖腳本（`paper_fig/`、`paper_fig_plot/`） |
| `reference/` | LaTeX 原稿 + PDF + 論文 |

MATLAB 校正包結構：`config/`（模型設定）、`function/` + `main/`（校正求解）、`data/`（校正結果 `.mat`，只留最終定案版）、`results/`（PDF）、`utils/scripts/<model>/`（額外計算的腳本）、`utils/data/<model>/`（額外計算的結果；Flux/Maxwell 已依模型分夾，其餘仍是 `utils/{scripts,data}/`）、`plot/`（只畫圖）、`figures/<model>/{current,voltage}/`。與校正無關的分析一律走 `utils/scripts/<model>` → `utils/data/<model>` → `plot`，不另開暫存資料夾。

Maxwell 專案在 repo 外：`D:\Maxwell_sim\<model>\{cad,project,scripts,export,materials}\`（`export/` 是匯出的 `.fld` 場）。

資料流：

```
CAD_model (STEP) ─┬─ apdl deck → MAPDL → ANSYS_data/*.dat ─────────┐
                  └─ D:\Maxwell_sim project → AEDT → export/*.fld ──┤
                                                                    ↓
              matlab/{Flux,Force}/…/main → data/*.mat → utils → plot → figures/*.png、results/*.pdf
```

## Quick Triggers（動手前先讀全文）

規則檔都在 `.claude/rules/`。

| 觸發 | 規則檔 |
|---|---|
| 清 db / 清 sim 副產物 / 清 ANSYS results / 磁碟滿 | `ansys-db-cleanup.md`（`geom/` 整層刪、`mesh/` 保主 `.db`+主 log、`sim/` 再加主 `.rmg`；**刪 `geom/` 前必先掃 >100MB 的錯置網格 db**；保留白名單、強制 dry-run 與使用者批准） |
| 動校正包的結構 / 輸出 / 座標系 / 編號 / 單位 | `calibration-shared-structure.md`（結構凍結、改先問）、`calibration-transfer-matrix-output.md`、`actuator-frame.md`、`pole-coil-numbering.md`、`unit-reference.md` |
| 畫圖 / 改圖 / 調字級 | `figure-style.md` |
| 寫 / 改任何 `plot/**` 腳本 | `plot-scripts-pure.md`（只准 load `.mat` → 畫圖） |
| 動 `results/` | `results-pdf-only.md`（只放 PDF） |
| 新建檔案 / 取檔名 / 命名 variant | `short-names.md` |
| 想複製一支腳本做變體 | `modify-existing-files.md`（改現有那支） |
| 腳本需要模型常數 / 看到 `addpath(backup)` | `no-backup-data.md`（一律走 `model_config(...)`） |

## ANSYS 可用性

執行 MAPDL 前先確認路徑存在：`G:\ANSYS Inc\v252\ansys\bin\winx64\MAPDL.exe`。不存在就先搜尋其他磁碟再跑。

## Hexapole Design Constraints (Mandatory)

These constraints apply to ALL hexapole designs in this repo. They are non-negotiable.

1. **Orthogonal pair axes**: 3 opposing pole pairs (P1-P2, P3-P4, P5-P6) must have mutually perpendicular connecting lines
2. **Tips on common sphere**: All 6 pole tips at distance R_norm from the sphere center (R_norm is adjustable)
3. **60-degree azimuthal offset**: Upper layer rotated 60 deg relative to Lower layer
4. **alpha = arctan(sqrt(2)) = 54.74 deg is FIXED**: derived from constraints 1-3, not a free parameter
   - `R_norm_xy = R_norm * sqrt(2/3)` and `R_norm_z = R_norm / sqrt(3)` — these formulas are locked
   - Lower poles at 0, 120, 240 deg; Upper poles at 60, 180, 300 deg

## Rules
- 6 Coil scripts are synchronized: only `CURR_ARRAY` values differ (one coil = 1, rest = 0)
- All code comments in English; explanations to user in Traditional Chinese
- Mark all APDL changes with `[ADDED]` or `[MODIFIED]` comments
- Always verify `D,ALL,MAG,0` boundary condition exists before `/SOLU`
- Preserve original commented-out code (prefixed `!****`) unless asked to remove
- Use tab indentation matching original style

## Prohibitions
- NEVER commit ANSYS output files (*.rst, *.db, *.full, etc.)
- NEVER change geometry parameters without explicit user approval
- NEVER modify element types (SOLID96, SOURC36) or material properties without approval
- NEVER remove boundary condition section (`[ADDED]` block near line 500)
- NEVER 跑任何 db / sim 清理（rm intermediates / rm result dir）前未先讀 `.claude/rules/ansys-db-cleanup.md` 全文；違反 = 違規
- NEVER 刪 `db/geom/` 前未先掃「>100MB 的 .db」——那是錯置的**網格** db（孤本、含指紋基準網格），必須先搬進 `db/mesh/`
- NEVER 寫 `rm -f <jobname>*` 清 ANSYS 檔（會連主 `.rmg` + `.db` 一起刪）——用規則裡的針對性 pattern
- NEVER change alpha (54.74 deg) or the R_norm_xy / R_norm_z formulas
- NEVER produce a pole configuration that violates pair-axis orthogonality

## Notation Standard
All symbols and terms follow Fei Long's 2016 dissertation (B, Phi, q, K_I, rho, R_a, g_I, N_c, ...).

- Use **paper pole names** (P1-P6) in all user-facing text, figures, and discussion; APDL coil indices (1-6) only in APDL code and raw data context
- **Coil → pole mapping is PER-MODEL, not a global rule** — it is each deck's build order, not physics.
  long2016 `[1,3,6,5,2,4]` / NTU `[1,3,6,5,2,4]` / **hung `identity`**; new models: identity (coil k = Pk).
  Never copy another model's map — see `.claude/rules/pole-coil-numbering.md`.
- **禁用 "WP" 這個字眼**（使用者拍板 2026-08-06）：圖、軸標、圖例、註解、對話一律不用。那一點就叫**原點**（六極尖共球球心 = 繪圖座標原點），圖上以黑點標示、不加文字。既有檔不強制回溯清理，動到哪個檔就順手改。
- Two meanings of rho: physical (500 um) vs fitted (900 um) — always clarify which
- Units: ANSYS outputs Tesla; figures use mT（完整單位表見 `unit-reference.md`）
