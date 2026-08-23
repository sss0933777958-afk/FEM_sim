# .claude/rules/ — path-scoped 編輯規則

**用途**：每個 `.md` 是一條 Claude Code 工作規則，依工作情境（觸發片語 / 涉及的目錄）自動載入；強制規則須「動手前先讀全文」。對應觸發片語總表見 `../../CLAUDE.md` 的 Quick Triggers。

**內容**（21 條，逐條一行用途）：
- `ansys-cad-alignment.md` — ANSYS 幾何尺寸必對齊 CAD（STEP/IGES = source of truth）；改幾何前必量 CAD、不一致必通報。
- `result-read-safety.md` — 讀 ANSYS 結果三層防呆（讀前回報消歧、讀後核指紋、以 RESULTS_MAP 為準），避免讀錯 case。
- `ansys-db-cleanup.md` — **（2026-08-17 由 `sim-cleanup.md` + `db-folder-retention.md` 合併）** `ANSYS_data/<model>/db/` 三層政策：`geom/` **整層刪**（可由 deck 重生）、`mesh/` 保主 `.db`+主 log、`sim/` 再加主 `.rmg`；⚠ **刪 `geom/` 前必先 `find -size +100M`** —— 實測 20 顆網格 db（29.2 GB、19 顆孤本、含指紋基準網格）被錯放在該層；含保留白名單、主檔 vs worker 判別、`rm -f <job>*` 陷阱、GUI scratch 夾整夾刪、強制 dry-run + 使用者批准。
- `apdl-editing.md` — APDL 腳本編輯規則（`[ADDED]`/`[MODIFIED]` 標記、6 coil 只差 CURR_ARRAY、保留 `D,ALL,MAG,0`）。
- `cad-import-ansys.md` — 複雜 CAD(STEP)→ANSYS 真實幾何匯入通用法（STEP→SpaceClaim `.x_t`→`ac4para` ANF→MAPDL `/INPUT`；別 primitive 硬拼、別 IGESIN/~PARAIN）。
- `deliver-step-for-check.md` — ANSYS 建完幾何**交付使用者疊 CAD 檢查一律出 STEP**（不出 ANSYS IGES：被 SW/OCC 讀成英吋 ×25.4，flag/name 都救不了）；用 OCC 產 mm STEP（`write.step.unit=MM`）。
- `fit-current-matches-sim.md` — 擬合電流必須等於 FEM 激發電流（目前 1A），不可把操作電流 0.6A 灌進 fit。
- `charge-model-source-convention.md` — 電荷模型符號慣例：每顆極激發時 B 一律從尖端射出（全 source）。
- `actuator-frame.md` — `Calibration_using_FEM_modeling/` 整包一律用 actuator（磁極軸）座標系；六極 `R_act` 旋轉、單極 +x 建模，不可停在 measure/WP frame。
- `calibration-transfer-matrix-output.md` — current_base / voltage_base 結果 PDF 除 K̄_I / D̄（無因次），一律加印 gauge 前物理轉移矩陣 ᴮH̄_I [mT/A] = ĝ_I·K̄_I 與 ᴮH̄_V [mT/mV] = ĝ_V·D̄；**附則一（2026-08-04 改）**：USE_BIAS 時偏移**只印無因次 `e/ℓ̂` 一張**（不再印 µm 版）。~~附則二（K̄(2,2) 固定 0.8340）~~ **已於 2026-08-17 依使用者指示徹底廢除、整段刪除**。**附則三（2026-08-15 立、2026-08-17 改定義）**：擬合誤差 PDF 一律印 **NMAE 向量範數版**（`[Σ_jΣ_i‖b_ij − S_i·ᴮĝ·M̄·u_j‖/N_p] / b̄ ·100`，先逐點取 3 維殘差長度再相加；舊的分量 L1 版與 RMSPE 都**停用**、只留 `.mat` 的 `NMAE_L1`/`RMSPE` 追溯，三代數字互不可比，實測新/舊 L1 比值 0.73~0.76）；四支 solve + 兩支 emit + 兩支 main 兩分支都要同步；含**換指標時的回填法**（從 `.mat` 設定重建取樣點、以「RMSPE 逐位相同」驗證）；另 sensor 距非 4.572mm 時 soff tag 要併進 `VAR_OUT`（否則 PDF 被覆蓋）。
- `calibration-shared-structure.md` — `matlab/APDL/Calibration_using_FEM_modeling/` 三模型通用校正的 canonical 結構**凍結**（config/function/main/data/results 依 model 分層 + 檔案切分固定）；新增/改名/移動/刪除資料夾或檔、重組切分、加校正模型**一律先問**；改既有檔內文不受限。
- `pole-coil-numbering.md` — 磁極/coil 編號：paper 名 P1–P6 全域通用，但 **APDL coil→paper 的 map 是 per-model**（long2016/NTU `[1,3,6,5,2,4]`、hung identity），**不可互抄**；驗證看 K̄_I 對角占優 + 全正。
- `gap-nogap-folder-convention.md` — 氣隙變體資料夾**只准 `gap_200um` + `no_gap` 兩個名字**（跨 data/db/apdl 各層），禁止增殖平行 gap 名夾。
- `modify-existing-files.md` — 要改腳本一律**改現有那支**（首選參數化），不要複製成 `_v2`/變體檔；產物落新夾則允許。
- `simulation-constraints.md` — APDL 幾何/材料/求解器硬約束（alpha=54.74° FIXED、tip 公式鎖定、元素型別、BC）。
- `unit-reference.md` — 單位統一慣例（ℓ̂ µm／b mT／ĝ_I mT/A／V mV／R_a A/Wb／力 pN），source of truth = reference/Unit Reference Sheet PDF。
- `no-backup-data.md` — **禁止使用 `backup/` 的資料與程式**（禁 `addpath(backup)`、禁讀其 `.m/.mat/.dat`）；模型設定一律走 live `model_config(...)`。起因＝`paper_fig_plot/` 有 12 支吃 `backup/.../mt_constants.m`，看不到 CAD 實測的真實錐體幾何。

**（2026-07-06 由 `magnetic_sim/ANSYS/main/rules/` 移入本層，改為全域自動載入）**：
- `matlab-output-layout.md` — MATLAB `.mat` 放產生它的程式旁 `data/`（`matlab/<model>/<activity>/data/`；`MATLAB_data/` 已移除）。
- `results-pdf-only.md` — Calibration 的 `results/<model>/{single,eighteen}/` 只放 `.pdf`（`.mat/.tex/.aux` 不留）。
- `figure-style.md` — 圖表風格 preset；**畫圖前必先問使用者要哪個風格選項**（①粗體框圖：大字粗體/box/無 grid/tick 減半/單位 `()`；3D 框體 A 立方 daspect / B 異質軸 pbaspect 兩變體依幾何選）；10^0/無單位不標；直方圖 nb=180；**同類比較圖共用 colorbar/clim（禁各自 auto-scale）**；**圖例一律照 `plot_sensor_B_hist.m` 範本**（northoutside + 寬度切齊框 + 黑粗框 + FS24 + 兩欄「系列｜統計值」）。
- `figure-output.md` — 圖一律**輸出實檔**到 `figures/`；要改**原地改腳本→重跑→覆蓋同檔**迭代到定案（不丟 temp）。

**位置**：本資料夾 2026-07-06 由 `FEM_sim/.claude/rules/` 移到 `magnetic_sim/.claude/rules/`（使用者拍板，改放磁學模擬層）。
**相關**：git-root Claude 設定見 `../../../.claude/README.md`；觸發片語總表與 repo 總覽見 `../../../CLAUDE.md`。
