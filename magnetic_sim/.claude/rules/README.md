# .claude/rules/ — path-scoped 編輯規則

**用途**：每個 `.md` 是一條 Claude Code 工作規則，依工作情境（觸發片語 / 涉及的目錄）自動載入；強制規則須「動手前先讀全文」。對應觸發片語總表見 `../../CLAUDE.md` 的 Quick Triggers。

**內容**（逐條一行用途）：
- `main-workspace.md` — `magnetic_sim/ANSYS/main/`（Kuo Quadrupole）工作目錄：各類產物的正規輸出位置表 + 命名規則。
- `main-workflows.md` — `main/doc/workflows/` 各 SOP 的自然語觸發對應（出 STEP / 跑 FEM / 抽場 / fit / 畫圖…）。
- `ansys-cad-alignment.md` — ANSYS 幾何尺寸必對齊 CAD（STEP/IGES = source of truth）；改幾何前必量 CAD、不一致必通報。
- `result-read-safety.md` — 讀 ANSYS 結果三層防呆（讀前回報消歧、讀後核指紋、以 RESULTS_MAP 為準），避免讀錯 case。
- `sim-cleanup.md` — 清 sim 副產物 SOP：6 項不可影響 + 2 項不可失去能力、強制 dry-run、預設 half-clean、`--full` 須明確同意。
- `apdl-editing.md` — APDL 腳本編輯規則（`[ADDED]`/`[MODIFIED]` 標記、6 coil 只差 CURR_ARRAY、保留 `D,ALL,MAG,0`）。
- `comsol-livelink.md` — COMSOL LiveLink 連線通用法（拆兩 process：獨立 server + `matlab -batch` 內 `mphstart(2036)`）。
- `cad-import-ansys.md` — 複雜 CAD(STEP)→ANSYS 真實幾何匯入通用法（STEP→SpaceClaim `.x_t`→`ac4para` ANF→MAPDL `/INPUT`；別 primitive 硬拼、別 IGESIN/~PARAIN）。
- `deliver-step-for-check.md` — ANSYS 建完幾何**交付使用者疊 CAD 檢查一律出 STEP**（不出 ANSYS IGES：被 SW/OCC 讀成英吋 ×25.4，flag/name 都救不了）；用 OCC 產 mm STEP（`write.step.unit=MM`）。
- `fit-current-matches-sim.md` — 擬合電流必須等於 FEM 激發電流（目前 1A），不可把操作電流 0.6A 灌進 fit。
- `charge-model-source-convention.md` — 電荷模型符號慣例：每顆極激發時 B 一律從尖端射出（全 source）。
- `actuator-frame.md` — `Calibration_using_FEM_modeling/` 整包一律用 actuator（磁極軸）座標系；六極 `R_act` 旋轉、單極 +x 建模，不可停在 measure/WP frame。
- `calibration-transfer-matrix-output.md` — current_base / voltage_base 結果 PDF 除 K̄_I / D̄（無因次），一律加印 gauge 前物理轉移矩陣 ᴮH̄_I [mT/A] = ĝ_I·K̄_I 與 ᴮH̄_V [mT/mV] = ĝ_V·D̄；**附則（2026-08-04 改）**：USE_BIAS 時偏移**只印無因次 `e/ℓ̂` 一張**（不再印 µm 版）。
- `calibration-shared-structure.md` — `matlab/APDL/Calibration_using_FEM_modeling/` 三模型通用校正的 canonical 結構**凍結**（config/function/main/data/results 依 model 分層 + 檔案切分固定）；新增/改名/移動/刪除資料夾或檔、重組切分、加校正模型**一律先問**；改既有檔內文不受限。
- `hung-docs.md` — 涉及 `magnetic_sim/ANSYS/backup/hung/` 時的必讀文件清單（troubleshooting、build workflow…）。
- `hexapole-build.md` — 「建 hexapole」觸發流程：依序收集 4 個參數、嚴格照 build-workflow 執行。
- `iges-model-id.md` — 貼 `model_check/<topic>/*.iges` 路徑時，從路徑識別物理模型，不問「這是哪個模型」。
- `simulation-constraints.md` — APDL 幾何/材料/求解器硬約束（alpha=54.74° FIXED、tip 公式鎖定、元素型別、BC）。
- `unit-reference.md` — 單位統一慣例（ℓ̂ µm／b mT／ĝ_I mT/A／V mV／R_a A/Wb／力 pN），source of truth = doc/Unit Reference Sheet PDF。

**（2026-07-06 由 `magnetic_sim/ANSYS/main/rules/` 移入本層，改為全域自動載入）**：
- `read-rules-first.md` — 每次在 `main/` 下動手前，全部規則須在 context（移入本層後已自動載入）。
- `db-folder-retention.md` — `ANSYS_data/<model>/db/` 子夾只留 `.db` + 主 `.rmg`（無 digit），殘留禁留。
- `matlab-output-layout.md` — MATLAB `.mat` 放產生它的程式旁 `data/`（`matlab/<model>/<activity>/data/`；`MATLAB_data/` 已移除）。
- `results-pdf-only.md` — `…/Hall_sensor_base_fix_dir/results/` 只放 `.pdf`（`.mat/.tex/.aux` 不留）。
- `figure-style.md` — 圖表風格 preset；**畫圖前必先問使用者要哪個風格選項**（①粗體框圖：大字粗體/box/無 grid/tick 減半/單位 `()`；3D 框體 A 立方 daspect / B 異質軸 pbaspect 兩變體依幾何選）；10^0/無單位不標；直方圖 nb=180；**同類比較圖共用 colorbar/clim（禁各自 auto-scale）**。
- `figure-output.md` — 圖一律**輸出實檔**到 `figures/`；要改**原地改腳本→重跑→覆蓋同檔**迭代到定案（不丟 temp）。

**位置**：本資料夾 2026-07-06 由 `FEM_sim/.claude/rules/` 移到 `magnetic_sim/.claude/rules/`（使用者拍板，改放磁學模擬層）。
**相關**：git-root Claude 設定見 `../../../.claude/README.md`；觸發片語總表與 repo 總覽見 `../../../CLAUDE.md`。
