# main/ — 活躍工作區權威規則（給 Claude）

`magnetic_sim/ANSYS/main/` 是 FEM 模擬的**活躍工作區**。第二層為 model topics：
`long2016_hexapole_halfcut`（主力）/ `kuo_quadrupole`（4-pole MEMS Quadrupole，原 kuo/）/ `zhang_quadrupole`（CAD 用 `long_fei`）。

**本檔 = 在 main/ 工作時的權威綱要**（cwd 在 main/ 或其子夾時自動載入）。

> 🗂 分層（都自動載入、由內而外疊加）：
> - **本檔 `main/CLAUDE.md`** — 活躍工作區完整規則（權威）。
> - `../../../CLAUDE.md`（`FEM_sim/`）— 容器架構 + **backup/ 樹**共通約束（跨設計規則兩邊鏡像，改一處要同步）。
> - `../../.claude/rules/*.md`（28 條細則）— 觸發式細則；見下「Quick Triggers」索引，動到特定東西前**先讀全文**。
> - memory（`~/.claude/projects/G--my-workspace/`）— Claude 跨 session 記憶（索引常駐、細節按需）。
>
> ⚠ **開新 session 請開在 `main/`（或子夾）**，才吃得到本檔 + 上面全部；開在更上層（FEM_sim / workspace 根）會漏掉本檔。

---

## 🔒 鐵則

### 1. 不擅自更動檔案架構
未經使用者明確指示，不得**移動 / 改名 / 刪除 / 新建資料夾 / 重組目錄 / 搬移檔案** → 一律**先問**。
- **例外**：改既有檔的**內文**（程式 / 文件內容）不受此限。
- 繪圖「新增功能組資料夾」也算架構變動 → 依繪圖規則，**開新組前先問**。

### 2. 改動同步 README
改了東西就把**受影響的** README 一起更新：改某夾內容 → 更新該夾 `README.md`；新增/改名/移動夾（須先問）→ 更新上層索引 README + 本檔「資料夾架構地圖」。範圍 = 受影響那幾份。

---

## Hexapole Design Constraints (Mandatory)
適用**所有** hexapole 設計、不可協商：
1. **Orthogonal pair axes**：P1-P2 / P3-P4 / P5-P6 三對連線互相垂直。
2. **Tips on common sphere**：6 尖端距 WP 中心 = R_norm（可調）。
3. **60° azimuthal offset**：上層相對下層轉 60°。
4. **alpha = arctan(√2) = 54.74° FIXED**（由 1–3 導出，非自由參數）：`R_norm_xy = R_norm·√(2/3)`、`R_norm_z = R_norm/√3` 公式鎖定；下極 0/120/240°、上極 60/180/300°。

## APDL / 通用 Rules
- 6 個 coil 腳本**同步**，只差 `CURR_ARRAY`（激發極=1、其餘=0）。
- code 註解一律**英文**；對使用者解說一律**繁體中文**。
- APDL 改動標 `[ADDED]` / `[MODIFIED]`；`/SOLU` 前必驗 `D,ALL,MAG,0` 邊界；保留 `!****` 註解碼（除非要求刪）；tab 縮排沿原樣。
- 未經批准不改幾何參數 / 元素型別（SOLID96、SOURC36）/ 材料。
- 用 dissertation notation；對外一律 **P1–P6**，APDL coil index 只在改 APDL code / raw 脈絡時提。

## 🎨 繪圖（強制，畫任何圖前先讀）
- **先問功能組 + 風格選項**（風格 preset 見 `figure-style.md`），不自己猜、不憑記憶。
- **一任務一腳本、原地改**：同一支改到定案；**定案前不另開**；使用者沒說「新增」就不開第二支（一組可多支，但各對一張已定案的圖）。
- **一律輸出實檔到該組 `figures/` → 覆蓋同檔迭代**到定案（不丟 temp preview 等定案才落地），見 `figure-output.md`。
- 不屬於任何現有功能組 → 要新增功能組資料夾（`matlab/<model>/<新組>/code/plot/` + `figures/`）→ **先問**。
- **場圖畫真實 FEM 節點原值，不用 scatteredInterpolant / 格點內插**（除非使用者明確要求 → 且圖說標「內插」）。減量用「每格挑最近 y=0 平面節點」（仍是節點原值）。見 memory `plot-real-nodes`。

## Prohibitions
- NEVER commit ANSYS 輸出（`*.rst / *.db / *.full` …）。
- NEVER 未經批准改幾何參數 / 元素型別 / 材料屬性。
- NEVER 改 alpha（54.74°）或 R_norm_xy / R_norm_z 公式；NEVER 產生違反對極正交的配置。
- NEVER 移除 BC 區塊（`[ADDED]` block）。
- NEVER 清 sim（rm intermediates / result dir）前沒讀 `sim-cleanup.md` 全文；NEVER 用 `--full` 除非使用者**明確同意**（預設 half-clean，保 `.db` + 主 `.rmg`）；NEVER 繞過 helper 手刻 `rm` ANSYS 檔。
- NEVER 改 ANSYS 幾何 / `mt_constants` 前沒量對應 CAD；不一致**不自己選值**，通報使用者拍板（`ansys-cad-alignment.md`）。
- NEVER 用內插畫場圖 / 把內插當 raw 呈現（除非明示）。

## Notation Standard
全依 Fei Long 2016 dissertation。canonical glossary：`backup/hexapole-long2016/docs/notation-glossary.md`。
- 對外一律 **P1–P6**；APDL coil index 只在 APDL code / raw 脈絡。
- **coil→paper map 是 per-model（非全域）**：long2016 `[1,3,6,5,2,4]` / NTU `[1,3,6,5,2,4]` / hung `identity`；新 model 用 identity。**別互抄**（`pole-coil-numbering.md`；舊全域宣稱是錯的、曾靜默污染 hung K̄_I）。
- **禁用 "WP" 這個字眼**（使用者拍板 2026-08-06）：圖、軸標、圖例、註解、對話一律不用（含「WP 框 / WP 區」）。那一點就叫**原點**（六極尖共球球心 = 繪圖座標原點），圖上以黑點標示、不加文字。既有檔不強制回溯清理，動到哪個檔就順手改；全 repo 清理要先問。詳見 `figure-style.md`「座標原點與『WP』字眼」。
- ρ 兩義：physical（500µm）vs fitted（900µm）—— 講清哪個。
- 單位：ANSYS 出 Tesla；WP 場圖用 mT；dissertation Fig 2.4 用 Gauss。（完整單位表見 `unit-reference.md`。）

---

## 資料夾架構地圖（9 個頂層夾）

每夾都有自己的 `README.md`，要動該夾前先讀。

| 資料夾 | 是什麼 / 放什麼 | 要找資料去這裡 |
|---|---|---|
| `CAD_model/` | SolidWorks 原檔 + STEP（幾何 **source of truth**） | 量尺寸、出圖前對齊 CAD |
| ~~`IGES/`~~ | **已刪除（2026-07-13）**：交付一律走 STEP（`deliver-step-for-check`）；IGES 中繼走 scratch（`ANSYS_data/<model>/db/geom_hexvariants/`）。⚠ 部分舊 export deck 仍 `IGESOUT` 到此路徑 → 重跑前需先建目錄或改指 scratch/model_check | 交付看 `model_check/` STEP |
| `model_check/` | 交付檢查用 mm STEP（+ 舊 IGES）；`IGESIN` 匯入用 | **交付幾何檢查、建 mesh 前匯入** |
| `apdl/` | APDL 腳本：`<model>/{geom,sim,postproc}/`（+ sweep） | 改幾何/參數/重跑 sim 的 input |
| `ANSYS_data/` | FEM 輸出 `<model>/<case>/`（`.dat` 場 / `.db` 模型 / `.cdb`） | **讀 FEM 結果**（.dat） |
| `matlab/` | MATLAB 分析碼 `<model>/<功能組>/code/...` + `figures/` + `results/` + **`data/`（`.mat` 成果）** | 跑分析、畫圖、**讀/寫 `.mat`** |
| ~~`MATLAB_data/`~~ | **已移除（2026-06-26）**：遷到各活動 `matlab/<model>/<activity>/data/`；`matlab_path()` deprecated | `.mat` 改去 `matlab/.../data/` |
| `doc/` | LaTeX 原稿 + 編譯 PDF + `workflows/`(SOP) | 推導、報告、流程 SOP |
| `.claude/` | Claude Code 本地設定（`settings.local.json`） | 非工作產物，通常不動 |

## 資料流 pipeline（一條龍）

```
CAD_model (SLDPRT/STEP)
   → apdl geom → (IGES 中繼走 scratch) → OCC → model_check/*.step (交付檢查)
   → apdl/<model>/geom + sim (APDL input)
   → [ANSYS MAPDL 求解]
   → ANSYS_data/<model>/<case>/*.dat (場) + *.db (模型)
   → matlab/<model>/<功能組>/code (讀 .dat 做 fit/矩陣/校正/畫圖)
   → matlab/<model>/<功能組>/data/*.mat (成果)  +  .../figures/*.png (圖)  +  .../results/*.pdf
   → doc/<主題>/ (LaTeX/PDF 報告)
```

## Resolver（路徑解析，不要硬寫絕對路徑）
在 `matlab/<model>/common/`，相對自身定位：
- `ansys_path('<model>'[, 'coilN', ...])` → `ANSYS_data/<model>/...`（讀 FEM `.dat/.db`）。
- ~~`matlab_path(...)`~~ **deprecated**：`.mat` 一律放各活動 `matlab/<model>/<activity>/data/`，用腳本自身 root 變數算 `fullfile(<本組夾>,'data')`。

## matlab/ 功能組 schema
`matlab/<model>/` 第二層是**功能組**（activity）：
- 單一主程式組 → `code/main/main.m`；多腳本組 → `code/scripts/`；純繪圖組 → `code/plot/`。
- 每組另有 `figures/`、`results/`（PDF）、**`data/`（`.mat` 成果，見 `matlab-output-layout.md`）**。

---

## ANSYS 可用性 + 典型指令
本機 MAPDL：`G:\ANSYS Inc\v252\ansys\bin\winx64\MAPDL.exe`（不存在則搜標準位置 `C:\Program Files\ANSYS Inc\<ver>\...`）。

```bash
ANSYS="G:\ANSYS Inc\v252\ansys\bin\winx64\MAPDL.exe"
# 單顆 coil（batch，無 GUI）；long2016 從 backup/hexapole-long2016/ 跑
cd magnetic_sim/ANSYS/backup/hexapole-long2016
"$ANSYS" -b -np 4 -m 24000 -dir "results/coil1" -j "coil1" \
  -i "$(pwd)/apdl/MT_Modeling_Geometry_Meshing_Solving_Coil1.txt" -o "results/coil1/solve.out"
# 6 顆連續：for i in 1..6，改 coil${i}
```
> ⚠ long2016 的 APDL deck 在 `backup/hexapole-long2016/`，在那工作時載的是 `FEM_sim/CLAUDE.md`（非本檔）。

---

## Quick Triggers → `.claude/rules/`（觸發式細則，動手前先讀全文）
（全文在 `../../.claude/rules/`，索引 `../../.claude/rules/README.md`）

| 觸發 | 規則檔 |
|---|---|
| 讀/載入 ANSYS 結果、抽 .dat、算矩陣/fit/畫場圖前 | `result-read-safety.md`（三層防呆 + `ANSYS_data/<topic>/RESULTS_MAP.md`） |
| 清 sim / 清 results / 磁碟滿 | `sim-cleanup.md`（先讀全文；預設 half-clean） |
| 清 `db/` 夾 | `db-folder-retention.md`（只留 `.db`+主 `.rmg`） |
| 寫/搬 `.mat` | `matlab-output-layout.md`（放程式旁 `data/`） |
| 動 `results/` | `results-pdf-only.md`（只放 PDF） |
| 畫任何圖 | `figure-style.md`（先問風格）+ `figure-output.md`（輸出實檔覆蓋迭代） |
| 改 ANSYS 幾何 / mt_constants / 對齊 CAD | `ansys-cad-alignment.md`（CAD=source of truth） |
| 匯入 CAD/STEP 進 ANSYS | `cad-import-ansys.md`（STEP→x_t→ac4para→/INPUT） |
| 交付幾何給使用者檢查 | `deliver-step-for-check.md`（一律出 STEP） |
| charge fit / 改 I_actual | `fit-current-matches-sim.md`（I = FEM 激發 1A） |
| Calibration_using_FEM_modeling 結構 | `calibration-shared-structure.md`（結構凍結、改先問）、`calibration-transfer-matrix-output.md`、`actuator-frame.md`、`charge-model-source-convention.md`、`pole-coil-numbering.md`、`unit-reference.md` |
| 出 STEP / 解析 STEP / 建 APDL 幾何 / 檢查模型 / 跑 COMSOL 等 SOP | `doc/workflows/`（入口 `workflows/README.md`；部分為 kuo 時代遺留） |
| 產物落點 / 資料夾架構 | 本檔「資料夾架構地圖」+ 各夾 `README.md` |

## Detailed Docs（long2016 深技術）
`backup/hexapole-long2016/docs/`：`fitting-methods.md`（[B-6x] 最終）、`model-validation.md`、`notation-glossary.md`、`coil-winding-sign-convention.md`、`charge-model-fitting.md`、`ansys-environment.md`、`simulation-parameters.md`、`workflow.md`、`troubleshooting.md`。

## Compact Instructions（context 壓縮時保留）
1. 6 腳本只差 `CURR_ARRAY`；2. `D,ALL,MAG,0` 對 DSP solver 必須；3. long2016 結果在 `backup/hexapole-long2016/results/coilN/`；4. 對使用者用繁體中文；5. 一律 P1–P6；6. notation 依 Long 2016（`notation-glossary.md`）；7. alpha=54.74° FIXED。
