# Calibration_using_FEM_modeling/ 結構凍結（強制讀取）

**使用者拍板（2026-07-22）**：`matlab/APDL/Calibration_using_FEM_modeling/` 已是**三模型（long2016 / hung / NTU）
通用校正管線**的 canonical 結構。**此結構凍結 —— 資料夾分層與檔案切分不得擅自更動；要改「結構」一律先問使用者。**

當工作涉及：
- 在 `Calibration_using_FEM_modeling/` 底下**新增 / 改名 / 移動 / 刪除**資料夾或 `.m` 檔
- 重組檔案切分（把某函式拆檔 / 併檔）、改 `config`/`function`/`main`/`data`/`results` 分層
- 為校正管線**新增一個模型**、或搬動校正檔

**動手前先讀完此規則。** 改既有檔的**內文**（見「何時不適用」）不受此限。

對應 memory：`feedback_calibration_shared_structure.md`
相關規則：`no-structure-change-without-ask`(通則)、`modify-existing-files.md`、`results-pdf-only.md`、`matlab-output-layout.md`。
相關 memory：[[project_calibration_shared_pipeline]]（管線內容與各檔職責）。

---

## 🔒 Canonical 結構（凍結，as of 2026-07-22）

```
Calibration_using_FEM_modeling/
  config/<model>/[<geom>/]mt_constants.m   per-model[/幾何變體] 設定 + 幾何（R_act/Pc_base/sensor）+ 方法旗標
                                       （long2016 分 tip40um/tip400um 子夾；hung/NTU flat。<model>=long2016.../hung.../NTU...）
  function/
    extract_ansys_data.m   純讀 .dat（raw, measure frame, Tesla）
    build_actuator_data.m  raw → actuator frame 資料包 ad（讀 cfg.R_act/Pc_base；**不再假設魔術角/canonical**；非 D̄ 矩陣）
    fitting.m              優化器（e 開關 USE_BIAS；build_S / make_Pc 為 local；幾何無關）
    solve_current.m        current 解 K̄_I / ĝ_I / 𝒞κ
    solve_voltage.m        voltage 解 D̄ / ĝ_V / 𝒞κ
    build_V_matrix.m       電壓提取（sensor 幾何 override>cfg>內建；V_METHOD csv-tet|scattered）
    model_config.m         dispatcher（model_config(model,geom)；variant 子夾+flat fallback；clear mt_constants 防快取）
    emit_tex.m / emit_results.m   低階 LaTeX helper / 讀自描述 .mat → PDF
    interp_field_to_points.m  變體 WP 場內插到參考 geom 點雲（INTERP_TO；2026-07-23 由 utils 移入）
    import_ansys_data.m / filter_iron_nodes.m   extract/build_actuator_data 的 helper（2026-07-23 由 backup copy 入）
  main/main.m              driver（頂部 per-run 調參 + BASE = current | voltage 開關）
  common_path/ansys_path.m           共用路徑 resolver（model-agnostic；2026-07-23 由 long2016/common 搬入）
  data/<model>/{.mat,csv}/           結果存放（.mat 自描述：結果 + 設定條件與參數；含保存的 single-pole .mat）
  results/<model>/{eighteen,single}/ PDF 輸出（eighteen = USE_BIAS true、single = false）
  utils/*.m                          **通用函式**（model-agnostic，`main.m` 不直接呼叫的第二層 helper）
                                     例：`pole_sensor_geometry.m`（sensor 幾何唯一來源；2026-08-08 使用者拍板放此層）
  utils/<model>/                     後處理 compute 腳本：讀 calib .mat → 進一步運算（如 SVD 分解）→ 存 computed .mat 到 data/（**不重跑校正**）
  plot/<model>/<base>/               plot 腳本：讀 .mat → 畫圖到 figures/（**不做重運算**；<base>=current|voltage，不分 param）
  figures/<model>/<base>/{single,eighteen,common}/   三模型圖檔（single=fix、eighteen=bias、common=比較/診斷/場）
```
> 2026-07-23（大搬家清理，per-model 三專案將刪）：**plot + figures** 併入本夾（`plot/`、`figures/`）、**全部 `.mat`**
> 併入 `data/<model>/.mat/`（80 個）、`ansys_path.m` 入 `common_path/`、`import_ansys_data`/`filter_iron_nodes` 入
> `function/`。**本夾已完全自足、live pipeline 不再 addpath `backup/`**（禁止跨專案）。刪除範圍限
> `matlab/APDL/{long2016,hung,NTU}_hexapole*/`；`backup/` 樹不刪（本夾已 copy 所需、不再依賴）。
>
> **2026-07-23 幾何通用化**：pipeline 不再假設魔術角。幾何（`R_act`/`Pc_base`/per-pole `sensor_pos/n`）+ 方法旗標
> （`v_method`/`interp_to`/`r_loc`/`sensor_r_loc`）**全由 `config/<model>/<geom>/mt_constants.m` 提供**，`build_actuator_data`/`build_V`
> 只消費。新增「幾何變體」= 加 `config/<model>/<geom>/`（照樣 mt_constants，提供 R_act/Pc_base/sensor/旗標）——**仍先問**。
> tip400 已由 bespoke driver（run_tip_calib，已刪）折成 `config/.../tip400um/` + main.m 旗標（`GEOM='tip400um'`）。
>
> **2026-07-24 後處理分工（使用者拍板）**：校正之後的後處理（如空間 SVD 分解 σ1/σ2/σ3、gain/iso…）走**三段式**、
> compute 與 plot 分離：
> ① **校正** = `main/main.m` + `function/`（頂部旋鈕切 MODEL/VARIANT/BASE/USE_BIAS）→ 產 `calib_*.mat` 到 `data/<model>/.mat/`。
> ② **後處理 compute** = `utils/<model>/<name>.m`：**讀** ①的 calib .mat → 逐點運算（SVD…）→ **存** computed `.mat`
>    到 `data/<model>/.mat/`（自描述）。**不重跑校正**（校正歸 ①）。
> ③ **畫圖** = `plot/<model>/<base>/<name>.m`：**讀** ②的 .mat → 提取 → 畫圖到 `figures/<model>/<base>/<param>/`。**不做重運算**。
> 資料流：`main+function (校正→data/.mat)` → `utils/<model>/ (compute→data/.mat)` → `plot/<model>/<base>/ (讀→figures)`。
> 範例：`utils/hung_hexapole/svd_sigma_hung_ws.m`（讀 calib→算 σ→存 `svd_sigma_hung_ws_R150.mat`）+
> `plot/hung_hexapole/current/plot_sigma_hist_hung_ws.m`（讀 σ.mat→σ1/σ2/σ3 疊圖）。
> ⚠ `utils/hung_hexapole/run_hung_current.m` 是**舊式**（utils 內自己重跑校正+inline 畫圖），不符此分工、暫留不動、勿沿用。
>
> **2026-08-08 使用者拍板：`utils/` 根層放 model-agnostic 通用函式**（`<model>/` 子夾維持後處理 compute 腳本，
> 兩者共存不衝突）。首例 = `utils/pole_sensor_geometry.m`（六顆 Hall sensor 位置與外法線的**唯一來源**，
> 兩分支各一份、內容相同）。`main.m` 已加 `addpath(fullfile(CAL,'utils'))`。
> ⚠ 曾一度建 `lib/` 放它，**已作廢刪除**——通用函式一律 `utils/` 根層，不要再開 `lib/`。
> （歷史對照：2026-07-23 曾把通用的 `interp_field_to_points.m` 由 utils **移入** `function/`；
> 該檔維持在 `function/` 不動，新的通用函式則放 `utils/` 根層。）

## 🔒 規則

1. **禁止擅自更動「結構」**：`function/` 的檔案切分、`config`/`data`/`results` 依 `<model>` 的分層、`main/` 的位置
   —— 任何**新增/改名/移動/刪除資料夾或檔案、重組切分**都要**先問使用者**。
2. **新增模型**：照既定 pattern 加 `config/<model>/mt_constants.m` + `data/<model>/{.mat,csv}/` +
   `results/<model>/{eighteen,single}/` —— 這是既定樣式，但**仍先跟使用者確認**再建。
3. **檔名 / 職責已定案**（別自作主張改名或搬邏輯）：`extract_ansys_data` 只純讀、`build_actuator_data` 只轉 actuator frame、
   優化器 = `fitting`（含 build_S/make_Pc local）、解參數分 `solve_current`/`solve_voltage`、電壓提取 = 一支
   `build_V_matrix`（sensor 幾何當內部 local）、輸出獨立 `emit_tex`+`emit_results`。
4. **落點固定**：`.mat` → `data/<model>/.mat/`、PDF → `results/<model>/{eighteen,single}/model_results_<base>_<variant>.pdf`。

## 何時不適用（不算「結構變動」，可逕改）

- 改**既有檔內的 code**（bug fix、在既有檔內加功能 / 分支 / 參數）—— 那是內文，不是結構。
- 產物（`.mat` / `.pdf` / `.csv`）落進既有 `data/`、`results/` 夾。
- 三模型各自舊的 `<model>/Calibration_using_FEM_modeling/{current_base,voltage_base}/` 副本（那是**別的**、待退役的
  per-model 舊碼，不是本 canonical 共用夾）。

## 觸發片語（任一即套用）

- 「動 / 改 Calibration_using_FEM_modeling 結構」「校正共用夾」「加校正模型」「搬 / 改名校正檔」
- 你正打算在 `Calibration_using_FEM_modeling/` 新增/改名/移動/刪除 `.m` 或資料夾時 → 停，先問。
