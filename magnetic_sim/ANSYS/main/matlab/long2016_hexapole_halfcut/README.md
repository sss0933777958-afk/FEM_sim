# matlab/long2016_hexapole_halfcut/ — Long Fei 下極半切六極（halfcut hexapole）的 MATLAB 分析母夾

**用途**：此 model（Long Fei 2016 下極半切 hexapole）所有 MATLAB 分析的根。FEM 在 ANSYS 求解、結果（`.dat`）放 `../../ANSYS_data/long2016_hexapole_halfcut/`，本夾的程式讀那些場做點電荷校正（charge calibration）、B_S 矩陣、sensor 配置、場視覺化等，並把成果 `.mat` 寫回 `../../MATLAB_data/`、auto-gen `.tex` 寫到各組 `results/`、圖寫到各組 `figures/`。

**內容**（第二層 = 功能組 activity）：
- `common/` — ★ resolver（`ansys_path.m` 讀 FEM `.dat/.db`、`matlab_path.m` 讀寫 `.mat`）；不要硬寫絕對路徑。
- `Calibration_using_FEM_modeling/` — 點電荷校正 + Hall-sensor d 主流程（單一主程式組 `{fix_dir, no_fix_dir, Hall_sensor_base_fix_dir, Hall_sensor_base_no_fix_dir}/code/main/main.m`）。
- `fixl_fit/`、`bias_fit/`、`bs_matrix/`、`validation/` — 多腳本組（`code/scripts/`）。
- `field_viz/`、`sensor_placement/` — 純繪圖組（`code/plot/`）。
- `field_cancellation/` — 掃 64 種 source/sink 電流組合，找工作空間磁場抵銷（`code/main/` driver + `code/plot/`，僅 `figures/`、無 `results/`）。
- 每組另含 `figures/`（`.png`/`.gif`）與 `results/`（auto-gen `.tex`）。

**資料來源 / 流向**：`ANSYS_data/<model>/<case>/*.dat`（用 `ansys_path`）→ MATLAB 算 fit/矩陣/場圖 → `MATLAB_data/<model>/<功能>/*.mat`（用 `matlab_path`）+ 各組 `results/*.tex` + 各組 `figures/*.png`。

**命名 / 慣例**：model-first + 功能組 schema；每組 `code/{main|scripts|plot}/` + `function/`；單一主程式組用 `code/main/main.m`。電荷模型一律全 source（B 從尖端射出）；擬合電流必須 = FEM 激發電流（1 A）。

**相關**：見上層 `../README.md`、`../../CLAUDE.md`（schema、resolver、繪圖腳本規則）。
