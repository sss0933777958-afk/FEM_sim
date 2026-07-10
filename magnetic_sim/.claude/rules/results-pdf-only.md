# 規則 #3：`Hall_sensor_base_fix_dir/results/` 只放 PDF

**使用者拍板（2026-06-26）**：
`matlab/long2016_hexapole_halfcut/Calibration_using_FEM_modeling/Hall_sensor_base_fix_dir/results/`
這個資料夾**只留最終 `.pdf`**（矩陣 / d 向量等 xelatex 排版輸出）。

## 🔒 規則
- 只留 `*.pdf`；**`.mat / .tex / .aux / .log / .txt` 一律不留**。
  - `.mat` 分析成果 → 改放同包 `data/`（見 `matlab-output-layout.md`）。
  - `.tex / .aux / .log` 排版中間檔 → xelatex 編完即刪。
- driver（`main.m` / `main_interp.m` / `main_Dmatrix.m` / `calib_gap100um.m` / `decouple/solve_DH_*.m`）若會寫 `.mat`/`.txt` 到 `results/`，**改寫到 `data/`**；PDF 才留 `results/`。
- 同步該夾 `results/README.md`。

## 觸發片語
- 動到 `Hall_sensor_base_fix_dir/results/` 內容時
- 「results 只放 pdf」/「清 results 中間檔」

## 何時不適用
- 其他功能組的 `results/`（如 `fix_dir/results`、`bs_matrix/results`）目前仍放 `.tex` 排版原稿——本規則只綁這一個 `Hall_sensor_base_fix_dir/results/`。

相關：`matlab-output-layout.md`、memory `reference_local_latex_compile`、`feedback_matlab_local_data_layout`。
