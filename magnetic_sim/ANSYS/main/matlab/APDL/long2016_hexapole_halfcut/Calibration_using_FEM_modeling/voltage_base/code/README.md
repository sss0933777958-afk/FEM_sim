# …/voltage_base/code/ — 電壓側（V̄, D̄, ^Bĝ_V）校正程式碼

**用途**：Hall-sensor 電壓校正主程式（USE_BIAS 統一 fix / 18-param bias），依角色分三組。
（2026-07-16 重整：`main.m` 末端直接產 PDF；`main_function/`=main pipeline、`function/`=與 main 無掛鉤的工具。）

**內容**：
- `main/` — driver `main.m`（六極電壓校正）+ `run_calib_D_bias_gap_calibrate.m`（gap_calibrate 的 18-param 變體 driver）。
  流程：`load_coils_actuator` → `select_ball` → all-source → `fitting` → `build_S_matrix`（profile G=D^v）
  → `build_V_matrix` → `solve_D_bar_gain` → 存 `.mat` → **`emit_model_results` + `gen_B_matrix` 產 2 份 PDF**。
- `main_function/` — **main 校正流程函式**：電壓 pipeline + emit（見該夾 README）。沿用 current_base 的 fitting/build_S_matrix。
- `function/` — **與 main 無掛鉤的工具/診斷函式**（`extract_Vmat`/`extract_Vmat_interp_center`/`extract_Vmat_elemB`/`make_scaled_coil1`/`compute_vmatrix`）。
- `plot/` — 繪圖 / 診斷腳本。

**資料流**：`main.m` 讀 `ANSYS_data/.dat` + sensor 局部網格 CSV → 存 `../data/calib_D_*.mat`
→ `emit_model_results`/`gen_B_matrix` 出 `../results/{single_param|eighteen_param}/*.pdf`。

**相關**：見上層 `../README.md`。
