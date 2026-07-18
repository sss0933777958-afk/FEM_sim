# …/current_base/code/ — 電流側（K̄_I）校正程式碼

**用途**：六極點電荷校正主程式（USE_BIAS 統一 fix / 18-param bias），依角色分三組。
（2026-07-16 重整：`main.m` 末端直接產 PDF；`main_function/`=main pipeline、`function/`=與 main 無掛鉤的工具。）

**內容**：
- `main/` — driver `main.m`（六極）+ 單極次要 driver `main_singlepole.m` / `main_singlepole_bias.m`。
  六極流程：`load_coils_actuator` → `select_ball` → `build_S_matrix` → `fitting` → `solve_KI_bar_gain`
  → `calc_ellipsoid` → 存 `.mat` → **`emit_model_results` 產 PDF**。
- `main_function/` — **main 校正流程函式**（main 直接呼叫）：pipeline + 輸出階段 emit（見該夾 README）。
- `function/` — **與 main 無掛鉤的工具函式**（`svd_SiHI_WP`、`sweep_singlepole_tipcut`）。
- `plot/` — 繪圖 / 診斷腳本（部分含自帶 local `make_Pc`）。

**資料流**：`main.m` 讀 `ANSYS_data/.dat`（經 `ansys_path`）→ `main_function/` 算 → 存 `../data/*.mat`
→ `emit_model_results` 讀 `.mat` 出 `../results/{single_param|eighteen_param}/*.pdf`；`plot/` 畫圖。

**相關**：見上層 `../README.md`。
