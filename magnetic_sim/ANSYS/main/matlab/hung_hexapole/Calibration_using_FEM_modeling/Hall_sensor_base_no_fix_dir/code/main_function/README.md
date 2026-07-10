# …/Hall_sensor_base_no_fix_dir/code/main_function/ — main 校正流程函式

**用途**：`main.m` 校正流程相關函式（本包）。
**內容**：
- `emit_model_results.m` — 輸出階段（獨立跑）：讀 `../data/calib_D_gap_200um.mat` 出結果 PDF（V/D̄/^Bĝ_V/Ĥ_V）到 `../results/`。

**慣例**：`main.m` 的 pipeline（`load_coils_actuator`/`select_ball`/`fit_bias`/`make_Pc`/`build_A` 等）**借用** `../../no_fix_dir/code/main_function` + `../../Hall_sensor_base_fix_dir/code/main_function`（`build_sensor_geometry`/`extract_Vmat_interp`），故本夾只放本包專屬的輸出階段 `emit_model_results`。

**相關**：見 `../README.md`。
