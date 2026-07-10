# …/Hall_sensor_base_fix_dir/code/main_function/ — main 校正流程函式

**用途**：`main.m` 校正流程相關函式（原在 `../function/`）：main 直接呼叫的 pipeline + 輸出階段 `emit_model_results`（獨立跑、讀 main 存的 `.mat` 出 PDF）。
**內容**：
- `build_sensor_geometry.m` — 6 顆 sensor 中心 `sensor_pos` + 法線 `sensor_n`（全域 WP 座標；下極磨平面、上極錐面 CAD 傾角）。
- `extract_Vmat_interp.m` — 內插版 Vmat：粗網格在 sensor 圓柱內均勻撒點（預設 1000）做 FEM tet 重心內插，對 B·n+ 平均 + all-source 翻號。
- `emit_model_results.m` — 輸出階段（獨立跑）：讀 `../data/calib_D_*.mat` 出結果 PDF（V/D̄/^Bĝ_V/Ĥ_V）到 `../results/`。

**慣例**：一檔一函式。**Hall_sensor_base_no_fix_dir 的 main 也 addpath 本夾**（借用這兩支）。main.m 另借 no_fix 的 `../../no_fix_dir/code/main_function`（load_coils_actuator/select_ball/build_A）。

**相關**：見 `../README.md`。
