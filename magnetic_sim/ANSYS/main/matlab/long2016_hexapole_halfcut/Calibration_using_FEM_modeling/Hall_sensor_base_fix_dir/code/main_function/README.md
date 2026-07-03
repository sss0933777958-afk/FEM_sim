# …/Hall_sensor_base_fix_dir/code/main_function/ — main.m 掛勾的函式

**用途**：`main.m` 掛勾、原在 `../function/` 的 pipeline 函式（2026-07-03 重構）。
**內容**：
- `build_sensor_geometry.m` — 6 顆 sensor 中心 `sensor_pos` + 法線 `sensor_n`（全域 WP 座標；下極磨平面、上極錐面 CAD 傾角）。
- `extract_Vmat_interp.m` — 內插版 Vmat：粗網格在 sensor 圓柱內均勻撒點（預設 1000）做 FEM tet 重心內插，對 B·n+ 平均 + all-source 翻號。

**慣例**：一檔一函式。**Hall_sensor_base_no_fix_dir 的 main 也 addpath 本夾**（借用這兩支）。main.m 另借 no_fix 的 `../../no_fix_dir/code/main_function`（load_coils_actuator/select_ball/build_A）。

**相關**：見 `../README.md`。
