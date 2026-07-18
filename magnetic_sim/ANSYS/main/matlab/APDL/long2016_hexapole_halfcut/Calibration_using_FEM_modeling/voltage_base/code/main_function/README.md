# …/voltage_base/code/main_function/ — main 校正流程函式（唯一函式夾）

**用途**：`main.m` 電壓側校正流程函式，一步一支；`main.m` 逐個呼叫，末端在此產 PDF。
（2026-07-16 重整：`function/` 已併入本夾，改名簡化。）

**內容**：
- `build_sensor_geometry` — 6 顆 sensor 中心 `sensor_pos` + 法線 `sensor_n`（全域 WP 座標）。
- `build_V_matrix` — sensor 電壓矩陣 V̄（sensor 圓柱內均勻撒點 tet 重心內插、B·n+ 平均 ×S_hall、all-source；原 `extract_Vmat_interp`）。
- `solve_D_bar_gain` — 由 G、V̄ 解 Ĥ_V → gauge D̄（D̄(1,1)=5/6）、^Bĝ_V（新增；與 current 側 `solve_KI_bar_gain` 對稱）。
- `emit_model_results` — **main 最後一步呼叫**：讀 `../data/calib_D_*.mat` 出主結果 PDF（D̄/G/V̄/ℓ̂/^Bĝ_V，18-param 加印 ê）到 `../results/{single_param|eighteen_param}/`。
- `gen_B_matrix` — **main 另呼叫**：出 B 矩陣（V_out/V_in）PDF 到同 `../results/` 子夾。
- `emit_tex` — 共用 LaTeX 排版 helper。

**沿用 current_base**：main 另 addpath `../../current_base/code/main_function`，借 `load_coils_actuator`/`select_ball`/`build_S_matrix`/`fitting`。
**已刪除**：`build_S`、`solve_d`、`sensor_residual(_bias)`、`sensor_cost_lhat`（死鏈）。**搬到 `../function/`**（與 main 無掛鉤的工具/診斷）：`extract_Vmat`、`extract_Vmat_interp_center`、`extract_Vmat_elemB`、`make_scaled_coil1`、`compute_vmatrix`。

**相關**：見 `../README.md`。
