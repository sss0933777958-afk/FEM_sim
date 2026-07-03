# …/Hall_sensor_base_fix_dir/code/ — 求 d 程式

**用途**：Hall-sensor 每極 `d` 的全部程式碼，分層：
- `main/` — 主程式：`main.m`（**唯一統一 driver**：single-parameter on-axis，fit ℓ̂ → G(=D^v) → `extract_Vmat_interp`（下極 −β、Ø0.3×0.1 圓柱、10k 內插）抽 V → 解 Ĥ_V/D̄/^Bĝ_V（論文 notation）；**存 `../data/calib_D_gap200um_mueq.mat`**（.mat field 名沿用 Dmat/g_V）+ console；**V/D̄/^Bĝ_V/Ĥ_V 的 PDF 由 `function/emit_model_results.m` 產生**（`../results/model_results_gap200um_mueq.pdf`）。已取代舊 `main_Dmatrix`/`main_Vmat`/`main_interp`）。
- `main_function/` — **main.m 掛勾**：`build_sensor_geometry`（6 sensor 位置/法線）、`extract_Vmat_interp`（內插版 Vmat，圓柱內均勻撒 n_uniform 點預設 1000、tet 重心內插）。Hall_no_fix 的 main 也 addpath 本夾借用這兩支。
- `function/` — 其餘模型/IO 輔助（`build_S`、`extract_Vmat`、`extract_Vmat_interp_center`、`extract_Vmat_elemB`、`compute_vmatrix`、`sensor_cost_lhat`、`sensor_residual`、`solve_d`、`emit_model_results`）。
- `decouple/` — 滿 6×6 解耦 `D_H`：`solve_DH_full.m`（真實節點 Vmat）、`solve_DH_interp.m`（**內插版 1000 點 Vmat → D_H 6×6 + V_j 對角 6×6×6**；cost_J = fix_dir 自由電荷下界，存 `calib_DH_interp.mat`）。
- `plot/` — 場視覺化 + sign 診斷（9 支：plot_P2sensor_Braw_P1exc / plot_P1P2_air_circuit_3d / plot_P2pole_circuit_2d / plot_P2sensor_tets_3d / plot_interp_tet_schematic / diag_P2_Bn_map / diag_P2P1_single / diag_Vmat_sign(_center)）。**逐支用途見 `plot/README.md`**。圖存 `../figures/`。

**相關**：見上層 `../README.md`。
