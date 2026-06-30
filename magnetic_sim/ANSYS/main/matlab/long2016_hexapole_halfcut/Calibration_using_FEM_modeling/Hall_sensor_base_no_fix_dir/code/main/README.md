# …/Hall_sensor_base_no_fix_dir/code/main/ — 主程式（18 參數 bias 版）

**用途**：
- `main.m` — 18 參數 bias 版 Hall-sensor per-pole `d` 的一條龍 driver。config 在頂部（R_select=150 µm、I=1 A、S_hall=130 V/T）。
- `main_Dmatrix.m` — **18-param bias 版 6×6 D 矩陣** driver（與 `Hall_sensor_base_fix_dir/main_Dmatrix.m` 的 on-axis 版並列，差別只在 ℓ̂-fit 用 `fit_bias`）。gap200um_mueq、100 點內插 Vmat（修正 sensor 位置、graded csv）、all-source（只翻 sink→D^v 對角全正）。產 `D^v, D̄, D, V, ê` PDF 到 `../../results/D_<variant>.pdf`、`.mat` 到 `../../data/calib_D_<variant>.mat`。ℓ̂≈0.857 mm、region err≈0.46%、recon~1e-15。沿用 `no_fix_dir/code/function`（load_coils_actuator/select_ball/fit_bias/make_Pc/build_A）+ `Hall_sensor_base_fix_dir/code/function`（build_sensor_geometry/extract_Vmat_interp）。

**流程**：
1. 載入 `calib_bias.mat`（no_fix_dir 校正）→ 取 `ell_hat, Pc_18, R, F, apdl_to_paper_idx`。
2. 重載 6-coil FEM（`wp`）→ R_select 近場球 → 旋進 actuator 框 → `build_S` 逐點 → `M=Σ S_iᵀS_i`、`c_j=Σ S_iᵀb_ij`（page-1）。
3. `build_sensor_geometry` → `extract_Vmat`（all-source，真實節點）→ `solve_d`（no-gain）→ `sensor_residual_bias`（actuator 框 cost J）（page-2）。
4. 存 `calib_sensor_d_no_fix_dir.mat`（`d`、`Vmat`、cost `J`；本組 `../../data/`，規則#2，不蓋 fix-ℓ 版）。

**輸出對齊** `Hall_sensor_base_fix_dir`：只回 `d`、`Vmat`、cost `J`（無 g_H / K_H / RMSE / LaTeX）。

**預期數值**：ell_hat≈0.857 mm（= calib_bias.mat）；輸出 cost `J`（無 RMSE）。注意：sensor 抽取已改真實節點，現網格每顆 sensor 僅命中 1 點（待網格加密）。

**命名 / 慣例**：單一主程式組 → `code/main/main.m`；模型數學一律在 `../function/`；不重跑 18 維 fminunc（載入結果）。

**相關**：見上層 `../README.md`、`../function/README.md`。
