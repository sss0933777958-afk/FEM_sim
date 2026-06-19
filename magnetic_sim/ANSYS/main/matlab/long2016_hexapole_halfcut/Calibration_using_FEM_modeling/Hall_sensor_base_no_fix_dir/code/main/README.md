# …/Hall_sensor_base_no_fix_dir/code/main/ — 求 d 主程式（18 參數 bias 版）

**用途**：`main.m` — 18 參數 bias 版 Hall-sensor per-pole `d` 的一條龍 driver。config 在頂部（R_select=150 µm、I=1 A、S_hall=130 V/T）。

**流程**：
1. 載入 `calib_bias.mat`（no_fix_dir 校正）→ 取 `ell_hat, Pc_18, R, F, apdl_to_paper_idx`。
2. 重載 6-coil FEM（`wp`）→ R_select 近場球 → 旋進 actuator 框 → `build_A` → `M=A'A`、`c=A'Bstack`（page-1）。
3. `build_sensor_geometry` → `extract_Vmat`（all-source）→ `solve_d` → `sensor_residual_bias`（actuator 框殘差）（page-2）。
4. 存 `calib_sensor_d_no_fix_dir.mat`（`MATLAB_data/.../charge_fit/fitting_d/`，不蓋 fix-ℓ 版）。
5. `write_d_tex`（d_v2/d_final）+ `compute_KH`+`write_KH_tex`（KH_v2/KH_final）→ `../../results/`。

**預期數值**：ell_hat≈0.857 mm（= calib_bias.mat）；sensor 模型相對 RMSE ≈8.4%。

**命名 / 慣例**：單一主程式組 → `code/main/main.m`；模型數學一律在 `../function/`；不重跑 18 維 fminunc（載入結果）。

**相關**：見上層 `../README.md`、`../function/README.md`。
