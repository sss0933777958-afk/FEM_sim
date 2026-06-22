# MATLAB_data/long2016_hexapole_halfcut/charge_fit/fitting_d — Hall-sensor 每極常數 d

**用途**：存放下極半切六極 **Hall-sensor 模型每極常數 `d`** 的求解結果。
與同層 `fit_KI_ball/`（各半徑 K_I 擬合）、`fitting_trend/`（R 掃描趨勢）並列。

**內容**：
- `calib_sensor_d.mat`：**fix_dir**（電荷鎖極軸）求出的**對角** d；模型 `b = S·V·d`（**no-gain，無 g_H**），含 `d/Vmat/exc_sign/ell_hat/J`。
- `calib_sensor_d_no_fix_dir.mat`：**no_fix_dir**（18-param bias、actuator 框）求出的 d（該版仍含 g_H）。
- `calib_DH_full.mat`：**滿 6×6 解耦矩陣 D_H**（`q=D_H·V`，由 `Hall_sensor_base_fix_dir/code/decouple/solve_DH_full.m` 產生）；cost_J = fix_dir 自由電荷下界(~0.0043)，含 `D_H/cost_J/ell_hat/Vmat/exc_sign`。

**資料來源 / 流向**：
- 產生：`matlab/.../Calibration_using_FEM_modeling/Hall_sensor_base_fix_dir/code/main/main.m`（→ `calib_sensor_d.mat`）、
  `…/Hall_sensor_base_no_fix_dir/code/main/main.m`（→ `calib_sensor_d_no_fix_dir.mat`）；
  legacy 等價實作 `fixl_fit/code/scripts/calib_fem.m` 也寫 `calib_sensor_d.mat`。
- 消費：`bias_fit/compare_models_nrmse.m`（讀 `calib_sensor_d.mat`）、
  `Hall_sensor_base_no_fix_dir/code/plot/verify_*`（讀 `calib_sensor_d_no_fix_dir.mat`）。

**相關**：見 [../README.md](../README.md)（charge_fit 總覽）與 [../../../../CLAUDE.md](../../../../CLAUDE.md)。
