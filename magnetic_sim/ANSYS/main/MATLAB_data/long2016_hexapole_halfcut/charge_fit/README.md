# MATLAB_data/long2016_hexapole_halfcut/charge_fit — 點電荷 K_I 模型擬合 / Hall-sensor 校正 / 驗證

**用途**：存放下極半切六極的**點電荷（point-charge）模型**成果：擬合 actuator 矩陣 K̂_I、長度 ℓ̂、增益 ĝ_B，校正 Hall-sensor 常數 d，以及電流組合的閉環驗證。

**內容**（全部歸進子夾，本層只有子夾 + 本 README，不放散檔）：
- `calibration/` — 校正與全域/聯合擬合結果（`calib_bias`、`calibration_final`、`fit_KI_full`、`joint_6coil_40um_fit`）。
- `fit_KI_ball/` — 各取樣半徑的 K_I 球擬合（`fit_KI_R040..R500` 等）。
- `fitting_trend/` — R 掃描趨勢與最佳 R* 選取。
- `fitting_d/` — Hall-sensor per-pole 常數 d（`calib_sensor_d`、`calib_sensor_d_no_fix_dir`）。
- `validation/` — 電流組合閉環驗證（`validate_combos_R150`）。

各子夾有自己的 README。

**資料來源 / 流向**：由 `matlab/long2016_hexapole_halfcut/`（Calibration using FEM modeling 的 `fix_l`/`no_fix_l` main 程式 + charge-fit scripts）讀 6 顆 coil 的 **1A** FEM `.dat`（注意：擬合電流必須對齊 FEM 激發電流 1A，見 `fit-current-matches-sim` 規則）產生。resolver = `matlab_path('long2016_hexapole_halfcut','charge_fit')`。

**相關**：見 [../README.md](../README.md)（本 model 功能總覽）、[../../README.md](../../README.md)（MATLAB_data 總覽）與 [../../../CLAUDE.md](../../../CLAUDE.md)。
