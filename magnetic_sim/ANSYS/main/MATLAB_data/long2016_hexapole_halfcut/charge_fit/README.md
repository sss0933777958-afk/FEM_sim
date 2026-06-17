# MATLAB_data/long2016_hexapole_halfcut/charge_fit — 點電荷 K_I 模型擬合 / Hall-sensor 校正 / 驗證

**用途**：存放下極半切六極的**點電荷（point-charge）模型**成果：擬合 actuator 矩陣 K̂_I、長度 ℓ̂、增益 ĝ_B，校正 Hall-sensor 常數 d，以及電流組合的閉環驗證。

**內容**（代表檔）：
- `fit_KI_full.mat`：document 版電荷模型擬合（K̂_I、ℓ̂、ĝ_B、R_a）。
- `calibration_final.mat`：最終校正（選定 R*=150 µm，ℓ̂≈0.856 mm、gB≈8.43e-3、‖K̂‖≈2.436）。
- `calib_sensor_d.mat` / `calib_bias.mat`：Hall-sensor per-pole 常數 d（含 g_H、殘差）/ bias 電荷版校正。
- `joint_6coil_40um_fit.mat`：R=40 µm 球的 6-coil 聯合擬合。
- `validate_combos.mat` / `validate_combos_R150.mat`：電流組合閉環驗證（R150 為選定半徑）。
- 子夾 `fit_KI_ball/`（各取樣半徑擬合）、`fitting_trend/`（R 掃描趨勢與選取）—各自有 README。

**資料來源 / 流向**：由 `matlab/long2016_hexapole_halfcut/`（Calibration using FEM modeling 的 `fix_l`/`no_fix_l` main 程式 + charge-fit scripts）讀 6 顆 coil 的 **1A** FEM `.dat`（注意：擬合電流必須對齊 FEM 激發電流 1A，見 `fit-current-matches-sim` 規則）產生。resolver = `matlab_path('long2016_hexapole_halfcut','charge_fit')`。

**相關**：見 [../README.md](../README.md)（本 model 功能總覽）、[../../README.md](../../README.md)（MATLAB_data 總覽）與 [../../../CLAUDE.md](../../../CLAUDE.md)。
