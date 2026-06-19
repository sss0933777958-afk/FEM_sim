# MATLAB_data/long2016_hexapole_halfcut/charge_fit/calibration — 點電荷模型校正 / 全域擬合結果

**用途**：存放下極半切六極點電荷模型的**校正與全域/聯合擬合**結果（與同層 `fit_KI_ball/` 各半徑掃描、
`fitting_trend/` 趨勢、`fitting_d/` Hall-sensor d 並列）。

**內容**：
- `calib_bias.mat`：bias 電荷版校正（no_fix_dir：ℓ̂、Pc_18、R、F；no_fix_dir 求 d 的輸入）。
- `calibration_final.mat`：最終校正（選定 R*=150 µm，ℓ̂≈0.856 mm、gB≈8.43e-3、‖K̂‖≈2.436）。**無 in-repo producer**（legacy 結果）。
- `fit_KI_full.mat`：document 版電荷模型全節點擬合（K̂_I、ℓ̂、ĝ_B、R_a）。
- `joint_6coil_40um_fit.mat`：R=40 µm 球的 6-coil 聯合擬合。

**資料來源 / 流向**：
- 產生：`bias_fit/.../calib_fem_bias.m`（→ `calib_bias.mat`）、`fixl_fit/.../fit_KI_full.m`（→ `fit_KI_full.mat`）、
  `fixl_fit/.../test_joint_6coil_fit_40um.m`（→ `joint_6coil_40um_fit.mat`）。
- 消費：`Hall_sensor_base_no_fix_dir` + `no_fix_dir` plot（讀 `calib_bias.mat`）、
  `validation/eval_validate_combos.m`（讀 `calibration_final.mat`）、`validation/validate_KI.m` + `doc/.../gen_KI_latex.m`（讀 `fit_KI_full.mat`）、
  `field_viz/plot_P1_fig25c_charge.m` + `doc/.../gen_fit_J_cube40_latex.m`（讀 `joint_6coil_40um_fit.mat`）。

**相關**：見 [../README.md](../README.md)（charge_fit 總覽）與 [../../../../CLAUDE.md](../../../../CLAUDE.md)。
