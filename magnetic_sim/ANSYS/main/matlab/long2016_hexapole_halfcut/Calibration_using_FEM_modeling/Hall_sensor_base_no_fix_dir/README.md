# …/Calibration_using_FEM_modeling/Hall_sensor_base_no_fix_dir/ — Hall-sensor d（18 參數 bias 版，一條龍）

**用途**：求 Long Fei 半切六極 Hall-sensor 模型 per-pole 常數 `d`，但 page-1 校正用 **18 參數 bias 模型**（`no_fix_dir`，actuator 框、電荷可離軸 `Pc_18`）而非 fix-ℓ。與姊妹 `Hall_sensor_base_fix_dir`（fix-ℓ 版）對稱。模型 `b_ij = g_H·S_i V_j d`（`g_H = 1/(4πℓ̂²)`），閉式解 `d`（含/無增益兩版），並依 derivation.pdf 求 `ĝ_F`/`K̂_H`/`R_a`。sensor 模型在 R≤150 µm 的相對 RMSE ≈ **8.4%**（比 fix-ℓ 版的 ~15.5% 低，因 18 參數場擬合更準）。

**內容**：`code/main/main.m`（求 d 一條龍 driver）、`code/function/`（模型/IO 輔助函式）、`code/plot/`（兩支 d-form vs Q-form 驗證）、`results/`（auto-gen `.tex`：d_v2 / d_final / KH_v2 / KH_final，直接放）、`figures/`（圖，目前空）。

**資料來源 / 流向**：**載入** `MATLAB_data/.../charge_fit/calib_bias.mat`（no_fix_dir 校正結果：`ell_hat≈0.857 mm, Pc_18, R, F`）→ 重載 6-coil FEM 場（1 A）旋進 actuator 框 → `build_A` 建 `M, c` → page-2 建 sensor 幾何、抽 all-source `Vmat`、解 `d`、算殘差 → 存 `calib_sensor_d_no_fix_dir.mat`（**不蓋** fix-ℓ 版的 `calib_sensor_d.mat`）→ 輸出 4 個 `.tex` 到 `results/`。

**命名 / 慣例**：架構比照 `../Hall_sensor_base_fix_dir/`（`code/{main,function,plot}` + `results/` + `figures/`，results/figures 不再多包子夾）；I=1 A 對齊 FEM 激發；sensor 號=物理 signed B·n+、all-source（翻下極 P1/P3/P6）；charge model 在 actuator 框、`Pc_18` 離軸。page-1 的 18 參數擬合本身在 `../no_fix_dir/` 與 `../../bias_fit/calib_fem_bias.m`（本專案只載入其結果，不重跑）。

**相關**：見上層 `../README.md`、姊妹 `../Hall_sensor_base_fix_dir/README.md`、`../no_fix_dir/README.md`、`../../../../CLAUDE.md`。
