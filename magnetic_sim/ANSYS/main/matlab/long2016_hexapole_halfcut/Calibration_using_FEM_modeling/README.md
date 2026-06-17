# matlab/long2016_hexapole_halfcut/Calibration_using_FEM_modeling/ — FEM 點電荷校正 + Hall-sensor d 功能組

**用途**：用 FEM 場校正 Long Fei 半切六極的「點電荷模型」（point-charge model）與其 Hall-sensor 延伸。每個子夾各放一個乾淨的單一主程式（single main program）交付夾：
- `fix_dir/` — document「固定方向（fix-ℓ）點電荷模型」；電荷固定在磁極軸 `pc_i = ℓ·d̂_i`，fit `{K̂_I, ℓ, gB}`。R=150 → ℓ≈0.856 mm, gB≈8.43e-3。
- `no_fix_dir/` — 「18-param bias 模型」（電荷可離軸 / 無固定方向）；actuator frame，含 1×17 bias e，每解 profile 出 6 個電荷量 g_j。R=150 → ℓ≈0.857 mm, gB≈9.50e-3。
- `Hall_sensor_base_fix_dir/` — Hall-sensor per-pole 常數 `d`，**基於 fix_dir**：把電荷綁到 sensor 讀數 `b_ij = g_H·S_i V_j d`，閉式解 d（含/無增益兩版）+ 求 ĝ_F/K̂_H/R_a + 驗證 + LaTeX。R=150 → 相對 RMSE ≈15.5%，R_a≈7.79e8。
- `Hall_sensor_base_no_fix_dir/` — 同上但 **基於 no_fix_dir（18-param bias，actuator 框、Pc_18 離軸）**；載入 `calib_bias.mat` 不重跑 18 維 fit。R=150 → 相對 RMSE ≈8.4%，R_a≈6.85e8。

**內容**：每個子夾 = `code/{main,function,plot}/` + `results/` + `figures/`。`code/main/main.m` 是唯一主程式（頂部有 config）；`code/function/` 放模型數學輔助函式；`code/plot/` 放繪圖/驗證；`results/` 直接收純結果 `.tex`（**不再多包子夾**）；`figures/` 放圖。

**資料來源 / 流向**：讀 `ANSYS_data/long2016_hexapole_halfcut/` 的 6-coil FEM 場（`.dat`，1 A 激發）→ 取 WP 半徑 R 內節點 → 擬合 / 求 d → 輸出純結果 `.tex` 到各子夾 `results/`。Hall_sensor_base_* 另把 d 結果存 `MATLAB_data/.../charge_fit/calib_sensor_d{,_no_fix_dir}.mat`。

**命名 / 慣例**：單一主程式組（不是多腳本組），故用 `code/main/main.m`；model-first；電荷模型全 source（all-source flip）；I_actual = 1 A 對齊 FEM 激發；`results/`、`figures/` 不再多包一層子夾。

**相關**：見上層 `../README.md`、`../../../CLAUDE.md`（matlab schema、繪圖腳本規則、fit-current 規則）。
