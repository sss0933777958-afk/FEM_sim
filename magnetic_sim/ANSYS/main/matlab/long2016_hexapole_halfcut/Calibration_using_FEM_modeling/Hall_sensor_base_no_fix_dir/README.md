# …/Calibration_using_FEM_modeling/Hall_sensor_base_no_fix_dir/ — Hall-sensor d（18 參數 bias 版，一條龍）

**用途**：求 Long Fei 半切六極 Hall-sensor 模型 per-pole 常數 `d`，但 page-1 校正用 **18 參數 bias 模型**（`no_fix_dir`，actuator 框、電荷可離軸 `Pc_18`）而非 fix-ℓ。與姊妹 `Hall_sensor_base_fix_dir`（fix-ℓ 版）對稱。模型 `b_ij = S_i V_j d`（**no-gain，無 g_H / K_H**），閉式解 `d`，**輸出與 fix-ℓ 版對齊：只回 `d`、`Vmat`、cost `J = Σ‖ε‖²`（無 RMSE）**。**唯一差異** = page-1 用 18 參數 bias 場（`Pc_18`、actuator 框）。

**內容**：`code/main/main.m`（求 d 一條龍 driver）、`code/function/`（模型/IO 輔助函式）、`code/plot/`（兩支 sensor 磁路箭頭圖，真實節點）、`data/`（`.mat` 成果，規則#2）、`figures/`（圖）。（已移除 g_H / K_H / LaTeX：無 `results/.tex`、無 verify 腳本，對齊 fix-ℓ 版的純 cost 輸出。）

**資料來源 / 流向**：**載入** `MATLAB_data/.../charge_fit/calibration/calib_bias.mat`（no_fix_dir 校正結果：`ell_hat≈0.857 mm, Pc_18, R, F`）→ 重載 6-coil FEM 場（1 A）旋進 actuator 框 → `build_S` 逐點建 `M, c` → page-2 建 sensor 幾何、抽 all-source `Vmat`、解 `d`、算殘差 → 存 `data/calib_sensor_d_no_fix_dir.mat`（本組 `data/`，規則#2；`d`、`Vmat`、cost `J`；**不蓋** fix-ℓ 版的 `calib_sensor_d.mat`）。`calib_bias.mat` 仍讀自 `MATLAB_data/.../charge_fit/calibration/`（外部 `bias_fit` 產，legacy）。

**命名 / 慣例**：架構比照 `../Hall_sensor_base_fix_dir/`（`code/{main,function,plot}` + `results/` + `figures/`，results/figures 不再多包子夾）；I=1 A 對齊 FEM 激發；sensor 號=物理 signed B·n+、all-source（翻下極 P1/P3/P6）；charge model 在 actuator 框、`Pc_18` 離軸。page-1 的 18 參數擬合本身在 `../no_fix_dir/` 與 `../../bias_fit/calib_fem_bias.m`（本專案只載入其結果，不重跑）。

**相關**：見上層 `../README.md`、姊妹 `../Hall_sensor_base_fix_dir/README.md`、`../no_fix_dir/README.md`、`../../../../CLAUDE.md`。
