# …/Calibration using FEM modeling/sensor_d/ — Hall-sensor 模型 per-pole 常數 d（一條龍）

**用途**：Long Fei 半切六極 **Hall-sensor 模型 per-pole 常數 `d` 的求解 → 後處理 → 驗證 → LaTeX 一條龍**。在 fix-ℓ 校正（R=150 µm → ℓ̂≈0.856 mm）基礎上，把電荷綁到 sensor 讀數 `v_jk`：模型 `b_ij = g_H·S_i V_j d`（`g_H = 1/(4πℓ̂²)`），閉式解 `d`，並依 derivation.pdf 求 `ĝ_F` / `K̂_H` / `R_a`。求出 d 含兩版：`d_final`（含增益，d~1e-8）與 `d_v2`（無增益，`d_v2 = d_final·g_H`，d~1e-3）。sensor 模型在 R≤150 µm 的相對 RMSE ≈ 15.5%（比 36-DOF 的 K̂_I fit 粗，因電荷被綁到 sensor 讀數）。

**內容**：`code/main/main.m`（求 d 一條龍 driver，config 在頂部）、`code/function/`（模型/IO 輔助函式）、`code/plot/`（兩支 d-form vs Q-form 驗證腳本）、`results/sensor_d/`（auto-gen `.tex`：d_v2 / d_final / KH_v2 / KH_final）、`figures/`（圖，目前空）。

**資料來源 / 流向**：讀 `ANSYS_data/long2016_hexapole_halfcut/` 6-coil FEM 場（1 A）→ page-1 擬合 ℓ̂、建 M/c → page-2 建 sensor 幾何（4.572 mm 沿軸 + 0.41 mm 離面、Ø0.3 mm disc 721 點面積平均）、抽 all-source `Vmat`、解 `d`、算殘差 → 存 `MATLAB_data/long2016_hexapole_halfcut/charge_fit/calib_sensor_d.mat`（同原路徑，供 `plot/verify_*` 與舊下游沿用）→ `write_d_tex` / `compute_KH`+`write_KH_tex` 輸出 4 個 `.tex` 到 `results/sensor_d/`。

**命名 / 慣例**：架構比照 `../fix_l/`（`code/{main,function,plot}` + `results/` + `figures/`）；I=1 A 對齊 FEM 激發；sensor 號=物理 signed B·n+、all-source（翻下極激發 P1/P3/P6）。求 d 的等價實作也在 `../../fixl_fit/code/scripts/calib_fem.m`（PAGE 2，保留並存）。

**相關**：見上層 `../README.md`、姊妹 `../fix_l/README.md`、`../../../../CLAUDE.md`。
