# …/Calibration using FEM modeling/Hall_sensor_base_fix_dir/ — Hall-sensor 模型求每極常數 d

**用途**：Long Fei 半切六極 **Hall-sensor 模型每極常數 `d` 的求解**（fix-ℓ：電荷鎖在極軸方向 d̂）。模型 `b_ij = g_H·S_i·V_j·d`（`g_H = 1/(4πℓ̂²)`），閉式解 `d`（minJ）。本包已收斂成單一條指定流程，不含後處理/LaTeX/驗證腳本。

**指定流程（四步）**：
1. **拿 ℓ̂**：從 fix_dir 校正結果載入 —— `MATLAB_data/long2016_hexapole_halfcut/charge_fit/fit_KI_ball/fit_KI_R<RRR>.mat` 的 `ell`（R=R_select，預設 150 µm → ℓ̂≈0.852 mm）。**不再用 fminbnd 自己 fit**。
2. **抽電壓 V=S·B（真實 FEM 節點）**：`build_sensor_geometry` 給 6 顆 sensor 中心/法線 n+；`extract_Vmat` 沿法線 n 開半徑 0.15 mm 圓柱選**真實節點**對 B·n 平均（**不內插**）；圓柱內無節點時取「in-plane 落在盤內、沿 n 最近」的節點（遠場 mesh 粗，實測多為每 sensor 命中 1 點）。含 all-source 翻號（翻下極激發 P1/P3/P6）。
3. **建模型**：殘差 ε = b_model − b_FEM、cost `J = Σ‖ε‖²`。
4. **minJ 閉式解 d**（`solve_d`）。

**內容**：`code/main/main.m`（driver，config 在頂部）、`code/function/`（5 支：`build_S` / `build_sensor_geometry` / `extract_Vmat` / `solve_d` / `sensor_residual`）。求解結果存 `MATLAB_data/long2016_hexapole_halfcut/charge_fit/fitting_d/calib_sensor_d.mat`（含 `d/gH/Vmat/exc_sign/ell_hat/J`）。

**資料來源 / 流向**：讀 `ANSYS_data/long2016_hexapole_halfcut/data/coil1..6` 6-coil FEM 場（1 A）→ 載 ℓ̂ → R≤R_select 球內 air 節點建 M/c → sensor 幾何 + 真實節點抽 `Vmat` → `solve_d` → `sensor_residual`（cost J）→ 存 `.mat`。

**預期數值（R=150 µm）**：ℓ̂≈0.852 mm、cost J≈0.142 T²（電壓改 real-node 單點抽樣）。

**命名 / 慣例**：I=1 A 對齊 FEM 激發；sensor 號=物理 signed B·n+、all-source（翻下極激發 P1/P3/P6）。

**相關**：見上層 `../README.md`、`../../../../CLAUDE.md`。
