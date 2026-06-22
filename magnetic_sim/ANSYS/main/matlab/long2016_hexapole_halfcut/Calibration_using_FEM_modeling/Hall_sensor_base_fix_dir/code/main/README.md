# …/Hall_sensor_base_fix_dir/code/main/ — 求 d 主程式

**用途**：`main.m` — Hall-sensor 每極 `d` 的指定流程 driver。config 在頂部（R_select=150 µm、I=1 A、S_hall=130 V/T）。

**流程（四步）**：
1. **拿 ℓ̂**：載入 fix_dir fit_KI_fixl 解出的 `MATLAB_data/.../charge_fit/calibration/fit_fixl_R<RRR>um.mat` 的 `ell`（≈0.856；不再 fminbnd、不用 fit_KI_ball sweep 版）。
2. **抽電壓 V=S·B**：載 6-coil FEM（`wp` 建殘差場 B / `all` 抽 sensor 電壓）→ `build_sensor_geometry` → `extract_Vmat`（真實節點、沿 n 圓柱、all-source）。
3. **解 d / 殘差**：`solve_d`（內外層雙重加總、閉式解 d，**無 g_H**：`d = (Σ_j V_j(Σ_i SᵀS)V_j)⁻¹(Σ_j V_j Σ_i Sᵀb)`）→ `sensor_residual` 回 cost `J = Σ‖ε‖²`（模型 b=S·V·d）。
4. **存解**：`calib_sensor_d.mat`（`MATLAB_data/.../charge_fit/fitting_d/`，含 d/Vmat/exc_sign/ell_hat/J，**無 gH**）。

**預期數值（R=150 µm）**：ℓ̂≈0.856 mm（fix_dir fit_KI_fixl）、cost J≈0.142 T²（全域座標、上極 CAD 傾角 36.59° 錐面 + 下極磨平面）、d 為 no-gain（~1e-2）。

**命名 / 慣例**：單一主程式組 → `code/main/main.m`；模型數學一律在 `../function/`。

**相關**：見上層 `../README.md`、`../function/README.md`。
