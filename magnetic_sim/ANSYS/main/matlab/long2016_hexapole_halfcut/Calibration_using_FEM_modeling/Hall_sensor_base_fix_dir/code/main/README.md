# …/Hall_sensor_base_fix_dir/code/main/ — 求 d 主程式

**用途**：`main.m` — Hall-sensor 每極 `d` 的指定流程 driver。config 在頂部（R_select=150 µm、I=1 A、S_hall=130 V/T）。

**流程（四步）**：
1. **拿 ℓ̂**：載入 `MATLAB_data/.../charge_fit/fit_KI_ball/fit_KI_R<RRR>.mat` 的 `ell`（不再 fminbnd）。
2. **抽電壓 V=S·B**：載 6-coil FEM（`wp` 建殘差場 B / `all` 抽 sensor 電壓）→ `build_sensor_geometry` → `extract_Vmat`（真實節點、沿 n 圓柱、all-source）。
3. **建模型/殘差**：在 ℓ̂ 建 M、c → `solve_d` 解 d → `sensor_residual` 回 cost `J = Σ‖ε‖²`。
4. **存解**：`calib_sensor_d.mat`（`MATLAB_data/.../charge_fit/fitting_d/`，含 d/gH/Vmat/exc_sign/ell_hat/J）。

**預期數值（R=150 µm）**：ℓ̂≈0.852 mm、cost J≈0.142 T²。

**命名 / 慣例**：單一主程式組 → `code/main/main.m`；模型數學一律在 `../function/`。

**相關**：見上層 `../README.md`、`../function/README.md`。
