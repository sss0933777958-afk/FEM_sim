# …/sensor_d/code/main/ — 求 d 主程式

**用途**：`main.m` — Hall-sensor per-pole `d` 的一條龍 driver。config 在頂部（R_select=150 µm、I=1 A、S_hall=130 V/T）。

**流程**：
1. 載 6-coil FEM 場（`wp` dataset）→ 取 R 球內節點 → `cost_J` + `fminbnd` 擬合 ℓ̂ → 在 ℓ̂ 建 M、c（page-1）。
2. `build_sensor_geometry` → `extract_Vmat`（all-source）→ `solve_d` → `sensor_residual`（page-2）。
3. 存 `calib_sensor_d.mat`（原 `MATLAB_data/.../charge_fit/` 路徑）。
4. `write_d_tex`（d_v2/d_final）+ `compute_KH`+`write_KH_tex`（KH_v2/KH_final）→ `../../results/sensor_d/`。

**預期數值**（R=150 µm）：ℓ̂≈0.856 mm、sensor 模型相對 RMSE ≈15.5%；與既有 `calib_sensor_d.mat` 一致。

**命名 / 慣例**：單一主程式組 → `code/main/main.m`；模型數學一律在 `../function/`。

**相關**：見上層 `../README.md`、`../function/README.md`。
