# matlab/long2016_hexapole_halfcut/Calibration_using_FEM_modeling/ — FEM 點電荷校正 + Hall-sensor d 功能組

**用途**：用 FEM 場校正 Long Fei 半切六極的「點電荷模型」（point-charge model）與其 Hall-sensor 延伸。每個子夾各放一個乾淨的單一主程式（single main program）交付夾：
- `fix_dir/` — document「固定方向（fix-ℓ）點電荷模型」；電荷固定在磁極軸 `pc_i = ℓ·d̂_i`，fit `{K̄, ℓ, ^Bĝ_I}`（論文 notation；.mat field 名沿用 Khat/gB）。**單位 Unit Sheet**：R=150 → ℓ̂≈**867 µm**, ^Bĝ_I≈**7.33 mT/A**（err 3.19%）。
- `no_fix_dir/` — 「18-param bias 模型」（電荷可離軸 / 無固定方向）；actuator frame，含 1×17 bias e，每解 profile 出 6 個電荷量 g_j（=G=D^v 的欄）。R=150 → ℓ̂≈**857 µm**, ^Bĝ_I≈**8.32 mT/A**（err 0.46%）。
- **單位慣例（Unit Reference Sheet；全組原生）**：b=**mT**、V=**mV**、ℓ̂=**µm**、^Bĝ_I=**mT/A**、^Bĝ_V=**mT/mV**、I=A。⚠ 擬合在 **SI 公尺**（well-scaled；ℓ̂ 於輸出 ×1e6 成 µm）；B 於 .dat 匯入即 ×1e3 成 mT。
- `Hall_sensor_base_fix_dir/` — Hall-sensor per-pole 常數 `d`，**基於 fix_dir**：把電荷綁到 sensor 讀數 `b_ij = g_H·S_i V_j d`，minJ 閉式解 d。**已收斂成單一指定流程**（載 ℓ̂ → 真實節點抽 V → 解 d；不含後處理/K̂_H/LaTeX/驗證）。電壓改用**真實 FEM 節點**（沿 n 圓柱選點、不內插）。R=150 → ℓ̂≈0.852 mm、cost J≈0.142 T²。
- `Hall_sensor_base_no_fix_dir/` — 同上但 **基於 no_fix_dir（18-param bias，actuator 框、Pc_18 離軸）**；載入 `calib_bias.mat` 不重跑 18 維 fit。R=150 → 相對 RMSE ≈8.4%，R_a≈6.85e8。

**內容**：多數子夾 = `code/{main,function,plot}/` + `results/` + `figures/`。`code/main/main.m` 是唯一主程式（頂部有 config）；`code/function/` 放模型數學輔助函式；`code/plot/` 放繪圖/驗證；`results/` 直接收純結果 `.tex`；`figures/` 放圖。**例外**：`Hall_sensor_base_fix_dir/` 已收斂成單一流程，只剩 `code/{main,function}/`（無 plot/results/figures、無 LaTeX 輸出），結果只存 `.mat`。

**資料來源 / 流向**：讀 `ANSYS_data/long2016_hexapole_halfcut/` 的 6-coil FEM 場（`.dat`，1 A 激發）→ 取 WP 半徑 R 內節點 → 擬合 / 求 d → 輸出純結果 `.tex` 到各子夾 `results/`。Hall_sensor_base_* 另把 d 結果存 `MATLAB_data/.../charge_fit/fitting_d/calib_sensor_d{,_no_fix_dir}.mat`。

**命名 / 慣例**：單一主程式組（不是多腳本組），故用 `code/main/main.m`；model-first；電荷模型全 source（all-source flip）；I_actual = 1 A 對齊 FEM 激發；`results/`、`figures/` 不再多包一層子夾。

**相關**：見上層 `../README.md`、`../../../CLAUDE.md`（matlab schema、繪圖腳本規則、fit-current 規則）。
