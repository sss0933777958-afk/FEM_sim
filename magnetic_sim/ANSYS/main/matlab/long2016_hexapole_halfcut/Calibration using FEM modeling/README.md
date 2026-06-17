# matlab/long2016_hexapole_halfcut/Calibration using FEM modeling/ — FEM 點電荷校正（charge calibration）功能組

**用途**：用 FEM 場校正 Long Fei 半切六極的「點電荷模型」（point-charge model）。同一個物理模型、兩種建模方式各放一個乾淨的單一主程式（single main program）交付夾：
- `fix_l/` — document「固定-ℓ 點電荷模型」（fix-ℓ）；電荷固定在磁極軸 `pc_i = ℓ·d̂_i`，fit `{K̂_I, ℓ, gB}`。R=150 → ℓ≈0.856 mm, gB≈8.43e-3。
- `no_fix_l/` — 「18-param bias 模型」；電荷可離軸（actuator frame，含 1×17 bias e），每解 profile 出 6 個電荷量 g_j。R=150 → ℓ≈0.857 mm, gB≈9.50e-3。

**內容**：每個子夾 = `code/{main,function,plot}/` + `results/` + `figures/`。`code/main/main.m` 是唯一主程式（頂部有 config：MODE single/sweep、R 範圍、I_actual=1、ball 取樣）；`code/function/` 放模型數學輔助函式；`code/plot/` 放該主程式的繪圖；`results/<fix_l|no_fix_l>/` 收純結果 `.tex`。

**資料來源 / 流向**：讀 `ANSYS_data/long2016_hexapole_halfcut/` 的 6-coil FEM 場（`.dat`，1 A 激發）→ `select_ball` 取 WP 半徑 R 內節點 → `lsqnonlin` 擬合 → `region_field_err` 算相對 RMS 誤差 → 輸出純結果 `.tex`（檔名 `fit_ball_R<R>um_<I>A.tex`）到 `results/<fix_l|no_fix_l>/`。

**命名 / 慣例**：單一主程式組（不是多腳本組），故用 `code/main/main.m`；model-first；電荷模型全 source（all-source flip：上極 P2/P4/P5）；I_actual = 1 A 對齊 FEM 激發。

**相關**：見上層 `../README.md`、`../../../CLAUDE.md`（matlab schema、繪圖腳本規則、fit-current 規則）。
