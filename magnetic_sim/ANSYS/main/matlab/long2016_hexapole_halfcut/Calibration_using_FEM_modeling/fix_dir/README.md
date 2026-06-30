# …/Calibration using FEM modeling/fix_l/ — 固定-ℓ 點電荷模型校正（fix-ℓ）

**用途**：document「固定-ℓ 點電荷模型」的乾淨單一主程式交付。電荷固定在磁極軸 `pc_i = ℓ·d̂_i`（無 bias），用 `lsqnonlin` fit `{K̂_I^FEM (6×6，K̂(1,1)=5/6 固定), ℓ, gB}`。R=150 µm → ℓ≈0.856 mm, gB≈8.43e-3, 相對 RMS 場誤差 ≈3.18%。

**內容**：`code/main/main.m`（主程式，config 在頂部）、`code/function/`（模型數學輔助函式）、`code/plot/`（該主程式繪圖）、`results/`（auto-gen `.tex`）、`data/`（`.mat` 成果，規則#2）、`figures/`（圖）。

**資料來源 / 流向**：讀 `ANSYS_data/long2016_hexapole_halfcut/` 6-coil FEM 場（1 A）→ `select_ball` 取 R 內節點 → `fit_KI_fixl` → `region_field_err` → `write_KI_tex` 輸出 `results/fit_ball_R<R>um_<I>A.tex`；**並存解** `data/fit_fixl_R<RRR>um.mat`（本組 `data/`，規則#2；`ell`/`gB`/`Khat`/`J`/`errpct`，供 `Hall_sensor_base_fix_dir/decouple` 載入 ℓ̂）。

**命名 / 慣例**：單一主程式組 → `code/main/main.m`；I_actual=1 A 對齊 FEM 激發；電荷模型全 source。

**相關**：見上層 `../README.md`、`../../../../CLAUDE.md`。
