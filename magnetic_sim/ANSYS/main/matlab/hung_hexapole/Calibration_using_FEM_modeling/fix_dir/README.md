# …/Calibration using FEM modeling/fix_l/ — 固定-ℓ 點電荷模型校正（fix-ℓ，**hung**）

**用途**：hung hexapole「固定-ℓ 點電荷模型」校正的乾淨單一主程式交付。電荷固定在磁極軸 `pc_i = ℓ·d̂_i`（無 bias），用 `lsqnonlin` fit `{K̄ (6×6，K̄(1,1)=5/6 固定), ℓ, ^Bĝ_I}`（論文 notation；code/.mat field 名沿用 `Khat`/`gB`）。**hung fit 結果**（R=150 µm、gap_200um）：ℓ̂≈**880.9 µm**, ^Bĝ_I≈**6.94 mT/A**, 相對 RMS 場誤差 ≈**2.42%**（C_mean=1348 (mT/A)³、κ_mean=0.829）。

**內容**：`code/main/main.m`（主程式，config 在頂部）、`code/main_function/`（載場/擬合/select_ball）、`code/function/`（模型數學輔助函式）、`code/plot/`（`plot_gain_iso_index.m` 直方圖、`plot_ref_planes_3d.m` 參考面、`plot_svd_heatmaps_2d.m` 三 actuator 面極座標熱圖）、`results/`（auto-gen `.tex`）、`data/`（`.mat` 成果，規則#2）、`figures/`（圖）。

**資料來源 / 流向**：讀 `ANSYS_data/hung_hexapole/data/coilN/gap_200um/` 6-coil FEM 場（1 A、hung mt_constants SPH_OFST=0、恆等 coil→pole 映射）→ `load_coils`（filter iron、all-source flip-sink、轉 actuator frame）→ `select_ball` 取 R 內節點 → `fit_KI_fixl` → `region_field_err` → 存解 `data/fit_fixl_R<RRR>um_gap_200um.mat`（`ell`/`gB`/`Khat`/`J`/`errpct`）。
> ⚠ 本組已 hung 化（原為 long2016 逐檔複製）；同層 `no_fix_dir/`、`Hall_sensor_base_*/` **仍是 long2016 複製、未 hung 化**。

**命名 / 慣例**：單一主程式組 → `code/main/main.m`；I_actual=1 A 對齊 FEM 激發；電荷模型全 source。

**相關**：見上層 `../README.md`、`../../../../CLAUDE.md`。
