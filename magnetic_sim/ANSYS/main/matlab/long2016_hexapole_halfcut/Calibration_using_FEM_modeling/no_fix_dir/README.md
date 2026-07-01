# …/Calibration using FEM modeling/no_fix_l/ — 18-param bias 點電荷模型校正（no-fix-ℓ）

**用途**：「18-param bias 模型」的乾淨單一主程式交付。電荷可離軸：`pc = ℓ·(Pc_base + E(ê))`（actuator frame，18 params = ℓ + 1×17 bias ê，e6z constrained）；每解 LS profile 出 6 個電荷量 g_j（=G=D^v 的欄），再 gauge（K̄(1,1)=5/6）得 ^Bĝ_I、K̄（論文 step 8；中間 Ĥ_I = G·Fᵀ(FFᵀ)⁻¹）。**單位 Unit Sheet**（b=mT、V=mV、ℓ̂=µm、^Bĝ_I=mT/A；擬合在 SI 公尺、ℓ̂ 輸出 ×1e6）：R=150 µm → ℓ̂≈**857 µm**, ^Bĝ_I≈**8.32 mT/A**（err 0.46%）。

**內容**：`code/main/main.m`（主程式，config 在頂部）、`code/function/`（模型數學輔助函式）、`code/plot/`（該主程式繪圖）、`results/`（auto-gen `.tex`）、`figures/`（圖）。

**資料來源 / 流向**：讀 `ANSYS_data/long2016_hexapole_halfcut/`（用 `../../common/ansys_path`）6-coil FEM 場（1 A、轉 actuator frame）→ `select_ball` → `fit_bias` → `gauge_KI` → `region_field_err` → `write_KbarI_tex` 出 `results/no_fix_l/fit_ball_R<R>um_<I>A.tex`。

**命名 / 慣例**：單一主程式組 → `code/main/main.m`；I_actual=1 A 對齊 FEM 激發；coil_sign=[1 -1 1 -1 -1 1] 全 source 顯示（翻上極 P2/P4/P5）。

**相關**：見上層 `../README.md`、`../../../../CLAUDE.md`。
