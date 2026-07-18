# …/Calibration using FEM modeling/no_fix_l/ — 18-param bias 點電荷模型校正（no-fix-ℓ）

**用途**：「18-param bias 模型」的乾淨單一主程式交付。電荷可離軸：`pc = ℓ·(Pc_base + E(ê))`（actuator frame，18 params = ℓ + 1×17 bias ê，e6z constrained）；每解 LS profile 出 6 個電荷量 g_j（=G=D^v 的欄），再 gauge（K̄(1,1)=5/6）得 ^Bĝ_I、K̄（論文 step 8；中間 Ĥ_I = G·Fᵀ(FFᵀ)⁻¹）。**單位 Unit Sheet**（b=mT、V=mV、ℓ̂=µm、^Bĝ_I=mT/A；擬合在 SI 公尺、ℓ̂ 輸出 ×1e6）：R=150 µm → ℓ̂≈**857 µm**, ^Bĝ_I≈**8.32 mT/A**（err 0.46%）。

**內容**：`code/main/main.m`（主程式，config 在頂部）、`code/function/`（模型數學輔助函式）、`code/plot/`（該主程式繪圖）、`results/`（auto-gen `.tex`）、`figures/`（圖）。

**資料來源 / 流向**：讀 `ANSYS_data/long2016_hexapole_halfcut/`（用 `../../common/ansys_path`）6-coil FEM 場（1 A、轉 actuator frame）→ `select_ball` → `fit_bias` → `gauge_KI` → `region_field_err` → `write_KbarI_tex` 出 `results/no_fix_l/fit_ball_R<R>um_<I>A.tex`。

**命名 / 慣例**：單一主程式組 → `code/main/main.m`；I_actual=1 A 對齊 FEM 激發。

**hung 專屬號誌 / 編號（2026-07-17 修正，別再抄 long2016）**：
- **`apdl_to_paper_idx = [1 2 3 4 5 6]`（identity）**：hung 的 deck 照 paper 序建（`mesh/_mesh_R700.txt` 的 `fc_ra(k)` 就是 P k 的方位）。long2016 的 `[1,3,6,5,2,4]` 是**它自己 deck 的建構順序**，套到 hung 會讓 K̄_I 欄位錯位。
- **不做 per-pole 翻號**：hung raw FEM = **六極全 sink**（deck 用 `+TURNS`），`load_coils_actuator.m:46` 的全域 negate 已讓六顆全 source。已移除從 long2016 抄來的 `coil_sign=[1 -1 1 -1 -1 1]`。
- **驗證判準**：K̄_I **對角占優 + 對角全正 + off-diag 全負 + 每列和≈0**、`K̄(1,1)=5/6`。
- 規則全文：`magnetic_sim/.claude/rules/pole-coil-numbering.md`、`charge-model-source-convention.md`。

**相關**：見上層 `../README.md`、`../../../../CLAUDE.md`。
