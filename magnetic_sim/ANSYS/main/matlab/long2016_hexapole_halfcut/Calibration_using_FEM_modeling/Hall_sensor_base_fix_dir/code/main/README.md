# …/Hall_sensor_base_fix_dir/code/main/ — 統一校正主程式

**用途**：`main.m` — 本包**唯一** driver（已統一，取代舊 `main_Dmatrix`/`main_Vmat`/`main_interp`）。
single-parameter（電荷在軸 d̂、無 bias）模型，求 6×6 Hall-sensor 校正並輸出 **V、D̄、^Bĝ_V、Ĥ_V**（論文 notation）。
config 在頂部（`VARIANT='gap200um_mueq'`、`R_select=150 µm`、`I=1 A`、`S_hall=130 V/T`、`n_uniform=10000`）。

> 論文↔code 變數：`G=Dv`（profiled charges=D^v）、`Ĥ_V=H_V`（舊 `Dmat`）、`^Bĝ_V=ghat_V_B`（舊 `g_V`）、`D̄=D_bar`。
> **單位（Unit Reference Sheet，原生）**：b=**mT**、V=**mV**、ℓ̂=**µm**、^Bĝ_V=**mT/mV**、Ĥ_V=**mT/mV**、G=**mT**。
> ⚠ 實作：擬合在 **SI 公尺**（well-scaled 最佳化）→ ℓ̂ 存/印前 ×1e6 成 µm；B 於 .dat 匯入即 ×1e3 成 mT；V=S_hall·B 直接 mV。

**流程（五步）**：
1. **載 6-coil FEM**（`load_coils_actuator`，`gap200um_mueq`、actuator 框）→ `select_ball` 選 R≤150 µm 球內 air 節點 → all-source（literal flip-sink：只翻下極 P1/P3/P6）。
2. **fit ℓ̂**（single-parameter 在軸 profiled field cost，`fminbnd`；Dv 已 profile 掉）。
3. **profile 電荷** `G = D^v = M⁻¹(Aᵀ·Bstack)`（每激發一欄）。
4. **抽 sensor 電壓 V**：`build_sensor_geometry`（**下極 −β 底錐面修正位置**）→ `extract_Vmat_interp`（**Ø0.3 mm × 0.1 mm 圓柱、均勻內插 10k 點平均**、graded sensor-local CSV、all-source）。
5. **解（論文 step 9）**：`Ĥ_V = G·Vᵀ(VVᵀ)⁻¹` → `D̄ = Ĥ_V·5/(6·h₁₁)`（gauge D̄₁₁=5/6）→ **`^Bĝ_V = (6/5)·Ĥ_V(1,1)`**（電壓側增益 T/V，Ĥ_V = ^Bĝ_V·D̄）。

**輸出**：
- console：先印 **V**（驗證）再印 D̄ / ^Bĝ_V / Ĥ_V + 摘要（ℓ̂、region err、recon）。
- PDF `../../results/D_gap200um_mueq.pdf`：V、D̄、^Bĝ_V、Ĥ_V（results/ 只留 PDF）。
- `.mat` `../../data/calib_D_gap200um_mueq.mat`（**.mat field 名沿用** `Dmat/D_bar/g_V/Dv_p/Vmat_p/ell_hat/...`，與下游 plot 腳本相容）。

**預期數值（gap200um_mueq, R=150 µm；原生單位）**：ℓ̂≈**867.45 µm**（on-axis）、recon `‖Ĥ_V·V−G‖/‖G‖`~1e-16、^Bĝ_V≈**6.14e-3 mT/mV**、V 對角≈**1000 mV**、region err 3.19%、G 對角全正。

**命名 / 慣例**：單一主程式 → `code/main/main.m`；模型數學一律在 `../function/`。

**相關**：見上層 `../README.md`、`../function/README.md`。
