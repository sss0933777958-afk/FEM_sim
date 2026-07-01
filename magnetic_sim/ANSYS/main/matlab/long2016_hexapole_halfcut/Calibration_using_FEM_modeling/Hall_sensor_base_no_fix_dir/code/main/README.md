# …/Hall_sensor_base_no_fix_dir/code/main/ — 統一校正主程式（18 參數 bias 版）

**用途**：`main.m` — 本包**唯一** driver（已統一，取代舊 `main_Dmatrix`）。18-param bias（電荷離軸
`ℓ̂·(Pc_base+E(ê))`）模型，求 6×6 Hall-sensor 校正並輸出 **V、D̄、^Bĝ_V、Ĥ_V**（論文 notation）。與
`../../Hall_sensor_base_fix_dir/code/main/main.m`（single-parameter 在軸）並列，**唯一差別 = ℓ̂-fit**
（本包多 17 個 ê 離軸自由度，`fit_bias`/lsqnonlin）。config 在頂部
（`VARIANT='gap200um_mueq'`、`R_select=150 µm`、`I=1 A`、`S_hall=130 V/T`、`n_uniform=10000`）。

> 論文↔code 變數：`G=Dv`（=D^v）、`Ĥ_V=H_V`（舊 `Dmat`）、`^Bĝ_V=ghat_V_B`（舊 `g_V`）、`D̄=D_bar`。
> **單位（Unit Reference Sheet，原生）**：b=**mT**、V=**mV**、ℓ̂=**µm**、^Bĝ_V=**mT/mV**、Ĥ_V=**mT/mV**、G=**mT**。
> ⚠ 實作：`fit_bias` 在 **SI 公尺**（well-scaled；µm 會讓 18-param 最佳化發散）→ ℓ̂ 存/印前 ×1e6 成 µm；B 匯入即 mT；V 直接 mV。

**流程（五步）**：
1. **載 6-coil FEM**（`load_coils_actuator`，`gap200um_mueq`）→ `select_ball` R≤150 µm 球 → all-source（literal flip-sink：只翻下極 P1/P3/P6）。
2. **fit {ℓ̂, ê(17)}**：`fit_bias`（lsqnonlin，profile Dv）→ `make_Pc` 離軸電荷位置。
3. **profile 電荷** `G = D^v = M⁻¹(Aᵀ·Bstack)`（每激發一欄）。
4. **抽 sensor 電壓 V**：`build_sensor_geometry`（**下極 −β 底錐面修正位置**）→ `extract_Vmat_interp`（**Ø0.3 mm × 0.1 mm 圓柱、均勻內插 10k 點平均**、graded sensor-local CSV、all-source）。
5. **解（論文 step 9）**：`Ĥ_V = G·Vᵀ(VVᵀ)⁻¹` → `D̄ = Ĥ_V·5/(6·h₁₁)` → **`^Bĝ_V = (6/5)·Ĥ_V(1,1)`**（Ĥ_V = ^Bĝ_V·D̄）。

**輸出**：
- console：先印 **V**（驗證）再印 D̄ / ^Bĝ_V / Ĥ_V + 摘要（ℓ̂、‖ê‖、region err、recon）。
- PDF `../../results/D_gap200um_mueq.pdf`：V、D̄、^Bĝ_V、Ĥ_V（results/ 只留 PDF）。
- `.mat` `../../data/calib_D_gap200um_mueq.mat`（**.mat field 名沿用** `Dmat/D_bar/g_V/Dv_p/Vmat_p/ell_hat/e_hat/E36/...`；ê 保留供下游）。

**預期數值（gap200um_mueq, R=150 µm；原生單位）**：ℓ̂≈**856.66 µm**（18-param bias）、region err≈0.46%、recon~1e-16、^Bĝ_V≈**7.12e-3 mT/mV**、V 對角≈**1000 mV**、G 對角全正。

**命名 / 慣例**：單一主程式 → `code/main/main.m`；模型數學一律在 `../function/`；不重跑 fminunc（fit_bias 用 lsqnonlin）。

**相關**：見上層 `../README.md`、`../function/README.md`。
