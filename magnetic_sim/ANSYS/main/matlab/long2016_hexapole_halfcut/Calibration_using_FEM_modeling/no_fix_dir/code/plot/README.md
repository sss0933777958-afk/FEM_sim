# …/no_fix_l/code/plot/ — no-fix-ℓ 主程式的繪圖

**用途**：18-param bias 校正主程式對應的繪圖腳本。
**內容**：
- `plot_nofixl_convergence.m`（fit 收斂）。
- `plot_P1_circuit_charge.m` / `plot_P2_circuit_charge.m`（單極側視磁路 + **離軸** 等效電荷，吃 `Rum`；本批 R=150）。charge 由 `calib_bias.mat` 的 `ℓ̂·(R'·Pc_18(:,k))` 算（18-param bias，actuator→measure），有別於 fix_dir 的在軸 `ℓ·d̂`。P1 'zoom' 裁到尖端/WP 強場區。
- `plot_bias_field_err_hist.m`（18-param bias 模型 vs FEM 場向量差誤差直方圖，gap200，選項①；沿用 load_coils_actuator/select_ball/fit_bias/make_Pc/build_A。場模型 = main_Dmatrix 的 A·g_j，故誤差與其一致；存 err 到 ../../data/field_err_hist_gap200um_mueq.mat）。
- `plot_fix_vs_nofix_err_hist.m`（fix-l vs 18-param bias 疊圖比較，純載兩個 err .mat、不 addpath code/function 以避同名函式撞 path，輸出本 dir figures/）。
- `plot_frames_lattice_bias_3d.m`（measure+actuator 兩座標框 + 固定 ℓ̂ 球殼 + 上下兩層三角面 + 6 顆 **bias 離殼磁荷**；出 `frames_lattice_bias_3d.png`＋多綠 R=150µm 球版 `frames_lattice_bias_R150_3d.png`。載 mat 的 ℓ̂/ê + make_Pc 重建 `pc=R_act'·(ell·Pc)`，不需 FEM；球殼固定不隨點縮放）。
- `plot_svd_gain_iso_3d.m`（逐點 `T=S_i·Ĥ_I` SVD → gain(‖T‖_F)/iso(σ₁/σ₃) z=0 平面山丘；同 fix_dir 風格但電荷用 bias `dhat_bias=R_act'·make_Pc(ê)`；iso 翻色、z 軸自動刻度+min/max。出 `svd_gain_3d.png`/`svd_iso_3d.png`）。

**資料來源 / 流向**：讀 FEM 場（`ANSYS_data/.../coil1|coil5/standard`）+ `MATLAB_data/.../charge_fit/calibration/calib_bias.mat`（R、Pc_18、ell_hat）→ 圖存 `../../figures/`。場圖一律畫真實 FEM 節點原值，不內插。all-source：P1（下極）`B=−B_FEM`、P2（上極）keep `B=+B_FEM`，兩極尖端皆射出。

**命名 / 慣例**：新增繪圖腳本前須依 `CLAUDE.md` 的繪圖腳本規則：先確認屬哪個功能組 → 一任務一腳本、原地反覆改到使用者定案 → 定案後才存最終圖（定案前用 MCP preview，不落地）；使用者沒說「新增」就不開第二支。

**相關**：見上層 `../README.md`、`../../../../../CLAUDE.md`（繪圖腳本規則）。
