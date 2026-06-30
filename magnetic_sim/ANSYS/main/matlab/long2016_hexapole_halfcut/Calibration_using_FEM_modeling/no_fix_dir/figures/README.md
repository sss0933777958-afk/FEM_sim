# …/no_fix_dir/figures/ — 圖檔

**用途**：18-param bias 校正主程式的繪圖輸出（直接放此層，不再多包子夾）。

**內容**（側視 xz 磁路箭頭 + **離軸** 等效電荷；真實 FEM 節點不內插、格點抽樣、turbo 依 |B|、真實極輪廓、洋紅小圓標 charge）：
- `P1_circuit_charge_R150_zoom.png` — P1（下極，coil1 自激）。all-source `B=−B_FEM`（尖端射出）；charge `q_P1=ℓ̂·(R'·Pc₁₈(:,1))`（**離軸 bias**，ℓ̂≈0.857 mm，Δy≈0 仍在 xz 平面）。視窗裁到尖端/WP 強場區 [-2,2]×[-14,-11.3]。
- `P2_circuit_charge_R150_zoom.png` — P2（上極，coil5 自激）。all-source 上極 keep `B=+B_FEM`（尖端射出）；charge `q_P2=ℓ̂·(R'·Pc₁₈(:,2))`（離軸 bias）；真實完整錐面傾斜。

**與 fix_dir 差別**：charge 位置為 18-param **離軸** bias（`ℓ̂·Pc_18` 在 actuator 框→旋回 measure），非 fix_dir 的在軸 `ℓ·d̂`；P1 偏移較明顯、P2 較小。

**場誤差直方圖**（選項①粗體框圖；圖標只標 N+median，無 region err）：
- `bias_field_err_hist_gap200um_mueq.png` — 18-param bias 模型 vs FEM（gap200）逐點逐激發向量差 |B_model−B_FEM|（6 激發合併，R=150 µm 球內 8774 點×6）。橫軸 mT、縱軸 count；median≈0.028 mT（參考 region err 0.46%，不顯示於圖；場模型同 Hall_sensor 的 main_Dmatrix，故數值一致）。由 `../code/plot/plot_bias_field_err_hist.m` 產生（同時存 `../data/field_err_hist_gap200um_mueq.mat` 供疊圖用）。
- `fix_vs_nofix_err_hist.png` — **single-parameter（= fix-l）vs 18-param bias 疊圖比較**（半透明 + 各自 median 虛線線段）：single-parameter median 0.183 mT、18-param bias median 0.028 mT；顯示 bias 模型誤差明顯更集中近 0。圖標 legend 標 `single-parameter` / `18-param bias`。由 `../code/plot/plot_fix_vs_nofix_err_hist.m` 產生（純載兩個 err .mat、不重算）。同檔亦存於 fix_dir/figures/。

**產生**：`../code/plot/plot_P1_circuit_charge.m('zoom',true,150)`、`plot_P2_circuit_charge.m(true,150)`；
讀 `ANSYS_data/.../coil1|coil5/standard` 場 + `MATLAB_data/.../charge_fit/calibration/calib_bias.mat`（R、Pc_18、ell_hat）。

**相關**：見上層 `../README.md`、`../code/plot/README.md`。
