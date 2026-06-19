# …/no_fix_dir/figures/ — 圖檔

**用途**：18-param bias 校正主程式的繪圖輸出（直接放此層，不再多包子夾）。

**內容**（側視 xz 磁路箭頭 + **離軸** 等效電荷；真實 FEM 節點不內插、格點抽樣、turbo 依 |B|、真實極輪廓、洋紅小圓標 charge）：
- `P1_circuit_charge_R150_zoom.png` — P1（下極，coil1 自激）。all-source `B=−B_FEM`（尖端射出）；charge `q_P1=ℓ̂·(R'·Pc₁₈(:,1))`（**離軸 bias**，ℓ̂≈0.857 mm，Δy≈0 仍在 xz 平面）。視窗裁到尖端/WP 強場區 [-2,2]×[-14,-11.3]。
- `P2_circuit_charge_R150_zoom.png` — P2（上極，coil5 自激）。all-source 上極 keep `B=+B_FEM`（尖端射出）；charge `q_P2=ℓ̂·(R'·Pc₁₈(:,2))`（離軸 bias）；真實完整錐面傾斜。

**與 fix_dir 差別**：charge 位置為 18-param **離軸** bias（`ℓ̂·Pc_18` 在 actuator 框→旋回 measure），非 fix_dir 的在軸 `ℓ·d̂`；P1 偏移較明顯、P2 較小。

**產生**：`../code/plot/plot_P1_circuit_charge.m('zoom',true,150)`、`plot_P2_circuit_charge.m(true,150)`；
讀 `ANSYS_data/.../coil1|coil5/standard` 場 + `MATLAB_data/.../charge_fit/calibration/calib_bias.mat`（R、Pc_18、ell_hat）。

**相關**：見上層 `../README.md`、`../code/plot/README.md`。
