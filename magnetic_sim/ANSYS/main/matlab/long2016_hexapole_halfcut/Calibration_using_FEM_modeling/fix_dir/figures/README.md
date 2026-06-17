# …/fix_dir/figures/ — 圖檔

**用途**：fix-ℓ 校正主程式的繪圖輸出（直接放此層，不再多包子夾）。

**內容**（側視 xz 磁路箭頭 + 在軸等效電荷；真實 FEM 節點不內插、格點抽樣、turbo 依 |B|、真實極輪廓、洋紅小圓標 charge）：
- `P1_circuit_charge_R150_zoom.png` — P1（下極，coil1 自激）。all-source `B=−B_FEM`（下極翻號→尖端射出 source）；charge `q_P1=ℓ·d̂₁`（在軸，ℓ≈0.852 mm）。視窗裁到尖端/WP 強場區 [-2,2]×[-14,-11.3]（後端 graded-mesh 稀疏帶已裁掉）。
- `P2_circuit_charge_R150_zoom.png` — P2（上極，coil5 自激）。all-source 上極 keep `B=+B_FEM`（尖端射出）；charge `q_P2=ℓ·d̂₂`（在軸）；真實完整錐面、沿極軸傾斜。

**產生**：`../code/plot/plot_P1_circuit_charge.m('zoom',true,150)`、`plot_P2_circuit_charge.m(true,150)`；
讀 `ANSYS_data/.../coil1|coil5/standard` 場 + `MATLAB_data/.../charge_fit/fit_KI_ball/fit_KI_R150.mat`（ℓ）。

**相關**：見上層 `../README.md`、`../code/plot/README.md`。
