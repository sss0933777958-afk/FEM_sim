# …/fix_dir/figures/ — 圖檔

**用途**：fix-ℓ 校正主程式的繪圖輸出（直接放此層，不再多包子夾）。

**內容**（側視 xz 磁路箭頭 + 在軸等效電荷；真實 FEM 節點不內插、格點抽樣、turbo 依 |B|、真實極輪廓、洋紅小圓標 charge）：
- `P1_circuit_charge_R150_zoom.png` — P1（下極，coil1 自激）。all-source `B=−B_FEM`（下極翻號→尖端射出 source）；charge `q_P1=ℓ·d̂₁`（在軸，ℓ≈0.852 mm）。視窗裁到尖端/WP 強場區 [-2,2]×[-14,-11.3]（後端 graded-mesh 稀疏帶已裁掉）。
- `P2_circuit_charge_R150_zoom.png` — P2（上極，coil5 自激）。all-source 上極 keep `B=+B_FEM`（尖端射出）；charge `q_P2=ℓ·d̂₂`（在軸）；真實完整錐面、沿極軸傾斜。

**場誤差直方圖**（選項①粗體框圖；圖標只標 N+median，無 region err）：
- `charge_field_err_hist_gap200um_mueq.png` — fix-ℓ 模型 vs FEM（gap200）逐點逐激發向量差 |B_model−B_FEM|（6 激發合併，R=150 µm 球內 8774 點×6）。橫軸 mT、縱軸 count；median≈0.183 mT（參考 region err 3.19%，不顯示於圖）。由 `../code/plot/plot_charge_field_err_hist.m` 產生（同時存 `../data/field_err_hist_gap200um_mueq.mat` 供疊圖用）。
- `fix_vs_nofix_err_hist.png` — **single-parameter（= fix-l）vs 18-param bias 疊圖比較**（半透明 + 各自 median 虛線線段）：single-parameter median 0.183 mT、18-param bias median 0.028 mT；顯示 bias 模型誤差明顯更集中近 0。圖標 legend 標 `single-parameter` / `18-param bias`。由 `../code/plot/plot_fix_vs_nofix_err_hist.m` 產生（純載兩個 err .mat、不重算）。同檔亦存於 no_fix_dir/figures/。

**SVD 致動 gain/iso 山丘 + 幾何示意圖**（model-derived、非 raw FEM；3D 框變體 B＝`box on`）：
- `svd_gain_3d.png` / `svd_iso_3d.png` — z=0 平面極座標底面、高度=`gain`(‖T‖_F，z 標籤 Gain)/`iso`(σmax/σmin=σ₁/σ₃，**翻色階**)、顏色=高度值、jet surf；z 軸自動刻度**+min/max**。逐點 `T_i=S_i·Ĥ_I` SVD。gain 碗(中心 18.23→外緣 18.9 mT/A)、iso 碗(中心 1.061→外緣 1.235)。`../code/plot/plot_svd_gain_iso_3d.m`。
- `frames_lattice_3d.png` / `frames_lattice_R150_3d.png` — measure[x_m..]+actuator[x_a..] 兩座標框 + ℓ̂≈867µm 球殼 + 6 極 P1..P6(=ℓ̂·d̂) + 上下兩層三角面 + ℓ̂ 虛線；R150 版多包綠 R=150µm 校正球。`plot_frames_lattice_3d.m`。
- `ref_planes_3d.png` — 獨立綠 R=150µm 球 + 中心(z=0)切面 + 紅 `ref`。`plot_ref_planes_3d.m`。
- `ellipsoids_yaxis_3d.png` — 沿 y 軸 7 點 (0,y,0) 的致動橢球(半軸 σ_k 沿 U(:,k)) + U 三主軸(紅σ₁/綠σ₂/藍σ₃)；daspect 等比不變形。`plot_ellipsoids_yaxis_3d.m`。
- `charge_positions_P1P2_3d.png` — P1/P2 磁荷位置(粉點+ℓ̂ 虛線到 WP)，3D 框變體 A 手動 `draw_box_edges`。`plot_charge_positions_3d.m`。

**產生**：`../code/plot/plot_P1_circuit_charge.m('zoom',true,150)`、`plot_P2_circuit_charge.m(true,150)`；
讀 `ANSYS_data/.../coil1|coil5/standard` 場 + `MATLAB_data/.../charge_fit/fit_KI_ball/fit_KI_R150.mat`（ℓ）。

**相關**：見上層 `../README.md`、`../code/plot/README.md`。
