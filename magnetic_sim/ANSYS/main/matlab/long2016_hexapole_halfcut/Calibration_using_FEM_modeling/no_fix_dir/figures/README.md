# …/no_fix_dir/figures/ — 圖檔

**用途**：18-param bias 校正主程式的繪圖輸出（直接放此層，不再多包子夾）。

**內容**（側視 xz 磁路箭頭 + **離軸** 等效電荷；真實 FEM 節點不內插、格點抽樣、turbo 依 |B|、真實極輪廓、洋紅小圓標 charge）：
- `P1_circuit_charge_R150_zoom.png` — P1（下極，coil1 自激）。all-source `B=−B_FEM`（尖端射出）；charge `q_P1=ℓ̂·(R'·Pc₁₈(:,1))`（**離軸 bias**，ℓ̂≈0.857 mm，Δy≈0 仍在 xz 平面）。視窗裁到尖端/WP 強場區 [-2,2]×[-14,-11.3]。
- `P2_circuit_charge_R150_zoom.png` — P2（上極，coil5 自激）。all-source 上極 keep `B=+B_FEM`（尖端射出）；charge `q_P2=ℓ̂·(R'·Pc₁₈(:,2))`（離軸 bias）；真實完整錐面傾斜。

**與 fix_dir 差別**：charge 位置為 18-param **離軸** bias（`ℓ̂·Pc_18` 在 actuator 框→旋回 measure），非 fix_dir 的在軸 `ℓ·d̂`；P1 偏移較明顯、P2 較小。

**場誤差直方圖**（選項①粗體框圖；圖標只標 N+median，無 region err）：
- `bias_field_err_hist_gap_200um.png` — 18-param bias 模型 vs FEM（gap200）逐點逐激發向量差 |B_model−B_FEM|（6 激發合併，R=150 µm 球內 8774 點×6）。橫軸 mT、縱軸 count；median≈0.028 mT（參考 region err 0.46%，不顯示於圖；場模型同 Hall_sensor 的 main_Dmatrix，故數值一致）。由 `../code/plot/plot_bias_field_err_hist.m` 產生（同時存 `../data/field_err_hist_gap_200um.mat` 供疊圖用）。
- `fix_vs_nofix_err_hist.png` — **single-parameter（= fix-l）vs 18-param bias 疊圖比較**（半透明 + 各自 median 虛線線段）：single-parameter median 0.183 mT、18-param bias median 0.028 mT；顯示 bias 模型誤差明顯更集中近 0。圖標 legend 標 `single-parameter` / `18-param bias`。由 `../code/plot/plot_fix_vs_nofix_err_hist.m` 產生（純載兩個 err .mat、不重算）。同檔亦存於 fix_dir/figures/。

**單極 bias 側視磁路箭頭圖**（風格①粗體框；真實 FEM 節點 grid-sample 不內插、turbo 依 |B|、極錐輪廓、WP `+`、bias 電荷粉點）：
- `singlepole_sideview_{filled,halfcut,tipcut}.png` — y=0 (xz) 面 |B| quiver + 該形狀極錐輪廓 + **bias 電荷粉點 @ (ℓ̂, ℓ̂·e_z)**（讀 `../data/singlepole_bias_fit.mat` 的 ell_um/e_z；halfcut 粉點偏下軸 ℓ̂·e_z≈−69µm）。**3 張共用同一 colorbar/clim**（共用 CAP=max prctile-92≈38mT，跨形狀可比、最大值數字不誤導；見 figure-style 同類比較共用色階規則）。**tipcut 有左上角尖端 zoom inset（純畫圖：dashed=完整尖端 / solid=切尖 + 場箭頭，無刻度/軸標/title，只留外框）**；只 inset 不另出獨立 zoom 圖。由 `../code/plot/plot_singlepole_sideview.m` 產生。

**單極 bias 場誤差直方圖**（選項①粗體框圖）：
- `singlepole_magerr_hist_R150.png` — 單極 **有-bias** 模型 vs FEM 逐節點向量差 `|B_model−B_FEM|/|B_FEM|·100%`，filled/halfcut/tipcut 三形狀疊圖（**直方圖規則 nb=180、FaceAlpha0.55、白邊、mean 虛線、headroom**）。電荷 `pc=ℓ·[1,e_y,e_z]`（讀 `../data/singlepole_bias_fit.mat`）。對照 fix_dir 無-bias 版（`fix_dir/figures/singlepole_magerr_hist_R150.png`）halfcut 孤立在 6.24%，此 bias 版三形狀收斂到 ~1%（halfcut 1.03%、filled/tipcut 0.96%）。由 `../code/plot/plot_singlepole_magerr_hist.m` 產生。

**SVD 致動指標逐節點直方圖**（原始資料真實計數、統一 135 bins、紅虛線在各自 mean、黑框 N/mean/CV%；bias 版）：
- `gain_index.png` / `iso_index.png` — R≤150µm 球內逐節點 C=∏σ_k（增益,(mT/A)³）/ κ=σ₃/σ₁（等向）分布（**藍**；8774 節點；C mean 1228/CV 5.19%、κ mean 0.876/CV 4.95%）。
- `gain_index_R500um.png` / `iso_index_R500um.png` — 同上但 **R≤500µm**（**棕** [0.55 0.35 0.17]；180423 節點；C mean 1361/CV ~10%、κ mean 0.463/CV 38.5%）。R500 已到點電荷模型邊緣（>~380µm 誤差增大），僅呈現該範圍分布。標準圖 x 軸各自尺度、乾淨整數刻度。
- `gain_index_overlay.png` / `iso_index_overlay.png` — **headline 疊圖**（= 掃描 R_k=150 frame）：R150（**淺藍** [0.45 0.68 0.90]）+ R500（**淺紅** [0.95 0.55 0.55]）；x 軸用 R500 尺度、共用 edges、半透明+白邊；legend 標半徑；「僅此」= 無統計框、無均值線。直觀對比：近中心 R150 集中且較等向（κ~0.88），R500 擴散、等向差（κ~0.46）。
- `gain_index_overlay_R{150,200,…,500}um.png` / `iso_index_overlay_R{…}um.png` — **掃描疊圖**（每 metric 8 張）：**R500 固定（淺紅）** + **R_k（淺藍）掃 150→500 step 50**；**x 軸固定用 R500 尺度**（8 張同軸可比）→ 逐張看 R_k 分布隨範圍變大、從高-κ 尖峰擴散到 = R500。legend 動態顯示該 R_k。
- 疊圖用**淺藍/淺紅**專屬配色；**標準圖維持 R150 深藍 navy [0.05 0.12 0.40]、R500 棕 [0.55 0.35 0.17]**（未被疊圖換色影響）。
- 由 `../code/plot/plot_gain_iso_index.m` 產生（效率：只 `select_ball(500e-6)` 算一次逐點 C/κ，各 R_k 用 `r<=R_k` 過濾巢狀子集；沿用 R150 fit 參數；bias 用 make_Pc(ê) 離軸電荷）。

**SVD 致動山丘 + 幾何示意（model-derived、非 raw FEM；bias 版）**：
- `svd_gain_3d.png` / `svd_iso_3d.png` — z=0 平面山丘，逐點 `T=S_i·Ĥ_I` SVD；gain(‖T‖_F,z 標籤 Gain)/iso(σ₁/σ₃,翻色)、z 軸自動刻度+min/max。電荷用 bias（`dhat_bias=R_act'·make_Pc(ê)`）。中心 gain≈18.21 mT/A、iso≈1.059。由 `../code/plot/plot_svd_gain_iso_3d.m`。
- `frames_lattice_bias_3d.png` / `frames_lattice_bias_R150_3d.png` — measure+actuator 兩框 + 固定 ℓ̂≈857µm 球殼 + 上下三角面 + **6 顆 bias 離殼磁荷**（下極 P1/P3/P6 偏移~140µm≈9°、上極 P2/P4/P5~19µm；R150 版多綠 R=150µm 校正球）。由 `../code/plot/plot_frames_lattice_bias_3d.m`。
- 結果 PDF `../results/model_results_gap_200um.pdf` 含 K̄_I/ℓ̂/G/F/ĝ_I/ê + **σ_tot=18.612 mT/A、iso_tot=1.144**（R≤150µm 8774 節點 mean）；K̄/G 已套 all-source 顯示翻號（對角全正）。

**產生**：`../code/plot/plot_P1_circuit_charge.m('zoom',true,150)`、`plot_P2_circuit_charge.m(true,150)`；
讀 `ANSYS_data/.../coil1|coil5/standard` 場 + `MATLAB_data/.../charge_fit/calibration/calib_bias.mat`（R、Pc_18、ell_hat）。

**相關**：見上層 `../README.md`、`../code/plot/README.md`。
