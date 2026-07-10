# …/fix_l/code/plot/ — fix-ℓ 主程式的繪圖

**用途**：hung fix-ℓ 校正主程式對應的繪圖腳本（gain/iso 性能量視覺化）。
**內容**：
- `plot_gain_iso_index.m`（R=500µm 球內真實 FEM 節點逐點 charge model SVD → C=∏σ_k、κ=σ₃/σ₁ **直方圖**；棕、135 bins、N/mean/CV 框、紅虛線=mean → `figures/{gain,iso}_index_R500um.png`）
- `plot_ref_planes_3d.m`（R=500µm 綠取樣球 + 藍 **ref_xy**(z=0) + 橘 **ref_xz**(y=0) 兩正交參考面，沿 x 軸相交 → `figures/ref_planes_3d.png`）
- `plot_svd_heatmaps_2d.m`（沿 **3 個 actuator 座標切面** `ref_{x_ay_a}/ref_{y_az_a}/ref_{x_az_a}`（基底 ea_x=P1/ea_y=P3/ea_z=P5）連續評估 charge model 的 C=∏σ / κ=σ₃/σ₁ **極座標熱圖**，r×θ 網格非 FEM 節點、jet 連續填色 + 半徑環 + colorbar，**rim 疊該面內磁極投影標記 P k**（落 0/90/180/270°，佐證 lobe 對齊磁極非 bug） → `figures/{gain,iso}_polar_{xaya,yaza,xaza}.png`）
- `plot_gain_iso_overlay.m`（**cross-model**：hung vs long 的 R500 C/κ **疊圖**，long=深紅/hung=深藍、共用 bins + 雙均值虛線 + 簡潔 legend → `{gain,iso}_hung_vs_long_R500um.png`，**同時存 long 與 hung 兩 figures/**）
- `plot_actuator_ellipsoids_3d.m`（**actuator 座標系（真實傾斜 d̂ 三軸）+ WP 中心 + 近 P1 tip 致動橢球**：measure frame 繪圖（box=x_m/y_m/z_m），三軸 x_a=d̂_P1/y_a=d̂_P3/z_a=d̂_P5（傾斜、指向 P1/P3/P5，同 frames_lattice）；橢球 = `T=(Dk/|Dk|³)·Ĥ_I` 的 SVD，WP(p=0，近球形/等向) 與 p=350µm·d̂_P1(沿 x_a 傾斜拉長)；半透明+淡網格+`camlight` 立體感；淡 R_norm=500µm 工作球殼（P1 tip 落在球面）；刻度標到 ±500µm；`view(120,25)`、box on+daspect；讀 `../../data/fit_fixl_R150um_gap_200um.mat` → `figures/actuator_frame_ellipsoids_3d.png`）

**資料來源 / 流向**：讀 `ANSYS_data/hung_hexapole/data/coilN/gap_200um/` 6-coil FEM（1 A、hung mt_constants）→ `load_coils`（air 節點、all-source flip-sink、actuator frame）→ 逐點 charge model（`data/fit_fixl_R150um_gap_200um.mat` 的 ℓ̂/gB/Khat；Pc_base=R_act·d̂）→ 圖存 `../../figures/`。直方圖用真實 FEM 節點原值不內插；極座標熱圖是 charge model 平滑擬合函數在 r×θ 網格上連續評估（非 FEM 節點），C/κ 為 rotation-invariant 純量場、沿 actuator 面取樣 `p=r(u·cosθ+v·sinθ)`、Pc 維持 measure 框不旋轉。

**命名 / 慣例**：新增繪圖腳本前須依 `CLAUDE.md` 的繪圖腳本規則：先確認屬哪個功能組 → 一任務一腳本、原地反覆改到使用者定案 → 定案後才存最終圖（定案前用 MCP preview，不落地）；使用者沒說「新增」就不開第二支。

**相關**：見上層 `../README.md`、`../../../../../CLAUDE.md`（繪圖腳本規則）。
