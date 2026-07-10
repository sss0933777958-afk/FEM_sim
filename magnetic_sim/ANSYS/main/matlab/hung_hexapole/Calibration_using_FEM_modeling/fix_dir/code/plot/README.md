# …/fix_l/code/plot/ — fix-ℓ 主程式的繪圖

**用途**：hung fix-ℓ 校正主程式對應的繪圖腳本（gain/iso 性能量視覺化）。
**內容**：
- `plot_gain_iso_index.m`（R=500µm 球內真實 FEM 節點逐點 charge model SVD → C=∏σ_k、κ=σ₃/σ₁ **直方圖**；棕、135 bins、N/mean/CV 框、紅虛線=mean → `figures/{gain,iso}_index_R500um.png`）
- `plot_ref_planes_3d.m`（R=500µm 綠取樣球 + 藍 **ref_xy**(z=0) + 橘 **ref_xz**(y=0) 兩正交參考面，沿 x 軸相交 → `figures/ref_planes_3d.png`）
- `plot_svd_heatmaps_2d.m`（R500 球內近 ref_xy/ref_xz 面真實節點的 C/κ **離散色塊熱圖**，scatter 方形 marker 不內插、gain=jet / iso=jet（紅=高κ 中心）+ colorbar → `figures/{gain,iso}_heatmap_{xy,xz}.png`）
- `plot_gain_iso_overlay.m`（**cross-model**：hung vs long 的 R500 C/κ **疊圖**，long=深紅/hung=深藍、共用 bins + 雙均值虛線 + 簡潔 legend → `{gain,iso}_hung_vs_long_R500um.png`，**同時存 long 與 hung 兩 figures/**）

**資料來源 / 流向**：讀 `ANSYS_data/hung_hexapole/data/coilN/gap_200um/` 6-coil FEM（1 A、hung mt_constants）→ `load_coils`（air 節點、all-source flip-sink、actuator frame）→ 逐點 charge model（`data/fit_fixl_R150um_gap_200um.mat` 的 ℓ̂/gB/Khat；Pc_base=R_act·d̂）→ 圖存 `../../figures/`。場/熱圖一律畫真實 FEM 節點原值，不內插（熱圖還原 measure 框 `Pm=Pa·R_act` 對齊 ref 面）。

**命名 / 慣例**：新增繪圖腳本前須依 `CLAUDE.md` 的繪圖腳本規則：先確認屬哪個功能組 → 一任務一腳本、原地反覆改到使用者定案 → 定案後才存最終圖（定案前用 MCP preview，不落地）；使用者沒說「新增」就不開第二支。

**相關**：見上層 `../README.md`、`../../../../../CLAUDE.md`（繪圖腳本規則）。
