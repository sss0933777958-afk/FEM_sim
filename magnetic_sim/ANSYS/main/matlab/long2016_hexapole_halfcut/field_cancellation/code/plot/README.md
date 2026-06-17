# …/field_cancellation/code/plot/ — 繪圖

**用途**：畫某個 source/sink 組合在 WP 中間 xy 切面的磁路圖。

**內容**：`plot_field_xy(sign_vec, tag, X,Y, BX,BY,BZ, scoreval)` —
- 給 sign 組合，疊加 `B=Σ s_j·B_FEM,j`，畫 xy 切面 quiver：箭頭=面內 (Bx,By) 單位方向、色=|B|(3D, turbo)、格點抽樣（一格一節點）；標 WP 中心 `+`、disc 邊界、6 極方位 ±（實心=+、空心=−）。
- `scoreval` = 該組合在 **WP 中心 R=50µm 球**的 **mean|B|**（標題顯示 `mean|B|_R50`）；繪圖節點是 500µm 切面（脈絡），與評分球刻意不同。場由 `../main` 傳入（真實節點不內插）。

**資料來源 / 流向**：場由 main 傳入 → 出圖 `../../figures/field_xy_<tag>.png`。

**命名 / 慣例**：依 `CLAUDE.md` 繪圖規則——真實節點不內插、一圖一腳本、定案才存最終圖。

**相關**：見上層 `../README.md`、`../../figures/README.md`。
