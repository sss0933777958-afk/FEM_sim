# long2016_hexapole_halfcut/field_viz/code/plot — 場視覺化繪圖腳本

**用途**：放 field_viz 的繪圖腳本（一張定案圖一支）：|B| contour、B 向量場（quiver）、磁路 streamline、磁通疊圖。

**內容**：例 `plot_P1_topview_Bcontour.m`（P1 平面 |B| contour）、`plot_Bvector_P1_full_circuit.m`（P1 全磁路向量）、`plot_P1_magnetic_circuit_streamlines.m`（streamline）、`plot_P1P2_circuit_with_flux.m`（磁路 + flux）、`plot_Bdensity_P1_P2_overlay.m`（密度疊圖）。

**資料來源 / 流向**：讀 ANSYS_data `coilN/*.dat`（`import_ansys_data` + `filter_iron_nodes`）→ WP frame + source 符號 → 畫真實節點 → 輸出 PNG 到 `../../figures/`。

**命名 / 慣例**：`plot_<物理>_<極/視角>.m`；輸出 PNG 放 `../../figures/`。**新增繪圖腳本前須依 `../../../CLAUDE.md` 繪圖腳本規則**：先跟使用者確認屬哪個功能組、一任務一支腳本（原地反覆改到定案、定案前不另開）、定案後才存圖、場圖一律真實 FEM 節點原值不內插。

**相關**：見 `../../README.md`、`../../../CLAUDE.md`（繪圖規則全文）。
