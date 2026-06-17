# magnetic_sim/ANSYS/backup/hexapole-long2016/figures/long — 發表圖檔（.png）子集（歸檔/archived）

> **本設計已歸檔（非活躍），保留供參考。** Long Fei 2016 dissertation 六極 hexapole（full）。現役設計請見 `../../../main/`。

**用途**：實際 `.png` 圖檔存放處——dissertation 圖重製、P2 cone 場方向分析、charge-fit 驗證、簡報 slide。

**內容**：`.png`：
- dissertation 圖重製：`fig2_3a/b/c.png`、`fig2_4a/b.png`、`fig2_6a/b.png`、`fig26_coil1_J_final.png`、`fig26_B6x_clean.png`、`fig26_B6x_validation.png`、`table_J_results.png`。
- P2 cone 場分析：`P2_Bmag_vs_s_sph10.png`、`P2_angle_vs_s_sph10.png`、`P2_cone_*_vs_d.png`、`P2_cone_iron_air_transition_vs_s.png`、`fig_P2_topview_Bvectors.png`、`cone_*` 方向定義圖。
- 驗證 / 簡報：`verify_superposition.png`、`ppt_slide1..3_*.png`。

**資料來源 / 流向**：由 `../../analysis/` 繪圖腳本（`generate_figures_2_*.m`、`plot_fig26_B6x.m`、`analyze_P2_cone_Bfield.m` 等）讀 `../../results/*.dat` + `../../data/*.mat` 產出。

**命名 / 慣例**：檔名對應 dissertation 圖號（fig2_3…）或分析主題（P2_cone…、cone_…）；`sph10` 後綴 = Coil5 球半徑 10 變體；場圖為真實 FEM 節點原值。

**相關**：見 ../README.md（figures 根）、../../README.md、../../analysis/。
