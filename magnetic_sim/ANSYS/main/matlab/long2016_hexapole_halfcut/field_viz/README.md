# long2016_hexapole_halfcut/field_viz — 磁場視覺化功能組（純繪圖）

**用途**：把 Long Fei 下極半切六極（long2016_hexapole_halfcut）的 FEM 磁場結果畫成圖 —— |B| contour、B 向量場（quiver）、磁路（magnetic circuit）streamline、磁通（flux）疊圖。**這是純繪圖功能組**，不做擬合 / 不算矩陣，只讀 FEM 場再出圖。

**內容**：
- `code/plot/` — 所有繪圖腳本（一張圖一支），例：
  - `plot_P1_topview_Bcontour.m`（P1 極座標平面 |B| contour）
  - `plot_Bvector_P1_full_circuit.m`、`plot_Bvector_sideview_P1.m`（B 向量場 / 側視）
  - `plot_P1_magnetic_circuit_streamlines.m`（磁路 streamline）
  - `plot_P1P2_circuit_with_flux.m`、`plot_P1_circuit_with_flux_steelonly.m`（磁路 + flux）
  - `plot_Bdensity_P1_P2_overlay.m`、`plot_pole_circuit_charge_3d.m`（密度疊圖 / 3D 電荷示意）
- `figures/` — 已定案輸出 PNG（例 `P1_circuit_charge_R50_zoom.png`、`Bdensity_P1_P2_overlay.png`）。

**資料來源 / 流向**：讀 ANSYS_data 的 `coilN/coilN_{coord,bfield}_*.dat`（透過 `ansys_path` / `import_ansys_data` + `filter_iron_nodes`）→ 轉 WP frame、套 source 符號慣例 → 直接畫真實 FEM 節點 → 輸出 PNG 到本組 `figures/`。FEM 在 1A 解，繪圖不縮放。

**命名 / 慣例**：功能組 schema = `code/plot/` + `figures/`（純繪圖組無 `scripts/` / `results/`）；一張已定案圖對應一支腳本；輸出 PNG 才放 `figures/`。場圖一律畫真實 FEM 節點原值，**不用 scatteredInterpolant / 格點內插**（除非使用者明確要求並在圖說標示）。source 符號慣例見 `.claude/rules/charge-model-source-convention.md`。

**相關**：見上層 `../README.md`（matlab/ 結構與 resolver）、`../../CLAUDE.md`（繪圖腳本規則全文）；讀結果防呆見 `.claude/rules/result-read-safety.md`。
