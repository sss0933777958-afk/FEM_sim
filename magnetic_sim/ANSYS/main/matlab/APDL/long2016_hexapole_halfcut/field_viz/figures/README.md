# long2016_hexapole_halfcut/field_viz/figures — 場視覺化定案輸出圖

**用途**：放 field_viz 已定案的輸出 PNG（場圖最終版）。

**內容**：例 `P1_topview_*`、`P1_circuit_charge_R50_zoom.png`、`P1_circuit_with_flux_steelonly.png`、`Bdensity_P1_P2_overlay.png`、`P2_circuit_with_flux.png`。

**資料來源 / 流向**：由 `../code/plot/*.m` 讀 FEM `.dat` 算後輸出至此。

**命名 / 慣例**：圖檔名對應其繪圖腳本；**定案後才存最終圖**（定案前用 MCP preview 討論，不落地），全是真實 FEM 節點原值（非內插）。

**相關**：見 `../README.md`、`../../../CLAUDE.md`（Figure Production 規則）。
