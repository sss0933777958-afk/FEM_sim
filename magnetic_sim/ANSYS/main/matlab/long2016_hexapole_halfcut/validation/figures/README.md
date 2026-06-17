# long2016_hexapole_halfcut/validation/figures — 驗證定案輸出圖

**用途**：放 validation 已定案的輸出 PNG（NRMSE 驗證圖）。

**內容**：`basis_nrmse_vs_R.png`、`basis_nrmse_6lines.png`、`combo_nrmse_C1..C5.png`、`combo_nrmse_overlay.png`、`testset_error_vs_R.png`、`validate_combos_nrmse_R130/R150.png`。

**資料來源 / 流向**：由 `../code/{scripts,plot}/*.m` 算 NRMSE 後輸出至此。

**命名 / 慣例**：圖檔名對應評估 / 繪圖腳本（`C#` = 組合、`R###` = 半徑）；定案後才存最終圖。

**相關**：見 `../README.md`、`../../../CLAUDE.md`。
