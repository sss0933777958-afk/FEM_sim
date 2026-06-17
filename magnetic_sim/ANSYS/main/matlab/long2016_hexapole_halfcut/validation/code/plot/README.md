# long2016_hexapole_halfcut/validation/code/plot — 驗證繪圖腳本

**用途**：把 validation 的 NRMSE 結果畫成疊圖。

**內容**：`plot_combo_nrmse_overlay.m`（多電流組合 NRMSE vs 取樣半徑疊在一張圖）。

**資料來源 / 流向**：讀 `../scripts/` 算出的 NRMSE 結果（/ MATLAB_data `.mat`）→ 出 PNG 到 `../../figures/`。

**命名 / 慣例**：一張定案圖一支腳本；PNG 放 `../../figures/`。**新增繪圖腳本前須依 `../../../../CLAUDE.md` 繪圖腳本規則**：先確認功能組、一任務一腳本（原地改到定案、定案前不另開）、定案後才存圖、真實節點不內插。

**相關**：見 `../../README.md`、`../../../../CLAUDE.md`。
