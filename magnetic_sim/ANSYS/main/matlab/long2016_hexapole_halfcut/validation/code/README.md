# long2016_hexapole_halfcut/validation/code — 程式碼層

**用途**：validation 功能組的程式碼容器，分驗證運算與繪圖。

**內容**：`scripts/`（NRMSE 評估 / 逐點誤差腳本）、`plot/`（NRMSE 疊圖）。

**資料來源 / 流向**：scripts 讀 ANSYS_data `.dat` + fixl_fit 擬合參數 → 算 NRMSE；plot 讀結果 → PNG（`../figures/`）。

**命名 / 慣例**：多腳本組分 `code/scripts/`（運算）與 `code/plot/`（繪圖）。

**相關**：見 `../README.md`、`../../../CLAUDE.md`。
