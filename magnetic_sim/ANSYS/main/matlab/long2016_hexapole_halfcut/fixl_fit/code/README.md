# long2016_hexapole_halfcut/fixl_fit/code — 程式碼層

**用途**：fixl_fit 功能組的程式碼容器，分擬合腳本與繪圖。

**內容**：`scripts/`（K_I 擬合 / 掃半徑 / 校正腳本）、`plot/`（收斂 + 參數 vs R 繪圖）。

**資料來源 / 流向**：scripts 讀 ANSYS_data `.dat` → 擬合 → `.mat`（MATLAB_data）+ `.tex`（`../results/`）；plot 讀擬合結果 → PNG（`../figures/`）。

**命名 / 慣例**：多腳本組分 `code/scripts/`（運算）與 `code/plot/`（繪圖）。

**相關**：見 `../README.md`、`../../../CLAUDE.md`。
