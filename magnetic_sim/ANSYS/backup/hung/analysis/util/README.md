# magnetic_sim/ANSYS/backup/hung/analysis/util — 小工具（helpers）（歸檔/archived）

> **本設計已歸檔（非活躍），保留供參考。** Hung hexapole（build-workflow design）。目前唯一活躍設計為 `../../../main/`。

**用途**：pipeline 周邊的小型輔助腳本（非主流程）。

**內容**：
- `print_J_positions.m` — 從 [J] ideal K_I 擬合結果印出各極 position 表。

**資料來源 / 流向**：讀 `../../data/KI_fit.mat`（或對應 [J] 結果 `.mat`）→ 輸出到 console，不落地檔案。

**相關**：見 `../README.md`（pipeline 與檔名對照）、上層 `../../README.md`。
