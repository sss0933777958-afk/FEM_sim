# magnetic_sim/ANSYS/backup/hung/analysis/core — 共用工具（shared utilities）（歸檔/archived）

> **本設計已歸檔（非活躍），保留供參考。** Hung hexapole（build-workflow design）。目前唯一活躍設計為 `../../../main/`。

**用途**：fitting / plot pipeline 共用的常數與資料載入工具，被 `../fit/`、`../plot/`、`../util/`、`../variants/` 以 `addpath('../core')` 引用。

**內容**：
- `mt_constants.m` — Hung hexapole 設計常數（R_norm、tip 座標、yoke 尺寸、pole 傾斜角）。
- `import_ansys_data.m` — 讀 ANSYS 匯出的座標 / B-field `.dat`（處理 MAPDL banner header 與黏連負號），回傳含 node_id, x,y,z, bx,by,bz, bsum 的 struct。

**資料來源 / 流向**：上游讀 `../../results/coilN/<變體>/*.dat`（ANSYS 場輸出）；本層只提供函式，不產出檔案。

**相關**：見 `../README.md`（pipeline 順序、path 機制、檔名對照）與上層 `../../README.md`。
