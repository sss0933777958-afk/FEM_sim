# magnetic_sim/ANSYS/backup/hung/scripts/run/variants — 舊版批次腳本（superseded run scripts）（歸檔/archived）

> **本設計已歸檔（非活躍），保留供參考。** Hung hexapole（build-workflow design）。目前唯一活躍設計為 `../../../../main/`。

**用途**：被取代的舊版 Coil2–6 批次執行腳本，保留供參考（非當前流程）。

**內容**：
- `run_coil2to6.sh` — v1：只求解、無 POST1 匯出。
- `run_coil2to6_refined.sh` — v2：加了 POST1，但用 `-m 4000`。

**資料來源 / 流向**：與當前版相同（呼叫 `../../../apdl/sim/` + `../../../apdl/postproc/`，輸出至 `../../../results/`）。

**命名 / 慣例**：`variants` = 被取代的舊版；當前版在上一層 `../`（`run_coil2to6_refined_v2.sh` = 主流程，`-m 8000` + 預清 `*.rmg`）。

**相關**：見 `../README.md`（版本沿革）、`../../README.md`、上層設計 README。
