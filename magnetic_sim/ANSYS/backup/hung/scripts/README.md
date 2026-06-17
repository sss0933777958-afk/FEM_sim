# magnetic_sim/ANSYS/backup/hung/scripts — 批次執行 / 工具腳本（automation & utilities）（歸檔/archived）

> **本設計已歸檔（非活躍），保留供參考。** Hung hexapole（build-workflow design）。目前唯一活躍設計為 `../../../main/`。

**用途**：包裝 MAPDL 呼叫的批次自動化，與 IGES→STEP 轉檔工具。這些不含模擬邏輯（邏輯在 `../apdl/`）。

**內容**：
- `run/` — bash 批次腳本（依序跑 Coil2–6 求解 + POST1 匯出）；`run/variants/` 為舊版（見 `run/README.md`）。
- `convert_iges_to_step.py` + `run_convert.bat` — 將 IGES 批次轉為 STEP 的工具。

**資料來源 / 流向**：`run/*.sh` 對每極呼叫 `../apdl/sim/*.txt` 與 `../apdl/postproc/*.txt`，輸出至 `../results/coilN/`（log 至 `../results/logs/`）；轉檔工具讀 `../IGES*/` 產 STEP。

**命名 / 慣例**：`run/` = 主批次；`run/variants/` = 被取代的舊版。`.sh` 為 build-automation（非 simulation），故與 `.txt` 分開放。

**相關**：見 `run/README.md`、`../apdl/README.md`、上層 `../README.md`。
