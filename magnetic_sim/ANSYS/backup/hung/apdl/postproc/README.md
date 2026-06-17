# magnetic_sim/ANSYS/backup/hung/apdl/postproc — 後處理 / 抽場（post-processing）（歸檔/archived）

> **本設計已歸檔（非活躍），保留供參考。** Hung hexapole（build-workflow design）。目前唯一活躍設計為 `../../../main/`。

**用途**：求解後抽取 / 匯出 B-field `.dat`，及產生幾何 / 模型診斷圖。

**內容**（代表檔）：
- `post_export_data.txt` + `post_export_data_coil[2-6].txt` — 匯出全模型 + WP 區的座標 / B-field `.dat`（各極）。`_l250_conv` / `_round_conv` = 收斂 / 變體版本。
- `post_extract_wp.txt` — 抽 WP（原點）BX/BY/BZ/BSUM 快速量測。
- `post_trace_tips.txt` / `post_trace_circuit.txt` — 6 極 tip / P1 磁路 8 點 BSUM 診斷。
- `post_plot_geometry.txt` / `post_plot_model.txt` — 幾何 / 模型 PNG 視圖。

**資料來源 / 流向**：apdl(`../sim/` 解出 `.db/.rmg`) → 本層 POST1 匯出 `.dat` 至 `../../results/coilN/<變體>/` → 供 `../../analysis/` 讀取。

**命名 / 慣例**：`coilN` 對應各極；`_conv` = NREFINE 收斂版；輸出 `coilN_{coord,bfield}_{all,wp}.dat`（SI 單位）。被 `../../scripts/run/` 批次腳本呼叫。

**相關**：見 `../README.md`、上層 `../../README.md`、`../../scripts/run/README.md`。
