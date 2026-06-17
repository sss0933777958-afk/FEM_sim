# magnetic_sim/ANSYS/backup/hung/apdl/variants — 替代極形設計（alternative pole designs）（歸檔/archived）

> **本設計已歸檔（非活躍），保留供參考。** Hung hexapole（build-workflow design）。目前唯一活躍設計為 `../../../main/`。

**用途**：非當前主 pipeline 的替代極形 / 診斷腳本，保持可運作，供日後研究復用。

**內容**（代表檔）：
- `MT_Hung_Simulate_Coil[1-6].txt` + `MT_Hung_SphereModel.txt` — baseline D-shape（無 fillet），6 顆極齊。
- `MT_Hung_Simulate_Coil1_filled.txt` + `MT_Hung_SphereModel_filled.txt` + `export_poles_filled.txt` — full-round 無 fillet（僅 Coil1）。
- `MT_Hung_Simulate_Coil1_round_filleted.txt` — full-round + fillet（僅 Coil1）。
- `mesh_steel_only.txt` — 無空氣域的鋼芯 mesh QA 診斷。

**資料來源 / 流向**：與主 sim 相同（→ `../../results/`，後處理共用 `../postproc/`）。

**命名 / 慣例**：`variants` = 替代極形變體；`filled` = full-round；`SphereModel*` = 對應幾何建模。多數變體僅 Coil1 存在，要做 6-coil sweep 需自建 Coil2-6。

**相關**：見 `../README.md`（variant 狀態表）、`../sim/README.md`、上層 `../../README.md`。
