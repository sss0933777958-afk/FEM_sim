# magnetic_sim/ANSYS/backup/hung/apdl/geom — 幾何建模 + IGES 匯出（geometry & IGES export）（歸檔/archived）

> **本設計已歸檔（非活躍），保留供參考。** Hung hexapole（build-workflow design）。目前唯一活躍設計為 `../../../main/`。

**用途**：建立 Hung hexapole 主幾何（D-shape + 40 µm fillet）並匯出 IGES。

**內容**（代表檔）：
- `MT_Hung_Assembly_Dfillet.txt` — **主幾何腳本**，建完整 hexapole 並匯出 `Full_Assembly_filleted.iges`（cone semi-angle 11.31°，junction 15.793 mm）。`_l250` = ℓ=250 µm 變體。
- `export_pole_filleted.txt` / `export_pole_RoundFillet.txt` — 匯出單一極 IGES。
- `export_parts.txt` — 各零件分別匯出至 `../../IGES/`（尺寸寫死，零件變更後須更新並重跑）。

**資料來源 / 流向**：APDL → 匯出 IGES 至 `../../IGES/`（再同步 `../../IGES_converted/`）；幾何供 `../sim/` 求解使用。

**命名 / 慣例**：APDL 乘號 `*` 前後不得有空格；IGES 匯出 MM=1/25.4（SolidWorks 相容）、模擬 MM=1e-3（MKS）。`Dfillet` = D-shape + fillet 主設計。

**相關**：見 `../README.md`（含 legacy 路徑對照）、上層 `../../README.md`。
