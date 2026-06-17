# magnetic_sim/ANSYS/backup/hexapole-long2016/docs — 技術文件（歸檔/archived）

> **本設計已歸檔（非活躍），保留供參考。** Long Fei 2016 dissertation 六極 hexapole（full）。現役設計請見 `../../../main/`。

**用途**：本設計的技術參考——擬合方法理論、符號對照、模型驗證、模擬參數、環境與除錯。動 pipeline 程式前先讀對應文件。

**內容**：Markdown `.md`：
- `fitting-methods.md` — [A]→[J]→[B-6x] 擬合方法（[B-6x] 為最終）。
- `notation-glossary.md` — **canonical** 符號／術語對照（dissertation alignment）。
- `charge-model-fitting.md` — 點電荷模型推導（point-charge model）。
- `coil-winding-sign-convention.md` — 極性與 coil_sign 修正。
- `model-validation.md` — APDL vs dissertation 對照。
- `simulation-parameters.md` — 幾何／材料／mesh／solver。
- `ansys-environment.md` — ANSYS 安裝、batch 模式、硬體。
- `workflow.md` — 4 階段 simulation-to-publication pipeline。
- `troubleshooting.md` — 已知錯誤與修法。

**資料來源 / 流向**：說明文件，非資料流節點；描述 `apdl → results → analysis → data → figures` 全流程的規範與理論。

**命名 / 慣例**：全部符號／術語以 Long 2016 dissertation 為準（見 notation-glossary.md）；論文極名 P1–P6 用於討論，APDL coil 索引僅用於 APDL／raw data context。

**相關**：見 ../README.md、../analysis/。
