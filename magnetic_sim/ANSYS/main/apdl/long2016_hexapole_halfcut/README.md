# apdl/long2016_hexapole_halfcut/ — Long Fei 六極下極半切 hexapole 的 APDL 腳本根

**用途**：這個 model（Long Fei 6 極 hexapole，下極 P1/P3/P6 磨平半錐、上極 P2/P4/P5 完整錐）所有 ANSYS MAPDL 腳本的 topic 根。涵蓋幾何建模 + IGES 匯出（`geom/`）、求解（`sim/`）、後處理抽場（`postproc/`）、GUI 檢視（`gui/`，非求解）四類。是重跑 sim 的 input。

**內容**：本層只放子資料夾，無散落 .txt：
- `geom/` — 幾何建模：`export/`（IGES 匯出）、`mesh/`（mesh-only 檢視）、`scripts/`（產生 APDL 的 .py）。
- `sim/` — 求解：`baseline/`（Long Fei verbatim 基準）、`gap200um_mueq/`（μ_r 等效氣隙）、`graded/`（漸變 mesh）、`mesh/`（mesh 抽取）、`resolve/`（重解）、`scripts/`（.py 產生器）。
- `postproc/` — 後處理：`dump/`（從存好的 .db 重抽場）。
- `gui/` — **GUI 檢視用（非求解）**：build-only（不 mesh、不解）建模存 `.db` 供 MAPDL GUI 開來檢視。`MT_AllSource.txt`＝6 極 all-source 電流方向（下極 CURR=−1、上極 +1），存 `ANSYS_data/.../db/allsource/allsource.db`。

**資料來源 / 流向**：input 幾何來自 `IGES_converted/long2016_hexapole_halfcut/`；解出的 `.dat`（場）/ `.db`（模型）存到 `ANSYS_data/long2016_hexapole_halfcut/<case>/`（如 `coilN/standard`、`coilN/graded`、`coilN/gap200um_mueq`）。

**命名 / 慣例**：`MT_Geom*` = 幾何、`MT_Sim_*` = 求解、`MT_Mesh_*` = mesh、`MT_Dump_*` = 後處理；6 顆 coil 腳本只差 `CURR_ARRAY`（coil N = 1，其餘 0）；改動標 `[ADDED]`/`[MODIFIED]`（英文註解）；`D,ALL,MAG,0` 邊界必存在於 `/SOLU` 前；改幾何尺寸先對齊 CAD（`ansys-cad-alignment.md`）；極一律用紙上名 P1-P6。

**相關**：見 `../README.md`、`../../CLAUDE.md`、`.claude/rules/apdl-editing.md`、`doc/workflows/apdl-fem-run.md`。
