# apdl/long2016_hexapole_halfcut/sim/ — 求解腳本（baseline + 變體 + mesh）

**用途**：半切 hexapole 的 ANSYS 求解 input。建幾何 + air domain + coil（SOURC36）+ MAG BC + /SOLU + 抽 WP 場。本層依「求解策略」分子資料夾。

**內容**：
- `baseline/` — `MT_Sim_P1..P6.txt`：Long Fei verbatim 基準（含 VADD 合併的正確 topology），6 顆極自激解。
- `gap200um_mueq/` — μ_r 等效氣隙變體（MAT_PROT μ_r=31 模 7mm post + 200µm gap 串聯磁阻）。
- `graded/` — region-graded mesh 的 P1-as-SOURCE 求解。
- `mesh/` — 從 baseline mesh 抽鐵 / mesh 檢視（`MT_Mesh_Iron`、`MT_Mesh_P1`）。
- `resolve/` — RESUME P1 graded master `.db`、只改 coil 電流重解 P2..P6（跳過 meshing）。
- `scripts/` — 產生上述 sim .txt 的 `.py` 輔助。

**資料來源 / 流向**：幾何取自 `hexapole-long2016` Long2016 source verbatim；解出 `.dat`（場）/ `.db`（模型）存到 `ANSYS_data/long2016_hexapole_halfcut/coilN/{standard,graded,gap200um_mueq}/`。

**命名 / 慣例**：`MT_Sim_P<N>*`，6 顆極只差 `CURR_ARRAY`（excited coil=1 或依 source 慣例 ±1，其餘 0），其餘內容同步；`D,ALL,MAG,0` 邊界必存在於 `/SOLU` 前；source 慣例＝每顆激發極場朝 WP 射出（下極 P1/P3/P6=-1、上極 P2/P4/P5=+1）；改動標 `[ADDED]`/`[MODIFIED]`；改尺寸先對齊 CAD。

**相關**：見 `../README.md`、`../../../CLAUDE.md`、`.claude/rules/{apdl-editing,fit-current-matches-sim,result-read-safety}.md`、`doc/workflows/apdl-fem-run.md`。
