# apdl/long2016_hexapole_halfcut/sim/resolve/ — 重用 master db 重解（跳過 meshing）

**用途**：6 顆極 mesh 完全相同、只差 coil 電流，所以 RESUME P1 graded master `.db`，只改 SOURC36 real-constant 電流逐極重解、抽 WP 場。**跳過（慢的）meshing**，省時。

**內容**：
- `MT_Resolve_P2toP6.txt` — RESUME `coil1/graded/3DMTmagneticfield.db`，依 source 慣例改電流，依序解 P2(coil5,+1)/P3(coil2,-1)/P4(coil6,+1)/P5(coil4,+1)/P6(coil3,-1)，dump 各極 WP `.dat`。

**資料來源 / 流向**：input = `ANSYS_data/long2016_hexapole_halfcut/coil1/graded/3DMTmagneticfield.db`（master，由 `../graded/MT_Sim_P1_graded.txt` 產生）；此腳本為一次 COMBINED solve，先寫到 `coil2/graded/`，跑完須把 coil3/4/5/6 的 `_wp.dat` 搬到各自 `coilN/graded/`。

**命名 / 慣例**：`MT_Resolve_<poles>.txt`；source 慣例＝下極 P1/P3/P6 電流 -1、上極 P2/P4/P5 電流 +1（場朝 WP 射出）；只改電流不改幾何/mesh；RESUME 後防禦性重定義常數；改動標 `[ADDED]`/`[MODIFIED]`。

**相關**：見 `../README.md`、`../graded/README.md`、`.claude/rules/{apdl-editing,charge-model-source-convention,result-read-safety}.md`。
