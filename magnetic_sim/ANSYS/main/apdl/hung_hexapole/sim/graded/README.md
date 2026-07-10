# apdl/long2016_hexapole_halfcut/sim/graded/ — 漸變 mesh 的 P1 求解（master db）

**用途**：用 region-graded mesh（極尖 + WP 核心加密）解 P1-as-SOURCE，產生供 charge-fit 用的密 mesh 場，並存出 master `.db` 讓 P2..P6 透過 `resolve/` 重用（mesh 對 6 顆極相同，只差電流）。

**內容**：
- `MT_Sim_P1_graded.txt` — graded-mesh、P1 當 source 的完整求解（建幾何 + graded mesh + coil + BC + /SOLU + 抽 WP 場），並 SAVE master db。

**資料來源 / 流向**：幾何同 baseline（Long2016 verbatim）、只改 mesh 控制；output → `ANSYS_data/long2016_hexapole_halfcut/coil1/graded/`（含 `3DMTmagneticfield.db` master + WP `.dat`，節點數與 baseline **不同**）；P2..P6 由 `../resolve/MT_Resolve_P2toP6.txt` RESUME 此 db 重解。

**命名 / 慣例**：`MT_Sim_P1_graded.txt`；source 慣例 P1（下極）電流 = -1（場朝 WP 射出）；`D,ALL,MAG,0` 必存在；mesh 控制改、尺寸不改；改動標 `[ADDED]`/`[MODIFIED]`。

**相關**：見 `../README.md`、`../resolve/README.md`、`../../postproc/dump/README.md`、`.claude/rules/{apdl-editing,charge-model-source-convention,result-read-safety}.md`。
