# apdl/long2016_hexapole_halfcut/sim/gap200um_mueq/ — μ_r 等效 200µm 氣隙求解

**用途**：用「材料屬性 trick」模擬 protrusion 與磁極間 200 µm 氣隙，**不動幾何**（避免 V_ID / mesh 失敗）。6 支腳本各激發一顆極，抽 WP 場，得到比 baseline 低 ~30% 的磁通。

**內容**：`MT_Sim_P1_gap200um_mueq.txt` … `MT_Sim_P6_gap200um_mueq.txt`（6 支，除 `CURR_ARRAY` 外相同）。

**資料來源 / 流向**：幾何完全沿用 Long2016 source（不改尺寸）；新增 `MAT_PROT` μ_r=31（= 7mm post + 200µm air gap 與 μ_r=280 鋼串聯的等效磁阻），mesh 後把 6 個 protrusion 區的元素由 `MAT_MT` 改 reassign 成 `MAT_PROT`（限 MAT=MAT_MT 元素，避免誤改 air）；output → `ANSYS_data/long2016_hexapole_halfcut/coilN/gap200um_mueq/`（節點數與 baseline **相同**，唯一判別靠 |B|max 低 ~30%）。

**命名 / 慣例**：`MT_Sim_P<N>_gap200um_mueq.txt`；6 支只差 `CURR_ARRAY`；此為 CAD-不對齊的**合法例外**（μ_r 等效法，非改幾何）須在註解標明；`D,ALL,MAG,0` 必存在；改動標 `[ADDED]`/`[MODIFIED]`。

**相關**：見 `../README.md`、`.claude/rules/{apdl-editing,ansys-cad-alignment,result-read-safety}.md`、memory `long-fei-mueq-gap-approach`。
