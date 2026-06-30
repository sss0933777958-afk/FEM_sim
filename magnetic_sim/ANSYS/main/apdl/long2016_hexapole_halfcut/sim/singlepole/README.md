# apdl/long2016_hexapole_halfcut/sim/singlepole/ — 單極模型 solve deck

**用途**：把抽出的**單極模型**（下極削平填回完整圓錐 + 4 塊支撐座 + 1 根 protrusion 鐵柱；無 yoke、無上極）做 FEM 求解。mesh/solve **分離**：mesh 由 `../../mesh/MT_Mesh_SinglePole.txt` 產 `.db`，本 deck `RESUME` 後求解。

**內容**：
- `MT_Sim_SinglePole.txt` — `RESUME db/singlepole/mesh_singlepole.db` → 設 coil **1A**（SOURC36、70 匝×1A）+ `D,ALL,MAG,0`（外圓柱側面+上下蓋）+ `magsolv,3`（DSP）→ 匯 B 場 4 檔到 `data/singlepole/`（`singlepole_{coord,bfield}_{all,tip}.dat`，tip=內球 R≤2mm @ 極尖）。

**對應 mesh deck**：`../../mesh/MT_Mesh_SinglePole.txt`（**均勻鐵件 0.3mm** + 空氣漸變 內球0.4/外圈4mm；材料 by location+size）。

**結果**（2026-06-25，1A）：mesh 704,747 節點 / 4,244,865 元素（鐵件 1.17M 均勻 0.3mm、內球 327k、外圈 2.74M）；solve ~8min；|B|max ≈ 147 mT @ 極尖，物理合理。輸出 `ANSYS_data/long2016_hexapole_halfcut/data/singlepole/`、解後 `.db` 在 `db/singlepole/sim_singlepole.db`。
> ⚠ 外圈空氣 ESIZE 4mm 被曲面 conform 過度細化成 2.74M（同 hexapole 現象）；鐵件均勻達成、解可跑，但總元素量偏大、解較慢。

**命名 / 慣例**：照 baseline sim 結構（air 域/材料/SOURC36/BC/輸出格式抄 `../baseline/MT_Sim_P1.txt`）；1A＝模型激發電流（per fit-current-matches-sim）；改動標 `[ADDED]`/`[MODIFIED]`。

**相關**：`../README.md`、`../../mesh/README.md`、`.claude/rules/{apdl-editing,result-read-safety,fit-current-matches-sim}.md`。
