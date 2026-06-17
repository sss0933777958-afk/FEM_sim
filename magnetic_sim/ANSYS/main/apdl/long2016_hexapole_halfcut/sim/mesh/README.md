# apdl/long2016_hexapole_halfcut/sim/mesh/ — mesh 抽取 / 子集（鐵 only）

**用途**：從求解用的完整 baseline mesh 抽出子集（例如只留鋼 = 6 極 + body）供檢視 / charge-fit，或單極 mesh。不求解，只處理 mesh 並存 .db/.cdb。

**內容**：
- `MT_Mesh_Iron.txt` — RESUME `mesh_baseline.db`，對非鋼（air）volume 用 `VCLEAR` 清 mesh、`EDELE` 刪 SOURC36 coil 元素，留下 494889 節點 baseline 的精確鋼子集；存 `mesh_baseline_iron.db` + `.cdb`。
- `MT_Mesh_P1.txt` — P1 相關單極 mesh。

**資料來源 / 流向**：input = `ANSYS_data/long2016_hexapole_halfcut/mesh/standard/mesh_baseline.db`；output → `ANSYS_data/long2016_hexapole_halfcut/mesh/standard_iron/`（`.db` + portable `.cdb`）。

**命名 / 慣例**：`MT_Mesh_<subset>.txt`；鋼 = volume 1-6（極）+ 8（body）= `MAT_MT`；air mesh 用 `VCLEAR`（非 `EDELE`，因 air 是 volume-attached）、coil 用 `EDELE`；改動標 `[ADDED]`/`[MODIFIED]`。

**相關**：見 `../README.md`、`../../geom/mesh/README.md`、`.claude/rules/apdl-editing.md`。
