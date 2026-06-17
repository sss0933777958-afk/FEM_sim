# apdl/long2016_hexapole_halfcut/geom/mesh/ — mesh-only 檢視腳本（不求解）

**用途**：對半切 hexapole 幾何做 mesh 並存 .db 供 GUI EPLOT 檢視。**無 coil、無 MAG BC、無 /SOLU、無 POST1** —— 純看 mesh 品質 / 密度。

**內容**：
- `MT_Mesh_Dense.txt` — 密 mesh 檢視。
- `MT_Mesh_Graded.txt` — region-graded mesh（極 0.3mm / WP sphere 0.3mm / yoke 1.5mm / outer air 4mm + 對 WP center NREFINE 多遍），幾何取自 `MT_Sim_P1.txt` verbatim、只改 mesh 控制。
- `MT_Mesh_IGES.txt` — 匯入 BREP/NURBS IGES（mm）後 free-tet（SOLID187，SMRTSIZE 6），存 .db。

**資料來源 / 流向**：幾何來自 baseline `MT_Sim_P1.txt`（verbatim）或 `IGES_converted/long2016_hexapole_halfcut/*.iges`；輸出 `mesh_*.db` + iron EPLOT PNG 到 `ANSYS_data/long2016_hexapole_halfcut/mesh/{graded,...}`。

**命名 / 慣例**：`MT_Mesh_<scheme>.txt`；mesh-only 腳本不得含 `/SOLU`；幾何 block 與 baseline 同步（只改 mesh 控制，不改尺寸）；改動標 `[ADDED]`/`[MODIFIED]`。

**相關**：見 `../README.md`、`../../sim/mesh/README.md`、`doc/workflows/apdl-geom-build.md`、`.claude/rules/apdl-editing.md`。
