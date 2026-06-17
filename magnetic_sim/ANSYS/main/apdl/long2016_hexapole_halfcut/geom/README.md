# apdl/long2016_hexapole_halfcut/geom/ — 幾何建模（含 IGES 匯出 + mesh 檢視）

**用途**：建 Long Fei 半切 hexapole 的幾何（1 yoke ring + 6 protrusion + 3 下極半錐 + 3 上極完整錐），匯出 IGES，以及純 mesh 檢視（不求解）。這層是「先有幾何」的階段，不含物理/求解。

**內容**：本層只放三個子資料夾：
- `export/` — `MT_Geom_Export*.txt`：幾何 only → IGESOUT + 3 張 PNG 視圖（含 mm/metre、gap200um、sphtip、HollowProt、WithCoil 等變體）。
- `mesh/` — `MT_Mesh_{Dense,Graded,IGES}.txt`：對既有幾何/IGES 做 free-tet mesh 存 .db 供 GUI 檢視，無 coil/無 MAG BC/無 solve。
- `scripts/` — `_generate_geom_export_gap.py`：產生 export APDL 的輔助 .py。

**資料來源 / 流向**：源幾何邏輯抄自 `long2016_hexapole_full` + `hexapole-long2016` verbatim；`export/` 寫 IGES 到 `IGES/`（metre）+ `IGES_converted/`（mm，flag 6→2）；mesh `.db` 與 PNG 存到 `ANSYS_data/long2016_hexapole_halfcut/{geom_export_*,mesh/*}`。

**命名 / 慣例**：`MT_Geom_Export*` 幾何匯出、`MT_Mesh_*` mesh-only；`_mm` 後綴 = SolidWorks 相容 mm 版、無後綴 = MKS metre 版；改動標 `[ADDED]`/`[MODIFIED]`；改尺寸（R_norm、POLE_R、YOKE…）前先量對齊 CAD，不可自選 round number。

**相關**：見 `../README.md`、`../../../CLAUDE.md`、`.claude/rules/apdl-editing.md`、`doc/workflows/{cad-export,step-to-apdl,apdl-geom-build}.md`。
