# apdl/ — APDL 腳本（FEM 的 input）

ANSYS MAPDL 的幾何 / 模擬 / 後處理腳本。這是**重跑 sim 的 input**，改參數/幾何/mesh 都在這裡改。

## 結構：`<model>/{geom,sim,postproc}/`
```
apdl/
├── long2016_hexapole_halfcut/  {geom, sim, postproc}
├── kuo_quadrupole/             {geom, sim, postproc}
└── zhang_quadrupole/           {geom, sim, postproc}
```
- `geom/` — `MT_Geom*.txt` 幾何建模（含 IGES 匯出腳本）。
- `sim/` — `MT_Sim_*.txt` / 各 coil 求解腳本（含 baseline / 變體 / mesh）。
- `postproc/` — 抽場 / dump / PATH 等後處理腳本。

## 檔案類型
- `*.txt` — APDL 腳本（純文字 input）。
- `*.py` — 產生/轉換 APDL 的輔助腳本。

## 規則
- 改動標 `[ADDED]` / `[MODIFIED]`（English comments），見 `.claude/rules/apdl-editing.md`。
- 多 coil 腳本**只差 `CURR_ARRAY`**，其餘同步；`D,ALL,MAG,0` 邊界必存在。
- 改幾何尺寸前對齊 CAD（`ansys-cad-alignment.md`）。
- 跑 FEM / 後處理 SOP：`doc/workflows/{apdl-fem-run,apdl-postproc}.md`。
- 解出來的結果存到 `../ANSYS_data/<model>/<case>/`。
