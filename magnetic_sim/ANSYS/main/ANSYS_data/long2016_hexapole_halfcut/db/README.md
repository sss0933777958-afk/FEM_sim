# ANSYS_data/long2016_hexapole_halfcut/db/ — ANSYS 模型 `.db`

**兩層結構**（2026-08-17 起；`geom/` 已依 `ansys-db-cleanup.md` 政策整層移除）：

| 層 | 內容 | 保留 | 用途 |
|---|---|---|---|
| `mesh/` | 幾何 ＋ 網格，**未解** | 主 `.db` ＋ 主 log | **求解的種子** —— sim deck `RESUME` 它再加激發求解 |
| `sim/` | 幾何 ＋ 網格 ＋ **場** | 主 `.db` ＋ 主 `.rmg` ＋ 主 log | GUI 互動檢視；不重解就抽新物理量（`RESUME` → POST1 → 匯出）|

> **為什麼 `mesh/` 不可失去**：節點編號與座標在切網格那一刻定下來，之後所有 `.dat` 都源自它。
> 網格重建 ⇒ 節點數變 ⇒ 舊 `.dat` 與新解**不可比**（`RESULTS_MAP.md` 的 494873 / 390579 節點指紋、
> `gap_200um` vs `no_gap` 同節點對照，都建立在此）。`geom/` 則完全由 deck 決定，可隨時重生。

## `mesh/`（12 case，56 GB）

| case | `.db` | 說明 |
|---|---|---|
| `mesh_baseline` | `mesh_baseline.db` | **基準網格** —— 產生 494873 / 390579 節點指紋 |
| `graded` | `mesh_graded.db` | charge-fit 用密網格 |
| `mesh_graded_basegap` | 13 顆（`mesh_gap*` / `mesh_protgap` / `exp*` / `gck*`）| 各 gap 變體的網格種子 |
| `singlepole` | 6 顆，分 `filled/` `filled_gap_x42a150/` `filled_gap_x42p5/` `filled_gap_x47p5/` `halfcut/` `tipcut/` | 單極三形狀＋gap 變體 |
| `sensor_refine` / `sensor_refine_p1p2` | 各 3 顆（`mesh_lv1/2/3.db`）| sensor 局部加密三等級 |
| `lower_filled` / `tip400um` / `vp` / `vp6` / `sensor_misplace` / `coarse` | 各 1 顆 | 其餘變體 |

## `sim/`（14 case，142 GB）

`graded`（六顆 coil）、`tip400um`（六顆）、`gap_basegap_solve`、`protgap_solve`、`gap_sweep`、
`gap400_solve`、`sensor_refine`(_p1p2)、`singlepole`、`lower_filled`、`graded_rsp`、`vp`、`vp6`。

## 慣例

- sim 解完**立刻**把 `.dat` 複製到 `../data/coil<N>/<variant>/` —— `.dat` 才是 MATLAB 讀的交付物。
- 副產物（`.esav` `.full` `.dbe` `.page` `.lock` `.DSP*` `.bat`、per-worker `.rmg`/log）**一律清掉**。
- 🔴 **worker 判別只有一個正確測試**：去掉 stem 尾數後，同夾必須存在該名的主檔。
  **不可只看「以數字結尾」** —— `sim_coil1..6.rmg`、`coil1_lv1/2/3`、`mesh_lv1.out` 都是**主檔**。
- ANSYS 輸出（`.db` / `.rmg` / `.dat`）**不進 git**（見 repo `.gitignore`）。

**相關**：`../data/README.md`、`../RESULTS_MAP.md`、`.claude/rules/ansys-db-cleanup.md`。
