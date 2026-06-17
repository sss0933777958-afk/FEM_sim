# long2016_hexapole_halfcut — 資料夾結構說明

Long Fei 2016「下極半切」六極 hexapole 的 **FEM 模擬輸出**（`.dat` 場 / `.db` 模型）。
同一幾何、每顆磁極 1 A 激發；paper 極名 P1–P6，APDL coil 索引對應
`{P1,P3,P6,P5,P2,P4} = coil{1,2,3,4,5,6}`。
（MATLAB 分析結果不在此——在頂層 `magnetic_sim/ANSYS/main/MATLAB_data/long2016_hexapole_halfcut/`。）

> **交付檔白名單原則**：本資料夾只保留「交付檔」——`.dat`（場結果）、`.db`/`.cdb`（mesh/模型）、`.md`。
> 所有 ANSYS 過度檔（`.rmg/.full/.esav/.rst/.out/.err/.log/.bat/.DSP*/.stat/.page*/scratch` …）一律不留。
> 詳見 `.claude/rules/sim-cleanup.md` §「歸檔資料夾保留原則」。

---

## `coil1/` … `coil6/` — 各磁極的求解場（FEM `.dat`）
每顆 coil 一個資料夾，底下用「**一個網格/sim 設定 = 一個子夾**」分類（三個子夾），
`import_ansys_data` 讀的就是這些 `.dat`：

### `coilN/standard/` — baseline 標準網格解（**主資料**）
| 檔 | 內容 |
|---|---|
| `coilN_bfield_all.dat` / `coilN_coord_all.dat` | **全域** B 場 + 節點座標（all dataset，~494,873 節點）|
| `coilN_bfield_wp.dat` / `coilN_coord_wp.dat` | **工作點區** B 場 + 座標（wp dataset，~390,581 節點）|

> **標準網格（Long2016 verbatim）** 的解，是主校正（R*=150µm、ℓ̂=0.856）、B̄ 矩陣、一般分析的
> 資料來源。腳本讀法：`import_ansys_data(fullfile(<root>, 'coilN', 'standard'), 'wp'/'all', 'coilN')`。

### `coilN/gap200um_mueq/` — μ_r 等效氣隙變體
μ_r 等效 200 µm 氣隙（不切幾何，protrusion 改 μ_r=31）。同標準網格、同節點數，
但 |B| 比 baseline **低約 30%**（B̄ matrix v4 的 gap 對照）。

### `coilN/graded/` — WP 核心加密的 charge-fit 場
**WP 核心加密網格**（~182,831 wp 節點；核心 R<40 µm 內加密、核心外放粗）。
專供 **R≤40 µm 小半徑 K_I 電荷擬合**（baseline 在那麼小的球內點不夠）。
`coil1/graded/` 另含 `_circuit` slab（磁路圖用）。

---

## `mesh/` — 純網格資料庫（只 `.db`/`.cdb`，無場）
| 子夾 | 內容 | .db 檔 | 生成腳本 |
|---|---|---|---|
| `standard/` | **標準求解網格**（494,889 節點，baseline 6-coil sim 用的就是這個）| `mesh_baseline.db` | `MT_Mesh_P1.txt` |
| `standard_iron/` | 上面的**鋼鐵-only 子集**（極+body，去空氣）| `mesh_baseline_iron.db/.cdb` | `MT_Mesh_Iron.txt` |
| `dense_wp/` | **WP 核心加密**網格（R<40 µm ≥130 節點）| `mesh_dense.db` | `MT_Mesh_Dense.txt` |
| `graded/` | **分區 graded** 網格（zoned ESIZE + WP NREFINE）| `mesh_graded.db` + `mesh_graded_iron.db` | `MT_Mesh_Graded.txt` |
| `from_iges/` | 從 **IGES 匯入幾何**建的網格（對照 primitive 建模）| `halfcut_iges_mesh.db` | `MT_Mesh_IGES.txt` |

---

## `RESULTS_MAP.md` — 讀結果防呆指紋表
載入哪個 dir / dataset 前先查此表（matched 節點數 + |B| max + case_tag 指紋）。
見 `.claude/rules/result-read-safety.md`。

---

## 路徑慣例（兩個資料根分開）
- **FEM 資料**（`.dat`/`.db`）：`magnetic_sim/ANSYS/main/ANSYS_data/<model>/...` → helper `kuo_paths('<model>','coil1','standard')`。
- **MATLAB 結果**（`.mat`）：`magnetic_sim/ANSYS/main/MATLAB_data/<model>/<功能>/` → helper `kuo_matlab_paths('<model>','charge_fit')`。
- 兩個 helper 都在 `magnetic_sim/ANSYS/main/analysis/common/`，相對解析、資料夾名各集中一處。
- 舊腳本仍有硬寫絕對路徑（已批次更新到新位置），逐步遷移到 helper 中。
