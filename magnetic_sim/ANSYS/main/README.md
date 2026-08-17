# magnetic_sim/ANSYS/main/ — 活躍工作區

六極電磁微探針的 FEM 模擬 + 集總參數（點磁荷）模型校正。目前**唯一活躍**的工作根。

**完整規則與導覽見本夾的 [`CLAUDE.md`](CLAUDE.md)**（架構地圖、資料流、繪圖規則、鐵則、Quick Triggers）。
本檔只是**人看的入口**：這裡有什麼、東西放哪、怎麼跑。

## 三個 model topic

| topic | 是什麼 | 狀態 |
|---|---|---|
| `long2016_hexapole_halfcut` | Long Fei 2016 論文六極（下極半切錐、上極全錐） | **主力** |
| `hung_hexapole` | Hung 六極（R300 / R500 / R700 三種工作半徑） | 對照 |
| `NTU_hexapole` | NTU 扁平板六極 | 對照 |

幾何鐵則（所有六極通用、不可協商）：三對極軸互相正交、六極尖共球（`R_norm`）、上下層方位差 60°
⇒ `alpha = arctan(√2) = 54.74°` **鎖死**，`R_norm_xy = R_norm·√(2/3)`、`R_norm_z = R_norm/√3`。
下極 P1/P3/P6 在 0/120/240°、上極 P2/P4/P5 在 60/180/300°。

## 頂層夾（各有自己的 README，動該夾前先讀）

| 夾 | 放什麼 | 什麼時候來這裡 |
|---|---|---|
| `CAD_model/` | SolidWorks 原檔 + STEP —— **幾何 source of truth** | 量尺寸、對齊 ANSYS 幾何 |
| `model_check/` | 交付檢查用 mm STEP（+ 舊 IGES） | 出檔給人疊 CAD 檢查、`IGESIN` 匯入 |
| `apdl/<model>/` | APDL 腳本 `{geom,mesh,sim,postproc,gui}/` | 改幾何 / 參數 / 重跑 sim |
| `ANSYS_data/<model>/` | FEM 輸出：`data/<variant>/coilN/*.dat`（場）、`db/{geom,mesh,sim}/`（模型） | **讀 FEM 結果** |
| `matlab/` | 分析碼，第一層分 `APDL/` 與 `Maxwell/` 兩個**求解器分支** | 校正、算矩陣、畫圖、讀寫 `.mat` |
| `figures/` | `paper_fig/Section*/`（論文圖）+ `paper_fig_plot/{plot,data}`（產生端）+ `.pptx` | 論文圖 |
| `doc/` | LaTeX 推導 / 報告 PDF + `workflows/`（操作 SOP） | 推導、跑流程前查 SOP |
| `.claude/` | Claude Code 本地設定 | 通常不動 |

> ⚠ `IGES/` 與 `MATLAB_data/` **已移除**：交付一律出 STEP（`model_check/`），`.mat` 一律放
> `matlab/<分支>/data/<model>/.mat/`。
> ⚠ `temp_figures/` 是歷史暫存夾、不在架構地圖內，勿把新產物寫進去。

## 資料流

```
CAD_model (SLDPRT/STEP，source of truth)
  → apdl/<model>/{geom,mesh}      建幾何 + 網格
  → apdl/<model>/sim              激發 coil1..6，[MAPDL 求解]
  → ANSYS_data/<model>/data/<variant>/coilN/*.dat        場（measure frame、Tesla）
  → matlab/{APDL,Maxwell}/.../main/main.m                校正管線
  → data/<model>/.mat  +  results/<model>/*.pdf  +  figures/
  → figures/paper_fig_plot/plot/*.m  →  figures/paper_fig/Section*/*.png
  → doc/<主題>/                                          報告
```

## 跑一次校正

```matlab
% 改 main.m 頂部旗標 → 執行；跑完自動存 .mat 並出 PDF
run('matlab/Maxwell/main/main.m')                                   % Maxwell 分支
run('matlab/APDL/Calibration_using_FEM_modeling/main/main.m')       % APDL 分支
```
輸出 = `ᴮĤ_I = ĝ_I·K̄_I`（mT/A）或 `ᴮĤ_V = ĝ_V·D̄`（mT/mV），連同 ℓ̂、G、V、𝒞/κ、RMSPE。

## 跑 ANSYS

本機 MAPDL：`G:\ANSYS Inc\v252\ansys\bin\winx64\MAPDL.exe`（換機先確認路徑存在）

```bash
ANSYS="G:\ANSYS Inc\v252\ansys\bin\winx64\MAPDL.exe"
"$ANSYS" -b -np 4 -m 24000 -dir "<結果夾>" -j "coil1" -i "<deck>.txt" -o "<結果夾>/solve.out"
```

## 最重要的幾條規則

1. **讀 FEM 結果前先查 `ANSYS_data/<model>/RESULTS_MAP.md`** 並核指紋（節點數 + |B| max）——
   `standard` 與 `gap_200um` 節點數相同，只能靠 |B| 差 ~30% 區分。
2. **不擅自更動檔案架構**（移動 / 改名 / 刪除 / 新建夾）——一律先問。改既有檔內文不受此限。
3. **禁止使用 `backup/` 的資料與程式**；常數一律走 `model_config(model, geom)`。
4. **改 ANSYS 幾何前先量 CAD**；不一致不可自己選值，通報使用者拍板。
5. **畫圖前先問風格選項**；場圖畫真實 FEM 節點、不內插（除非明示並標示）。
6. **清 db / sim 前先讀 `ansys-db-cleanup.md`**（`geom/` 整層刪、`mesh/`+`sim/` 留主檔）。
7. 對外一律用 **P1–P6**（paper 名）；APDL coil index 只在改 deck / raw 脈絡時提。

完整清單見 `CLAUDE.md` 的 Quick Triggers 與 `../../.claude/rules/`（21 條）。
