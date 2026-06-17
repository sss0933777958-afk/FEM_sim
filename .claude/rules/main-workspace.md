# Kuo Quadrupole 工作目錄規則

當工作涉及 `magnetic_sim/ANSYS/main/` 目錄（cwd 在 `magnetic_sim/ANSYS/main/`、討論 Quadrupole、編輯 `magnetic_sim/ANSYS/main/*` 檔案）時，**所有新建的模擬結果、程式碼、圖檔都必須寫進 `magnetic_sim/ANSYS/main/` 下對應子目錄**，不得寫到 `magnetic_sim/ANSYS/backup/hexapole-long2016/`、`magnetic_sim/ANSYS/backup/hung/`、外部 `G:\my_workspace\report\...` 或 git root。

## 兩層結構：canonical 階段 + topic 子資料夾

`magnetic_sim/ANSYS/main/` 下每個 canonical 子目錄都用 **topic 子資料夾**進一步分群：

- `kuo_quadrupole/` — Kuo 自己的 4-pole MEMS Quadrupole 工作（**預設 topic**）
- `<other_topic>/` — 其他 topic（如 `long2016_h1h2/` 借 Long Fei 幾何做的驗證）

## 正規輸出位置

| 類型 | 目的地（Kuo Quadrupole 工作） |
|---|---|
| APDL 幾何腳本 | `magnetic_sim/ANSYS/main/apdl/kuo_quadrupole/geom/` |
| APDL 模擬腳本 | `magnetic_sim/ANSYS/main/apdl/kuo_quadrupole/sim/` |
| APDL 後處理腳本 | `magnetic_sim/ANSYS/main/apdl/kuo_quadrupole/postproc/` |
| APDL sweep / batch | `magnetic_sim/ANSYS/main/apdl/kuo_quadrupole/sweep/` |
| ANSYS 結果（.rst, .rmg, .db, .out, .dat） | `magnetic_sim/ANSYS/main/ANSYS_data/kuo_quadrupole/<case_tag>/` |
| MATLAB 分析腳本 | `magnetic_sim/ANSYS/main/analysis/kuo_quadrupole/{core,fit,plot}/` |
| Fitting 資料（.mat, .csv） | `magnetic_sim/ANSYS/main/ANSYS_data/kuo_quadrupole/` |
| COMSOL 腳本（.m, .ps1） | `magnetic_sim/ANSYS/main/comsol/kuo_quadrupole/` |
| COMSOL 模型（.mph） | `magnetic_sim/ANSYS/main/mph/kuo_quadrupole/` |
| SEMulator 流程與輸出 | `magnetic_sim/ANSYS/main/semulator/kuo_quadrupole/` |
| 幾何匯出 IGES | `magnetic_sim/ANSYS/main/IGES/kuo_quadrupole/` + `magnetic_sim/ANSYS/main/IGES_converted/kuo_quadrupole/`（必須同步） |
| SolidWorks / STEP 原檔 | `magnetic_sim/ANSYS/main/CAD/kuo_quadrupole/{SLDPRT,STEP}/` |
| 所有圖檔（含報告圖） | `magnetic_sim/ANSYS/main/figures/<topic>/<case_tag>/`（多 variant 設計，如 `kuo_quadrupole/Lp046_T55_R0500/`）或 `magnetic_sim/ANSYS/main/figures/<topic>/`（單一 configuration，如 `long2016_p1_only/`）；`<case_tag>` 不可帶方法後綴（如 `_jfit`、`_fit`、`_post` — 方法不是模型）；**同一 case 不論 fit/sim/post/plot 圖都進同一個 `<case_tag>/`** |
| LaTeX 腳本（.tex 原始檔） | `magnetic_sim/ANSYS/main/doc/<analysis>/<topic>/scripts/`（**analysis-first** schema） |
| 編譯後 PDF / 文件 | `magnetic_sim/ANSYS/main/doc/<analysis>/<topic>/pdf/`（所有後續文件相關產物都進這裡） |
| **跨 analysis 流程 SOP**（Claude 可 follow 的通用文件） | `magnetic_sim/ANSYS/main/doc/workflows/`（cross-topic, **不分 model**；目前 14 份:README + Round 1 4 份 + Round 2 9 份） |

### `<topic>` 命名規則（doc 第二層）

- **第二層 `<topic>` 名稱必須跟 `magnetic_sim/ANSYS/main/figures/` 頂層 topic 完全一致**（例:`long2016_hexapole_halfcut` / `long2016_hexapole_full` / `kuo_quadrupole` / `zhang_quadrupole` / ...）
- **`<topic>` 必須是「模型」命名**（從檔名能看出哪個物理模型）— 不可用「活動」命名（如 `long2016_charge_fit` / `long2016_h1h2`），活動應歸到 model 的 sub-folder
- 若 `magnetic_sim/ANSYS/main/figures/` 還沒有對應 topic 夾，**先去 `magnetic_sim/ANSYS/main/figures/` 建一個**，兩邊命名必須同步
- **唯一例外**:真的跨變體 / 純理論文件可用 `general/`（小寫、不加底線；figures/ 不必有對應夾）

### 活動（activity）放在 model 之下（適用 analysis / figures / data）

- `magnetic_sim/ANSYS/main/analysis/<model>/<activity>/` 跟 `magnetic_sim/ANSYS/main/figures/<model>/<activity>/` 跟 `magnetic_sim/ANSYS/main/ANSYS_data/<model>/<activity>/`
  都採用 **model-first、activity 為 sub-folder** schema
- **不可**用活動名當 topic（避免從檔名看不出模型）
- 例:H1/H2 比值分析在 Long Fei 完整 hexapole 上做 → `magnetic_sim/ANSYS/main/analysis/long2016_hexapole_full/h1h2/`，
  不要建 `magnetic_sim/ANSYS/main/analysis/long2016_h1h2/`
- 常見活動 sub-folder 名:`charge_fit` / `h1h2` / `bs_matrix` / `eddy_current` / ...
| 參考論文 / PDF | `magnetic_sim/ANSYS/main/reference/kuo_quadrupole/<paper-folder>/` |

## 命名規則

- 在 `kuo_quadrupole/` 之內，**檔名不必再帶 `MT_Kuo_Quadrupole_` 前綴**（dir 已表明 topic）。例：
  - APDL sim：`MT_Sim_Dipole.txt`、`MT_Sim_PlanarCoil.txt`、`MT_Sim_TURNS6_Final.txt`
  - APDL postproc：`MT_Post_Fig5.txt`、`MT_Diagnose_Energy_Shells.txt`、`MT_Sweep_Post.txt`
  - APDL geom：`MT_Geom.txt`、`MT_Geom_R0500.txt`
- `<case_tag>` 命名沿用 `Lp<L_P>_T<thickness>_<excitation>` pattern（例：`Lp0p46_T55_TURNS6`、`Lp0p50_T100_dipole`）

## 強制規則

1. 不可把新 kuo 工作產物寫到其他設計目錄（如 `magnetic_sim/ANSYS/backup/`）
2. 不可把圖檔輸出到 `G:\my_workspace\report\` 或其他 `magnetic_sim/ANSYS/main/` 之外的絕對路徑；報告需要圖時改 reference `magnetic_sim/ANSYS/main/figures/kuo_quadrupole/<topic>/<file>.png`
3. `magnetic_sim/ANSYS/main/IGES/kuo_quadrupole/` 與 `magnetic_sim/ANSYS/main/IGES_converted/kuo_quadrupole/` 必須同步：
   ```bash
   cp magnetic_sim/ANSYS/main/IGES/kuo_quadrupole/Part.iges magnetic_sim/ANSYS/main/IGES_converted/kuo_quadrupole/Part.iges
   sed -i "s/,1.0,6,,/,1.0,1,,/" magnetic_sim/ANSYS/main/IGES_converted/kuo_quadrupole/Part.iges
   ```
   刪除或重命名 `.iges` 時兩邊都要同步處理
4. **新 topic（不是 Kuo Quadrupole 自身工作）必須開新 topic 子資料夾**，不可丟到 `kuo_quadrupole/` 下混雜。例：借 Long Fei 2016 幾何的 H1/H2 驗證放 `magnetic_sim/ANSYS/main/<canonical>/long2016_h1h2/`

## 例外處理

需要寫到 `magnetic_sim/ANSYS/main/` 之外時，**必須先向使用者明確說明原因並取得同意**，例如：
- 修改 git root 的 `CLAUDE.md` / `.gitignore`（屬於 repo 設定）

## 不涉及 magnetic_sim/ANSYS/main/ 時忽略此規則
