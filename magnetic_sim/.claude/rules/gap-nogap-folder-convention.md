# gap/nogap 變體資料夾命名（強制、別再增殖）

**使用者拍板（2026-07-10，動怒點）**：每個 hexapole model 的「有氣隙 / 無氣隙」FEM 變體，**跨 `data/coilN/`、`db/`、`apdl/sim/`、`apdl/geom/export/` 一律只用兩個名字**：

| 意義 | 資料夾/檔名 tag |
|---|---|
| 有氣隙（canonical = **最近定案結果**） | **`gap_200um`** |
| 無氣隙（同網格 slab 翻鋼、與 gap 同節點對照） | **`no_gap`** |

對應 memory：`feedback_gap_nogap_folder_convention`
相關規則：`modify-existing-files.md`、`no-structure-change-without-ask`、`ansys-db-cleanup.md`。

---

## 🔒 核心規則

1. **只准 `gap_200um` + `no_gap` 兩個變體夾**（每 model、每層）。canonical gap = 最近拍板的幾何/結果（例：hung 目前 = 縮鐵棒38 + radial40 + 100mm² partial）；換了定案就**覆蓋** `gap_200um`，不另開名。
2. **禁止增殖平行 gap 名夾**。歷史上亂長過的（**別再犯**）：`gap200um_mueq`、`gap100um_mueq`、`lowfwd_gap`、`nogap_graded`、`nogap_samemesh`、`gap_basegap_solve` … → 全部收斂成 `gap_200um` / `no_gap`（已於 2026-07-10 清理 long+hung）。
3. **要開「新 variant 夾」之前先問使用者**。跑實驗/換位置/換方法時：
   - 首選**覆蓋現有 canonical**（`gap_200um` / `no_gap`），產物落新夾只在「要保留對照」時、且**先問**。
   - throwaway 探索用 scratchpad 或臨時 sed，**跑完清掉、不落地成第三個 gap 夾**。
4. **冗餘一律刪**：被取代的舊 gap/nogap 變體（data 夾、db solve scratch、apdl/sim 夾、geom deck 檔）清乾淨，只留 canonical + 還在用的非-gap 分析（如 `baseline`/`standard`、`singlepole`、`graded`）。刪前 dry-run 列清單、存疑先問。

## 命名一致性（跨層）
- `data/coilN/{gap_200um,no_gap}/`、`apdl/<m>/sim/{gap_200um,no_gap}/`、fit `.mat` = `fit_fixl_R*_{gap_200um,no_gap}.mat`、MATLAB `VARIANT='gap_200um'|'no_gap'`、deck `/CWD`+POST 輸出路徑 = 同夾名。**改一個就同步全部**（deck 內路徑 = 它的夾名）。

## 觸發片語（任一即套用）
- 「跑 gap / no-gap sim」「存 gap 結果」「另開變體夾」「換氣隙位置再跑」
- 你正打算建第三個含 `gap`/`nogap` 的夾時 → 停，覆蓋 canonical 或先問。

## 何時不適用
- 非 gap/nogap 的獨立分析（`singlepole`、`sensor` 等本來就別的模型/主題）——各自命名，不受此兩夾限制。
- 純 mesh `.db` 夾（`mesh_graded*`）——mesh 不綁 variant 名，另循 `ansys-db-cleanup`。
