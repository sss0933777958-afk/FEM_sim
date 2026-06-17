# Sim 副產物清理規則(強制讀取)

當「**清理 ANSYS sim 結果 / 副產物 / intermediate / disk space / G: 滿了 / 整理 results**」相關工作觸發時,**動手前必須先讀完此規則全文**。違反即視為違規。

對應 memory:`feedback_ansys_sim_cleanup_sop.md`
Helper script:`magnetic_sim/ANSYS/main/apdl/common/clean_sim_dir.sh`

---

## 🔒 不可影響的 6 項工作(criteria for what's safe to delete)

任何清理動作**必須保證**以下 6 項工作完全不受影響。砍前 mental check 走過一遍:

| 工作 | 為什麼不受影響 |
|---|---|
| 1. MATLAB 讀 sim 結果畫圖 / fit / 算矩陣 | 全部讀 `.dat`,跟 intermediates 無關 |
| 2. 改 `.txt` 改參數 / 改幾何 / 改 mesh 重跑 | `.txt` 才是 sim input,intermediates 不是 |
| 3. 重新跑已存在的 plot / fit / analysis script | 它們讀 `.dat` + `.mat`,都還在 |
| 4. 已產出的 `.mat / .png / .tex / .pdf` | 都是獨立檔案,**從未動過** |
| 5. APDL source / IGES / postproc 腳本 | 跟結果無關,**從未動過** |
| 6. 新 sim(新 case_tag)| 新解新存,跟舊 intermediates 無關 |

→ 結論:**`.dat` + `solve.out` + `magsolv.out` + 所有 `analysis/` / `data/` / `figures/` / `doc/` / `apdl/` / `IGES/` / `CAD/` 必須一個不動**。

## ⚠ 也不可影響的 2 項能力(half-clean 必保)

使用者 2026-05-30 拍板:這兩項即使省 5 GB 也不接受失去:

| 能力 | 需要的檔 |
|---|---|
| A. 重抽新物理量(H 場、energy、新 PATH)不想 re-solve | `3DMTmagneticfield.rmg`(主,no digit)+ `.db` |
| B. ANSYS GUI 開來互動檢視 mesh / contour / quality | `3DMTmagneticfield.db` |

→ 預設 clean 模式**必須保留 `.db` + 主 `.rmg`**(half-clean)。

---

## 清理模式對照

| 模式 | 保留 | 刪除 | per sim | 失去能力 |
|---|---|---|---|---|
| **half-clean(預設)** | `.dat` + `solve.out` + `magsolv.out` + `.db` + 主 `.rmg` | per-worker `.esav/.rmg`、主 `.esav`、`.full`、log、lock | 16-20 GB → **~5 GB** | 只失去 Resume 中斷的解(我們 sim 短,不需要) |
| `--full` | 只 `.dat` + `solve.out` + `magsolv.out` | 所有 intermediates 含 `.db` + `.rmg` | 16-20 GB → **145 MB** | **失去 A + B**(GUI Resume + 重 postproc) |

**預設一律 half-clean。`--full` 要使用者明確同意才可用**。

---

## 🗄 歸檔資料夾保留原則(archived,使用者 2026-06-15 拍板,比 half-clean 更嚴)

當資料夾是**已整理歸檔的 `magnetic_sim/ANSYS/main/ANSYS_data/` 結構**(coil1..6 / coilN/gap變體 / mesh/ / geom/ /
各 case_tag / MATLAB_data/),保留原則改為**白名單**:**只留交付檔,其餘 ANSYS 過度檔一律刪**。

| | 副檔名 | 說明 |
|---|---|---|
| **保留(交付檔)** | `.dat` `.db` `.cdb` `.mat` `.csv` `.npz` `.md` `.iges` | 場結果 / mesh·模型 / MATLAB 資料 / 說明 / 幾何 |
| **刪除(過度檔)** | `.rmg .full .esav .rst .rth .emat`(重求解中間檔,最大宗)、`.out .err .log .bat .BAT .stat .DSP*`(DSP/DSPtriU/...)`.page* .lock .mntr .ans .inp .dbb .rdb .dbe .png` + `scratch` | 全部 ANSYS 副產物 |

- **跟 half-clean 的差異**:此模式**連主 `.rmg` + log(solve.out/magsolv.out)都刪**。
  代價 = 失去「不 re-solve 重抽新物理量(主 .rmg)」與「Resume 中斷解」;
  **保留** = `.db` 重開模型/mesh + `.dat` 既有場 + `.mat` 分析結果。我們 sim 短,可接受。
- **適用**:整理 / 歸檔 `magnetic_sim/ANSYS/main/ANSYS_data/` 既有結構時(本次 2026-06-15 重組即此模式)。
- **不適用**:剛跑完、可能要 Resume 或重 postproc 新物理量的「活躍」sim → 仍走預設 half-clean。
- **安全性**:白名單含 `.mat/.csv/.npz` → `MATLAB_data/`、`all_coil/` 的資料天然不受影響;
  whitelist 清理可直接 walk 整個 `magnetic_sim/ANSYS/main/ANSYS_data/`,免特別排除。
- **對應**:上面「不可影響的 6 項工作」仍成立(MATLAB 全讀 `.dat`+`.mat`,過度檔本就不讀)。

---

## 強制流程(no exceptions)

### Step 1:讀此規則 + 對應 memory

```
Read .claude/rules/sim-cleanup.md
Read memory feedback_ansys_sim_cleanup_sop.md
```

### Step 2:Dry-run 列清單

任何 rm 前**必先 dry-run**,輸出每個目標 dir + 大小 + 模式(half / full)+ 預期釋出。
用 `find ... -type f \( ... \)`(不加 `-delete`)列檔案。

### Step 3:Sanity check 每 dir 有 `.dat`

無 `.dat` 的 dir **拒絕清**(可能 sim 沒完成 / 沒跑 postproc)。
使用者要 force 才能 override(`--force`)。

### Step 4:使用者**明確批准**才執行

不可自動 rm。dry-run 報告後等使用者說「OK / 執行」。
若使用者**未指定模式**,**預設 half-clean**,清前再說一句「half-clean,保 .db + 主 .rmg(GUI + 重 postproc 都還能用)」。

### Step 5:依保留清單手動清理（helper 已移除）

舊的 `clean_sim_dir.sh` helper 已於 2026-06-15 刪除（apdl 只留 .txt）。現在清理**手動**進行,
嚴格照本規則的保留/刪除清單(見上「清理模式對照」與「🗄 歸檔資料夾保留原則」):
- **歸檔模式(預設)**:只留白名單交付檔(`.dat .db .cdb .mat .csv .npz .md .iges`),其餘 ANSYS 過度檔全刪。
- 用 Python walk(可靠、跨平台)或逐 dir `rm`,**刪前先 dry-run 列清單**(Step 2)、**每 dir 確認有 `.dat`**(Step 3)、**使用者明確批准**(Step 4)。

不要憑記憶亂 `rm -rf` ANSYS 檔 — 一律照本規則白名單。

### Step 6:Report 釋出 + df

清完報告 before / after / saved + `df -h /g`。

---

## 砍整個 result dir(rm -rf entire dir)的額外規則

不只是清 intermediates,而是要**砍整個 dir**時(例如歷史 broken sim、被取代的舊版本):

1. 每個要砍的 dir 必須有 **memory 或 .claude rule 明確標記為 obsolete / deprecated / broken**(例如 [[long-fei-b-bar-matrix-v4]] 標 v3 obsolete)
2. **不可**根據檔名猜測(例如「`_BAD_`」聽起來像廢的就砍)— 必須有書面證據
3. 砍前必須 `grep -rl <dirname> analysis/ apdl/` 確認**沒有 active script 引用**
4. 整個 dir rm 後 commit 一筆紀錄(若有 git tracking)

---

## 不可清的清單(永遠不要砍)

不論任何模式、任何理由,以下**絕對禁止**清:

- 任何 `.dat` 檔
- `solve.out`、`magsolv.out`、`post.out` 等 log
- `<design>/data/<topic>/*.mat`、`*.csv` — 分析結果
- `<design>/figures/` 全部 — 圖檔
- `<design>/analysis/` 全部 — MATLAB 腳本
- `<design>/apdl/` 全部 — APDL source(geom + sim + postproc)
- `<design>/IGES*/` 全部 — 幾何匯出
- `<design>/CAD/` 全部 — SolidWorks 原檔
- `<design>/doc/` 全部 — LaTeX、PDF
- `<design>/comsol/`、`mph/`、`semulator/`、`reference/` 全部
- 任何 `.git/` 內容

---

## 觸發片語(以下任一即啟動此規則)

- 「清理 sim 副產物」/「清 ANSYS 結果」/「清 results」
- 「磁碟滿了」/「G: 滿了」/「D: 滿了」/「整理磁碟」
- 「清掉沒用的 sim」/「砍舊 sim」/「清歷史檔」
- 「跑 cleanup」/「整理 sim dir」
- 「rm -rf 結果資料夾」(危險 — 強制走 SOP)

## 何時不適用

此規則只管 **sim intermediate / result dir** 清理。**不**管:
- 編譯產物如 `*.aux`、`*.toc`、`__pycache__/`(屬一般軟體 cleanup)
- 使用者明確要求砍**某具體檔**(非清理場景)
