# RESULTS_MAP — long2016_hexapole_halfcut

**Single source of truth：這個 topic 底下每個 result dir 是什麼、能不能信、載入後該看到什麼指紋。**

讀任何結果前**先查這張表**（規則：`.claude/rules/result-read-safety.md`）。
本表凌駕 memory —— 若 memory 舊指標與本表衝突，以本表為準並更新 memory。
最後核對：2026-06-12（用 `wc -l` + `import_ansys_data` + bsum 欄掃描，read-only）。

---

## 怎麼用（層②核指紋）

`import_ansys_data(<dir>, <dataset>, <coilN>)` 讀
`<dir>/<coilN>_coord_<dataset>.dat` + `<coilN>_bfield_<dataset>.dat`，並印
`Matched %d nodes`。載入後**對照下表的「matched 節點」與「|B| max」**：

- 節點數可區分 **mesh 種類**（baseline vs graded）。
- 節點數**無法**區分 baseline vs `gap200um_mueq`（同 mesh）→ **必須再看 |B| max**：
  gap 版整體低約 30%。

`dataset` 也是指紋的一部分：同一個 baseline 幾何有 `all`（全域）與 `wp`（WP 區域子集）兩種匯出，節點數不同，**回報時要連 dataset 一起講**。

---

## CANONICAL（可信，正式分析用）

| dir | coil/dataset | 物理意義 | matched 節點 (wp) | `_all` 行數 | WP 區 \|B\| max 指紋 |
|---|---|---|---|---|---|
| `coil1` | coil1, wp/all | **Long2016 verbatim baseline**（v4 重建，下極半切）P1 激發 | 390579（air 341428） | 494873 | ~1.02 T |
| `coil2` | coil2 | baseline，P3 激發 | 390579 | 494873 | ~1.14 T |
| `coil3` | coil3 | baseline，P6 激發 | 390579 | 494873 | ~1.14 T |
| `coil4` | coil4 | baseline，P5 激發 | 390579 | 494873 | ~0.75 T |
| `coil5` | coil5 | baseline，P2 激發 | 390579 | 494873 | ~0.45 T |
| `coil6` | coil6 | baseline，P4 激發 | 390579 | 494873 | ~0.74 T |

> apdl coil 索引 ↔ 紙上極名：{1,2,3,4,5,6} = {P1,P3,P6,P5,P2,P4}（見 notation-glossary）。
> **2026-06-12 的泛化誤差分析、KI trend、B̄ matrix v4 baseline 用的就是這 6 個 dir 的 `wp` dataset。**

---

## VARIANT（可信，但**不是 baseline**，別跟上面混）

| dir | 物理意義 | matched 節點 | \|B\| max 指紋（判別關鍵） | 用途 |
|---|---|---|---|---|
| `coil1_gap200um_mueq` … `coil6_gap200um_mueq` | **μ_r 等效 200µm 氣隙**（不切幾何，protrusion 改 μ_r=31）| 與 baseline **相同**（~390579/494873）| **比 baseline 低約 30%**（coil1 ~0.71 / coil2 ~0.83 / coil5 ~0.25 T）| B̄ matrix v4 的 gap 對照 |
| `P1_graded` | graded 密 mesh，**P1 only**，charge-fit 用 | 與 baseline **不同**（graded mesh）| ~1.20 T | KI 電荷擬合密網格 |
| `P2toP6_graded` | graded 密 mesh，P2–P6 | graded（不同）| ~1.21 T | KI 電荷擬合密網格 |

⚠ **最易混點**：`coilN` vs `coilN_gap200um_mueq` 節點數一模一樣，**只有 \|B\| 能分**。要 baseline 卻載到 gap 版 → \|B\| 會低 ~30%，層②必須在此攔下。

---

## NON-RESULT（**沒有 .dat，不是模擬結果資料**，不可當結果讀）

幾何 / mesh 匯出與 log，**不含 `*_bfield_*.dat`**：

- `geom_export_metre`, `geom_export_metre_HollowProt`, `geom_export_metre_gap200um`, `geom_export_metre_sph3`
- `geom_export_mm`, `geom_export_mm_HollowProt`, `geom_export_mm_HollowProt_Plain`, `geom_export_mm_gap200um`, `geom_export_mm_sph3`
- `geom_mm_hp_split`, `geom_mm_withcoil`
- `mesh_dense`, `mesh_graded`, `mesh_iges`
- `*.log`（batch / resume log）

被要求「讀結果」卻指向以上任一 → **判定非結果資料，停下來問使用者**。

---

## OBSOLETE / BROKEN（歷史，不可用於分析）

| 對象 | 狀態 | 證據 |
|---|---|---|
| `TwBError` | broken（0 個 .dat，名稱即 error run）| 目錄掃描無 .dat |
| `coil1_pre_fine_mesh`（**已不存在**）| 歷史：舊 B̄ v1/v2 曾以此為「正確 halfcut coil1」，現已移除 | memory `project_long_fei_B_bar`（過時指標）；現行正確 coil1 = 上方 CANONICAL `coil1`（494873）|
| B̄ matrix v1 / v2 / v3 | 被 v4 取代 | memory [[long-fei-b-bar-matrix-v4]]、[[halfcut-missing-vadd-bug]] |
| smrt4 全域跑 halfcut | 不可靠（MAPDL crash）| memory [[smrt4-unreliable-halfcut]]；標準是 smrt5 + 0.6A |

> 關於 `coil1_pre_fine_mesh`：這正是「讀錯結果」的典型成因 —— memory 指向一個**已換版/已刪除**的 dir。現在以**本表**為準：baseline 正確 coil1 就是 `coil1`（494873 節點）。

---

## 維護
- 新增 / 重跑 / 刪除 result dir 時**同步更新本表**。
- 新 topic（full / kuo_quadrupole / zhang…）要各自建 `RESULTS_MAP.md`（本次僅 halfcut）。
