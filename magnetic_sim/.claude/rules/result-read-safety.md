# 讀 ANSYS 結果防呆規則（強制讀取）

**這份規則解決一個問題：使用者要讀 A 結果，Claude 卻讀成 B 結果。**

當任何工作涉及「**載入 / 讀取 ANSYS 模擬結果**」時 —— 抽 `.dat`、`import_ansys_data`、載入 `coilN`、跑 postproc 前載入、算矩陣 / fit / 畫場圖前載入 result —— **動手前必須先讀完此規則全文**，並照三層流程執行。違反即視為違規。

對應 memory：`feedback_result_read_safety.md`
對應清單（single source of truth）：各 topic 的 `magnetic_sim/ANSYS/main/ANSYS_data/<topic>/RESULTS_MAP.md`
（目前已建：`magnetic_sim/ANSYS/main/ANSYS_data/long2016_hexapole_halfcut/RESULTS_MAP.md`）

相關規則：`.claude/rules/iges-model-id.md`（識別物理模型）、`.claude/rules/main-workspace.md`、`.claude/rules/sim-cleanup.md`
相關 memory：[[verify-before-act]]、[[long-fei-b-bar-matrix-v4]]、[[smrt4-unreliable-halfcut]]

---

## 為什麼會讀錯（根因，先認清）

1. **資料夾名字會騙人 / 太相近**：`coil1` vs `coil1_gap200um_mueq` vs `coil1_pre_fine_mesh`、`geom_*`、`mesh_*`、`smrt4/smrt5`。
2. **憑記憶猜路徑**：看到「讀 A」就自己推一個 dir，沒先回報、沒讓使用者確認。
3. **載入後不核指紋**：沒檢查「這份資料的特徵是否符合預期」，於是讀到 uncut / gap / graded / obsolete 也照算。
4. **memory 指標會過時**：真相換版（v1→v4）後，舊 memory 仍指向已刪除 / 已取代的 dir。

---

## 🔒 三層防呆（每次讀結果都要走完）

### 層 ① 讀之前：先回報 + 消歧（不自己猜）

動手讀任何結果前，**先輸出一行**：

> 「我要讀：`<絕對路徑>` 的 `<dataset>`，它是 `<topic>` 的 `<物理意義>`，期望指紋 ≈ matched `<N>` 節點 / \|B\| max `<值>`。」

- 路徑、dataset、物理意義**先查 `RESULTS_MAP.md`**，不靠記憶。
- 只要有 **≥2 個合理候選**（例如要 baseline 但同時存在 `coilN` 與 `coilN_gap200um_mueq`），**列出候選請使用者選，不自己挑一個**。
- 候選若指向 **NON-RESULT**（`geom_*` / `mesh_*` / `.log`，無 `*_bfield_*.dat`）→ 判定「非結果資料」，停下來問。

### 層 ② 讀之後：核指紋（多訊號，對不上就停）

`import_ansys_data` 會印 `Matched %d nodes`。載入後**立刻比對**下列三訊號與 `RESULTS_MAP.md` 期望值：

| 訊號 | 來源 | 能分辨 |
|---|---|---|
| **matched 節點數** | import 印的 `Matched N nodes` | mesh 種類（baseline 390579/wp vs graded） |
| **\|B\| max（air 節點）** | 算 `max(bsum)` | **baseline vs gap_mueq**（gap 低 ~30%，節點數相同分不出）|
| **路徑 case_tag** | dir 名字串 | 變體 / obsolete |

⚠ **節點數相同 ≠ 同一份**：`coilN` 與 `coilN_gap200um_mueq` 節點數一樣，**必須再看 \|B\| max** 才能確認。
任一訊號對不上預期 → **立刻停、回報差異、不繼續計算 / 畫圖 / fit**。

### 層 ③ 全程：以 RESULTS_MAP 為準（凌駕 memory）

- 「哪個 dir 是對的」一律查 `RESULTS_MAP.md`，**不靠 memory 的舊指標**。
- 發現 memory 與清單衝突（例如 memory 指向已刪除的 `coil1_pre_fine_mesh`）→ 以清單為準，並**更新該 memory**。
- 清單沒涵蓋的 topic / 新 dir → **先 inspect 補進清單**再讀，不裸讀。

---

## 強制規則

1. 讀結果前**先查 `RESULTS_MAP.md`**，再走層①回報 + 消歧。
2. **≥2 候選不自己猜**，列給使用者選。
3. 載入後**必核 matched 節點數 + \|B\| max + case_tag** 三訊號；對不上**停**。
4. 報告載入時**連 dataset（all / wp / circuit）一起講**，不只報一個節點數。
5. 清單未涵蓋者**先 inspect 補表**，不裸讀。
6. memory 與清單衝突 → **清單為準 + 更新 memory**。

---

## 已知陷阱速查（long2016_hexapole_halfcut）

- **baseline 正確 coil1–6** = `coilN`，`all`=494873 行 / `wp`=390579 matched（air 341428），\|B\| max ~1.0–1.14 T。
- **`coilN_gap200um_mueq`** = μ_r 等效氣隙，**同節點數**，\|B\| max **低 ~30%**（唯一判別靠 \|B\|）。
- **`P1_graded` / `P2toP6_graded`** = charge-fit 密 mesh，節點數**不同**。
- **`geom_*` / `mesh_*` / `*.log` / `TwBError`** = 無 `.dat`，**不是結果**。
- **`coil1_pre_fine_mesh`** = 已刪除的舊「正確 coil1」；現以 `coil1`（494873）為準。

---

## 觸發片語（任一即啟動此規則）

- 「讀結果」/「載入 coilN」/「抽 .dat」/「import_ansys_data」
- 「跑 postproc」/「算 B 矩陣 / fit / 畫場圖」前需載入 result
- 「讀 baseline / gap / graded」/「讀那個模擬」/「讀 A 結果」

## 何時不適用

- 不載入 result 的純幾何 / mesh / IGES / CAD 工作。
- 已在記憶體中的資料二次處理（同一次對話已核過指紋、未換 dir）。
- COMSOL / 非 ANSYS 結果（COMSOL 有自己的座標 / 單位坑，見 `.claude/rules/comsol-livelink.md`）。
