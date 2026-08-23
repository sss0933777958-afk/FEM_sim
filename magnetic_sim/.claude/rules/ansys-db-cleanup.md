# ANSYS `db/` 清理規則（強制讀取）

**使用者拍板（2026-08-17）**：`ANSYS_data/<model>/db/` 三層依「重生成本」分級處理：

| 層 | 內容 | 政策 |
|---|---|---|
| `db/geom/` | 純幾何（點/線/面/體，無網格、無場）| **整層刪** —— 跑幾何 deck 幾分鐘可重生 |
| `db/mesh/` | 幾何 ＋ 網格，未解 | **保主 `.db` ＋ 主 log**，其餘刪 |
| `db/sim/` | 幾何 ＋ 網格 ＋ 場 | **保主 `.db` ＋ 主 `.rmg` ＋ 主 log**，其餘刪 |

為什麼這樣分：`mesh` 的**節點編號與座標在切網格那一刻定下來**，之後所有 `.dat` 都源自它 ——
網格重建 ⇒ 節點數變 ⇒ 舊 `.dat` 與新解**不可比**（`gap_200um` vs `no_gap` 同節點對照、
`result-read-safety` 的 494873 節點指紋都建立在此）。`sim` 的場已抽成 `.dat`，主 `.rmg` 只多買
「不重解就抽新物理量（H 場 / 能量 / 新 PATH）」。`geom` 完全由腳本決定，最便宜。

（本規則 2026-08-17 由 `sim-cleanup.md` ＋ `db-folder-retention.md` 合併，兩者已刪除。）

---

## ⚠ 刪 `geom/` 前的強制前置檢查（不可略過）

**`geom/` 這層的名字會騙人。** 2026-08-17 實測：三個 model 共 **20 顆網格 db（29.2 GB）被錯放在
`geom/` 底下**（`mesh_baseline/`、`mesh_graded*/`、`mesh_graded_basegap/`、`singlepole/`、
`mesh_graded_upper/`），其中 **19 顆是孤本**，含 `result-read-safety` 用來核指紋的基準網格。
整層刪會毀掉它們。

**動手前必跑**：

```bash
find <model>/db/geom -name "*.db" -size +100M -printf "%10s  %P\n" | sort -rn
```

- **有輸出 ⇒ 那是網格 db**（純幾何只有 4–22 MB，**大小本身就是判準**）
  → **先搬進 `db/mesh/<case>/`**，搬完才刪 `geom/` 剩餘。
- **無輸出 ⇒ 全是純幾何**，可整層刪。
- 搬移後**必須同步更新指向舊路徑的 APDL deck**（`RESUME` / `/CWD` / `/OUTPUT`）：
  `grep -rl "db[/\\]geom[/\\]<case>" apdl/`（2026-08-17 實測有 20 支 deck 指向 `db/geom/`）。

---

## 🔒 `mesh/` 與 `sim/` 的保留白名單

| | 檔案 | 說明 |
|---|---|---|
| **保留** | `*.db` | 模型 ＋ 網格（GUI `RESUME`、重解都靠它）|
| **保留** | 主 `*.rmg`（**僅 `sim/`**）| 不重解就抽新物理量的結果庫 |
| **保留** | 主 run log（`solve.out` / `magsolv.out` / `<deck>.out`）| 節點數 / 警告 / 收斂 audit |
| **刪除** | per-worker `*.rmg`、per-worker log（`*0.out` / `*1.err` …）| DMP 分散 worker 副本 |
| **刪除** | `*.esav` `*.full` `*.DSP*` `*.emat` `*.rst` `*.rth` | 重求解中間檔（最大宗）|
| **刪除** | `*.err *.stat *.lock *.page* *.bat *.tmp *.txt *.dbb *.ldhi` `scratch` `menust.tmp` | 鎖 / 暫存 / 散落 deck |
| **搬移** | `*.png`（render / 檢視圖）| → `figures/<model>/<activity>/`，**不留 db/** |

### 🔴 「主檔 vs worker」判別 —— 只有一個正確測試

**worker ⟺ 把 stem 的尾數去掉後，同夾存在該名的主檔。** 兩個條件都要成立，缺一不可。

```python
stem = os.path.splitext(fname)[0]
base = re.sub(r'_?\d+$', '', stem)
is_worker = re.search(r'\d$', stem) and base != stem and base in 同夾其他stem集合
```

⚠⚠ **絕對不可只用「stem 以數字結尾」當判準** —— 2026-08-17 實測：那會把
`sim_coil1.rmg` … `sim_coil6.rmg`（**六顆不同 coil，各自都是主檔**）全判成 worker，
一次誤刪 **163 GB** 不可重生的結果。同類受害者：`coil1_lv1/lv2/lv3`、`coil1_gap400`、
`coil2to6_gap200`（gap200 是 case 名不是 worker index）。

真 worker 長這樣：同夾同時有 `sim_vp.rmg`（主）＋ `sim_vp0/1/2/3.rmg`（DMP 副本）。

**同一測試也適用 log**（`.out` / `.log`）—— `mesh_lv1.out` / `mesh_lv2.out` / `mesh_lv3.out`
是三個 refine 等級的**主 log**，不是 worker。實測套錯判準會誤刪 104 個主 log。

**判別不準時一律當主檔留著** —— 留錯只是佔空間，刪錯不可逆。

### ⚠ `rm` 指令陷阱（2026-06-29 踩過：誤刪主 `.rmg`）

**絕對不要寫 `rm -f <jobname>*`** —— 會把主 `.rmg` ＋ `.db` 一起刪掉。正解是針對性 pattern：

```bash
rm -f <job>_*.rmg  <job>*.esav  <job>*.full  <job>*.DSP*  <job>*.page*  <job>*.stat  <job>*.err  <job>*.ldhi
# 留 <job>.db、<job>.rmg（主）、<job>*.out
```

（`<job>_*.rmg` 只中 worker `_0/_1…`，不中主 `<job>.rmg`。）

### 🖥 GUI scratch 夾：關掉 GUI 後整夾刪

互動 MAPDL GUI（`-g`）一律用獨立 scratch 夾當 `-dir`（`db/{mesh,geom}/_gui_view_<tag>/`），
**絕不可用真正的 mesh/geom 夾** —— GUI 關閉時使用者按 Save 會**覆蓋掉 canonical db**
（hung R700/R300 兩顆 `mesh_graded.db` 曾被覆蓋成 steel-only，只能重跑還原）。

使用者關掉 GUI 後，**Claude 主動把整個 scratch 夾刪掉**（不必問）：`rm -rf .../_gui_view_<tag>/`。
夾內全是拋棄式副本（`view.db` 一顆就 ~4 GB）。刪前先核 canonical db 的時間戳/大小未變；
**若真被覆蓋了，`view.dbb` 或真夾的 `.dbb` 是覆蓋前備份 —— 先救援再刪。**

---

## 🔒 絕對不可清（任何模式、任何理由）

- 任何 `.dat`（FEM 場，MATLAB 的交付物）、`.mat` `.csv` `.npz` `.cdb`
- `matlab/` 全部（腳本 ＋ `data/`）、`figures/` 全部、`apdl/` 全部
- `CAD_model/` `model_check/` `comsol/` `reference/`（原 `doc/` 已於 2026-08-23 併入 `reference/`）
- 任何 `.git/` 內容

## 強制流程

1. **Dry-run**：列每個目標夾的「將刪 / 將留」分類 ＋ 大小 ＋ 預期釋出。用 `find`（不加 `-delete`）。
2. **Sanity check**：刪完每個 `mesh/` `sim/` 夾仍須 **≥1 顆 `.db`**（否則無法 GUI resume）；
   刪 `sim/` 相關前確認對應 `.dat` 已抽出（`find <model>/data -name "*.dat" | wc -l`）。
3. **等使用者明確批准**才執行 —— 不可自動 `rm`。
4. **回報** before / after / freed ＋ `df -h /g`。

## 觸發片語

- 「清 db」/「整理 db 資料夾」/「清 sim 副產物」/「清 results」
- 「磁碟滿了」/「G: 滿了」/「整理磁碟」/「砍舊 sim」
- 跑完 sim 要歸檔、或**使用者關掉 GUI**（`-g` 檢視結束）時

## 何時不適用

- `ANSYS_data/<model>/data/`（`.dat` ＋ `.mat` 交付檔白名單）—— 那是交付物，見「絕對不可清」。
- 還在跑 / 可能要 Resume 中斷解的「活躍」sim → 跑完再套。
- 編譯產物（`*.aux`、`__pycache__/`）—— 一般軟體 cleanup，非本規則。

相關：`result-read-safety.md`（節點數指紋）、`matlab-output-layout.md`、`modify-existing-files.md`。
