# 規則 #1：`ANSYS_data/<model>/db/` 只留 `.db` + 主 `.rmg` + 主 log

**使用者拍板（2026-06-26；2026-07-06 加主 log）**：`ANSYS_data/<model>/db/` 底下**每個子資料夾**只保留 **`.db`（模型/網格）+ 主 `.rmg`（無 digit 結尾的那顆）+ 主 run log（`solve.out`/`magsolv.out`/`<deck>.out`，無 digit 結尾）**；其餘 ANSYS 殘留檔**一律刪、禁止留**；**render/檢視 `.png` 搬到 `figures/<model>/<activity>/`、不留 db/**。

> 此為 db/ 資料夾的**保留白名單**。2026-07-06 更新後**與全域 `…/FEM_sim/.claude/rules/sim-cleanup.md` half-clean 一致**（都保 `.db` + 主 `.rmg` + solve/magsolv log）；差別只在本規則另清 per-worker log + 散落中間檔、且圖歸 figures/。保 `.db`+主`.rmg` = GUI resume + 不重解重抽物理量；保主 log = 節點數/警告/收斂 audit。

## 🔒 db/ 先分 `geom/`、`mesh/`、`sim/` 三層（使用者拍板 2026-07-15；2026-07-15 加 mesh/ 層）
`ANSYS_data/<model>/db/` 底下**先分三層再放 `<case>/`**：`db/geom/<case>/`、`db/mesh/<case>/`、`db/sim/<case>/`。
**依 db 性質分**（呼應「純幾何 db vs 已 mesh db（網格種子）vs 解完 db（含場、可 RESUME 重抽資料）是不同東西」）：
- **`db/geom/<case>/`** — 純幾何 db（不含網格、不含場）：`geom*.db`、`exp*.db`、`gck*.db`、`from_parasolid*.db`、
  `*.iges`（幾何匯出）+ 對應主 log（`geom*.out`/`run.out`，無 digit 結尾）。
- **`db/mesh/<case>/`** — 已 mesh db（幾何+網格、不含場）：`mesh*.db`、`<case>.db`（mesh-only deck 產）+ 主 log
  （`mesh*.out`/`meshfull.out`，無 digit 結尾）。**（2026-07-15 使用者拍板新增；mesh db 從 geom/ 改放此層）**。
  ⚠ 既有 `db/geom/mesh_graded*` 為舊置（暫不搬、待使用者定）；**新的 mesh-only 產物一律進 `db/mesh/`**。
- **`db/sim/<case>/`** — 解完 db（含場）：`sim*.db`、**主 `.rmg`**（無 digit 結尾）、`gap_sweep*`、`coil2to6*` +
  主 log（`magsolv.out`/`solve.out`）。
- **混合 case**（同 case 有 mesh 幾何 db 又有 solve 結果 db+rmg，如 `singlepole/<shape>/`）→ **依檔案拆兩邊**
  （`mesh_*.db`→geom、`sim_*.db+.rmg`→sim，保留 `<shape>/` 子結構）。
- **空夾一律刪**。每 case 內仍照下方 half-clean 白名單（清 per-worker + 中間檔）。
- **deck 連動**：`/CWD`、`/OUTPUT`、`RESUME` 路徑須指 `db/geom/<case>` 或 `db/sim/<case>`（依該 case 性質；
  一個 sim deck 常同時 RESUME `db/geom/<mesh_case>`（種子）+ 寫 `db/sim/<solve_case>`（結果））。

## 🔒 保留 / 刪除清單

| | 副檔名 | 說明 |
|---|---|---|
| **保留** | `*.db` | 模型 + 網格（GUI `RESUME`、重解都靠它） |
| **保留** | 主 `*.rmg`（**stem 不以 digit 結尾**，如 `sim_singlepole.rmg`、`coil2to6_gap200.rmg`） | 不重解就重抽新物理量（H/energy）的結果庫 |
| **保留** | 主 run log：`solve.out` / `magsolv.out` / `<deck>.out`（如 `mesh_filled.out`；**stem 不以 digit 結尾**） | 節點數/警告/收斂 audit（2026-07-06 使用者拍板加入；與 sim-cleanup half-clean 一致） |
| **刪除** | per-worker `*.rmg`（stem 以 digit 結尾，如 `sim_singlepole0.rmg`、`*_0.rmg`） | DMP 分散 worker 副本 |
| **刪除** | `*.esav` `*.full` `*.DSP*` | 最大宗中間檔 |
| **刪除** | **per-worker log**（stem 以 digit 結尾：`*0.out`/`*1.out`/`*0.err`/`*0.log`）、`*.err` `*.stat` `*.lock` `*.page*` `*.bat` `*.tmp` `*.txt` `scratch` `menust.tmp` | worker log / 鎖 / 暫存 / 散落 deck（**主 log 已保，見上**） |
| **搬移** | `*.png`（render / 檢視圖） → **搬 `figures/<model>/<activity>/`，不留 db/** | 圖歸 figures/（figure-output 規則）；db/ 不放圖 |

**「主 .rmg」判別**：stem（去 `.rmg`）**不以數字結尾**＝主；以數字結尾（含 `_0`、或無底線的 `0/1/2/3`）＝ worker，刪。
⚠ 注意命名陷阱：`coil2to6_gap200.rmg` 的 stem 以 `200` 結尾但**是主檔**（gap200 是 case 名，不是 worker index）——別被尾數字誤判。worker 一定是「同夾存在去掉尾數的同名主檔」那一組 `0/1/2/3`。

## ⚠ rm 指令陷阱（2026-06-29 踩過：誤刪主 .rmg）
**用 `rm` 腳本清理時，絕對不要寫 `rm -f <jobname>*`**——那會把主 `.rmg` + `.db` 一起刪掉（變成 full-clean、失去免重解能力）。曾因 `rm -f gap200_sb_P$k*` 把整顆 .rmg/.db 刪了 → 之後要元素 B 被迫重解。
**正解**＝用**針對性 pattern 只刪該刪的**，保住主 `.rmg` + `.db`：
```
rm -f <job>_*.rmg  <job>*.esav  <job>*.full  <job>*.DSP*  <job>*.page*  <job>*.stat  <job>*.err  <job>*.ldhi
# 留 <job>.db、<job>.rmg（主）、<job>*.out、*.dat
```
（`<job>_*.rmg` 只中 worker `_0/_1…`，不中主 `<job>.rmg`。）

## 🖥 GUI scratch 夾：關掉 GUI 就整夾刪（使用者拍板 2026-07-17）

開互動 MAPDL GUI（`-g`）檢視 mesh/model 一律用**獨立 scratch 夾**當 `-dir`（`db/mesh/_gui_view_<tag>/`、`db/geom/_gui_view_<tag>/`），
**絕不可用真正的 mesh/geom 夾**——GUI 關閉時使用者按 Save，MAPDL 會把當下畫面模型寫回 `-dir`，用真夾就**覆蓋掉 canonical db**
（hung R700/R300 兩顆 `mesh_graded.db` 被覆蓋成 steel-only，只能重跑還原）。

**使用者關掉 GUI 後，Claude 必須主動把整個 scratch 夾刪掉**（不必問、不留任何一檔）：

```
rm -rf <db/mesh|geom>/_gui_view_<tag>/
```

- **整夾刪、不是只清殘留**：夾內全部都是拋棄式副本 —— `view.db`（GUI Save 覆蓋過的畫面模型）、`view.dbb`（覆蓋前備份）、
  `!T!.BAT`、`view[0-3].{err,out,log,page}`、`cleanup-ansys-*.bat`。**沒有一個是交付物**（canonical db 在真 mesh 夾、原封未動）。
- **一顆 view.db ≈ 4 GB、加 .dbb 就 8 GB+**，留著純浪費（本次 `_gui_view_sleeve` 佔 8.4 GB）。
- 刪之前**先核 canonical db 的時間戳/大小未變**（證明 GUI Save 沒寫到真夾），再刪 scratch、回報釋出量。
- ⚠ 唯一例外：**若 canonical db 真的被覆蓋了**（時間戳變了、`VOLU/節點數`對不上），`view.dbb` 或真夾的 `.dbb` 是覆蓋前的備份 —— **先救援再刪**。

對應 memory：`feedback_gui_save_overwrites_db.md`。

## 強制流程
1. **動手前先 dry-run** 列出每個子夾「將刪 / 將留」分類（照上表），確認每夾刪完仍 **≥1 個 `.db`**（否則該夾不能 GUI resume，停下問）。
2. 使用者**明確批准**才刪（或本規則已是常駐授權的例行清理時，dry-run 後逕行並回報）。
3. 清完報告 before/after/freed + `df -h /g`。

## 觸發片語（任一即套用）
- 「清 db」/「整理 db 資料夾」/「db 只留 db 檔」
- 跑完 sim 後歸檔 `.db` / mesh、要清 db/ 夾殘留時
- 「db/ 怎麼留」/「db 殘留檔」
- **使用者關掉 GUI（`-g` 檢視結束）→ 自動觸發上面「GUI scratch 夾整夾刪」**

## 何時不適用
- `ANSYS_data/<model>/data/`（FEM `.dat` 場 + `.mat` 等**交付檔**白名單）——那是另一條（見 `matlab-output-layout.md` 與全域 sim-cleanup「歸檔資料夾保留原則」）。
- 還在跑 / 可能要 Resume 中斷解的「活躍」sim → 暫不清，跑完再套。

相關：全域 `…/FEM_sim/.claude/rules/sim-cleanup.md`、memory `feedback_ansys_sim_cleanup_sop`、`feedback_matlab_local_data_layout`。
