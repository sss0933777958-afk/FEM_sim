# 規則 #1：`ANSYS_data/<model>/db/` 只留 `.db` + 主 `.rmg`

**使用者拍板（2026-06-26）**：`ANSYS_data/<model>/db/` 底下**每個子資料夾**只保留 **`.db`（模型/網格）+ 主 `.rmg`（無 digit 結尾的那顆）**；其餘 ANSYS 殘留檔**一律刪、禁止留**。

> 此為 db/ 資料夾的**保留白名單**，比全域 `…/FEM_sim/.claude/rules/sim-cleanup.md` 的 half-clean 再嚴一點（half-clean 還會留 log）。兩者一致處：**都保 `.db` + 主 `.rmg`**（GUI resume + 不重解就重抽物理量的能力都不失）。本規則只是把 db/ 夾的 log/中間檔也清掉。

## 🔒 保留 / 刪除清單

| | 副檔名 | 說明 |
|---|---|---|
| **保留** | `*.db` | 模型 + 網格（GUI `RESUME`、重解都靠它） |
| **保留** | 主 `*.rmg`（**stem 不以 digit 結尾**，如 `sim_singlepole.rmg`、`coil2to6_gap200.rmg`） | 不重解就重抽新物理量（H/energy）的結果庫 |
| **刪除** | per-worker `*.rmg`（stem 以 digit 結尾，如 `sim_singlepole0.rmg`、`*_0.rmg`） | DMP 分散 worker 副本 |
| **刪除** | `*.esav` `*.full` `*.DSP*` | 最大宗中間檔 |
| **刪除** | `*.out` `*.log` `*.err` `*.stat` `*.lock` `*.page*` `*.bat` `*.tmp` `*.txt` `scratch` `menust.tmp` | log / 鎖 / 暫存 / 散落 deck |

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

## 強制流程
1. **動手前先 dry-run** 列出每個子夾「將刪 / 將留」分類（照上表），確認每夾刪完仍 **≥1 個 `.db`**（否則該夾不能 GUI resume，停下問）。
2. 使用者**明確批准**才刪（或本規則已是常駐授權的例行清理時，dry-run 後逕行並回報）。
3. 清完報告 before/after/freed + `df -h /g`。

## 觸發片語（任一即套用）
- 「清 db」/「整理 db 資料夾」/「db 只留 db 檔」
- 跑完 sim 後歸檔 `.db` / mesh、要清 db/ 夾殘留時
- 「db/ 怎麼留」/「db 殘留檔」

## 何時不適用
- `ANSYS_data/<model>/data/`（FEM `.dat` 場 + `.mat` 等**交付檔**白名單）——那是另一條（見 `matlab-output-layout.md` 與全域 sim-cleanup「歸檔資料夾保留原則」）。
- 還在跑 / 可能要 Resume 中斷解的「活躍」sim → 暫不清，跑完再套。

相關：全域 `…/FEM_sim/.claude/rules/sim-cleanup.md`、memory `feedback_ansys_sim_cleanup_sop`、`feedback_matlab_local_data_layout`。
