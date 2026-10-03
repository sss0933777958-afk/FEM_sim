# ANSYS `db/` 清理

## 三層政策（`ANSYS_data/<model>/db/`）

| 層 | 內容 | 政策 |
|---|---|---|
| `geom/` | 純幾何 | 整層刪 |
| `mesh/` | 幾何 + 網格 | 保主 `.db` + 主 log，其餘刪 |
| `sim/` | 幾何 + 網格 + 場 | 保主 `.db` + 主 `.rmg` + 主 log，其餘刪 |

## 刪 `geom/` 前必做

```bash
find <model>/db/geom -name "*.db" -size +100M -printf "%10s  %P\n" | sort -rn
```

有輸出 ⇒ 那是放錯層的網格 db，先搬進 `db/mesh/<case>/`，並同步改指向舊路徑的 APDL deck（`RESUME` / `/CWD` / `/OUTPUT`），搬完才刪其餘。

## 保留與刪除

- **保留**：`*.db`、主 `*.rmg`（僅 `sim/`）、主 run log（`solve.out` / `magsolv.out` / `<deck>.out`）。
- **刪除**：worker 的 `*.rmg` 與 log、`*.esav *.full *.DSP* *.emat *.rst *.rth`、`*.err *.stat *.lock *.page* *.bat *.tmp *.txt *.dbb *.ldhi`、`scratch`、`menust.tmp`。
- **搬走**：`*.png` → `figures/`。

**worker 判別**：stem 尾數去掉後，同夾存在該名的主檔，兩條件都成立才是 worker（`sim_vp0.rmg` 旁有 `sim_vp.rmg`）。只看「stem 以數字結尾」會把 `sim_coil1..6`、`mesh_lv1..3` 這類主檔誤判。判不準一律當主檔留著。

**rm 寫法**：不可寫 `rm -f <job>*`。用
```bash
rm -f <job>_*.rmg <job>*.esav <job>*.full <job>*.DSP* <job>*.page* <job>*.stat <job>*.err <job>*.ldhi
```

**GUI**：互動 MAPDL 一律用獨立 scratch 夾 `db/{mesh,geom}/_gui_view_<tag>/` 當 `-dir`；使用者關掉 GUI 後整夾刪（先確認 canonical db 的時間戳與大小沒變）。

## 絕對不可清

`.dat`、`.mat`、`.csv`、`.npz`、`.cdb`；`matlab/`、`figures/`、`apdl/`、`CAD_model/`、`model_check/`、`reference/`；`.git/`。

## 流程

1. Dry-run：列每個目標夾的將刪 / 將留、大小、預期釋出。
2. 確認刪完每個 `mesh/`、`sim/` 夾仍有 ≥1 顆 `.db`，且對應 `.dat` 已抽出。
3. 等使用者明確批准才執行。
4. 回報 before / after / 釋出量。

還在跑或可能要 Resume 的 sim 不清。
