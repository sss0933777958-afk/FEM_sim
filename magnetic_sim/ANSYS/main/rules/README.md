# `main/rules/` — main/ 本地規則

`magnetic_sim/ANSYS/main/` 專屬的資料夾佈局/保留/工作規則（使用者 2026-06-26 拍板）。
**每次在 main/ 開工前先讀本資料夾全部規則**（規則 #4）。動到對應路徑前讀該規則全文。
全域規則仍見 `…/FEM_sim/.claude/rules/`；這裡放只綁 main/ 的局部規則。

| 規則檔 | 一句話 | 綁的路徑 / 時機 |
|---|---|---|
| `read-rules-first.md`（#4） | **每次開工先讀本資料夾全部規則** | main/ 下任何工作 |
| `db-folder-retention.md`（#1） | `db/` 子夾只留 `.db` + 主 `.rmg`，殘留禁留 | `ANSYS_data/<model>/db/**` |
| `matlab-output-layout.md`（#2） | `.mat` 放產生它的程式旁 `data/`（`MATLAB_data/` 已移除） | `matlab/<model>/<activity>/data/` |
| `results-pdf-only.md`（#3） | 該 `results/` 只放 PDF | `…/Hall_sensor_base_fix_dir/results/` |
| `figure-style.md`（繪圖風格） | 圖表視覺風格 preset 目錄；**畫圖前必先問使用者要哪個選項**。選項①＝粗體框圖（大字粗體/box/無 grid/tick 減半/單位 `()`）；**3D 框體＝手動 `draw_box_edges`（12 邊省略最遠角 3 條）+ box off + `daspect([1 1 1])` + 三軸同刻度（不要 axis equal／BoxStyle full／單純 box on）**；數值標註：**10^0 不標、無單位不標**（圖+結果 PDF 通用） | 畫任何圖之前 |
| `figure-output.md`（繪圖輸出） | 圖一律**輸出實檔**到 `figures/`；要改**原地改腳本→重跑→覆蓋同檔**，迭代到定案（不丟 temp 等定案） | 畫任何圖 |

相關 memory：`feedback_matlab_local_data_layout`、`feedback_ansys_sim_cleanup_sop`、`reference_local_latex_compile`、`feedback_field_quiver_style`。
