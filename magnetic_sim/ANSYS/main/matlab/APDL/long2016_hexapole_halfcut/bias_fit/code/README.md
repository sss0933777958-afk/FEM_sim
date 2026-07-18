# long2016_hexapole_halfcut/bias_fit/code/ — bias_fit 程式碼層

**用途**：bias_fit 功能組的程式碼容器，依性質分 `scripts/`（計算）與 `plot/`（繪圖）。

**內容**：
- `scripts/` — fit / 半徑掃描 / 模型比較 / 電荷位置診斷（無圖輸出，產 `.mat` + `.tex`）。
- `plot/` — 收斂、參數 vs R、電荷動畫等繪圖腳本（讀 `.mat` 產 `.png`/`.gif`）。

**資料來源 / 流向**：見上層 `../README.md`（讀 `ansys_path` 的 FEM `.dat` → 算 → `matlab_path` 的 `.mat` / `results` 的 `.tex` / `figures` 的圖）。

**命名 / 慣例**：`code/{scripts,plot}` schema；計算與繪圖分流，不混放。

**相關**：見上層 `../README.md`、`../../common/README.md`、`../../../../CLAUDE.md`。
