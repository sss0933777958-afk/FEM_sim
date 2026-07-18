# long2016_hexapole_halfcut/bias_fit/code/plot/ — bias_fit 繪圖腳本

**用途**：no_fix_l bias 模型的視覺化腳本（讀 `.mat` 產 `.png`/`.gif`）。

**內容**：
- `plot_nofixl_convergence.m` — fit cost（殘差平方和）vs lsqnonlin 疊代數，每半徑一條曲線。
- `plot_nofixl_params_vs_R.m` — `ℓ̂` / `gB` / `e_hat` 等參數隨取樣半徑 R 的趨勢圖。
- `anim_charge_side.m` — 等效電荷側視動畫（P1/P2 隨 R 變化的 `.gif` + 代表幀 `.png`）。

**資料來源 / 流向**：讀 `matlab_path(model,'charge_fit',...)` 的 `.mat` → 圖存 `../../figures/` 與 `main/figures/long2016_hexapole_halfcut/`。

**命名 / 慣例**：⚠ **新增繪圖腳本前須先依 `../../../../../CLAUDE.md`「繪圖腳本規則」**：先確認屬哪個功能組、一任務一腳本原地改、定案後才存最終圖（定案前用 MCP preview）；場圖一律真實 FEM 節點原值、不內插。

**相關**：見 `../../README.md`、`../../../common/README.md`、`../../../../../CLAUDE.md`（繪圖規則 + Figure Production）。
