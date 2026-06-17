# long2016_hexapole_halfcut/fixl_fit/code/plot — 擬合結果繪圖腳本

**用途**：把 fixl_fit 的擬合過程 / 結果畫成圖（收斂曲線、參數 vs 取樣半徑 R）。

**內容**：`plot_KI_convergence.m`、`plot_fixl_convergence.m`（cost 收斂）、`plot_fixl_params_vs_R.m`（ℓ̂/ĝ_B/NRMSE vs R）。

**資料來源 / 流向**：讀擬合產出（MATLAB_data `.mat` / scripts 計算結果）→ 出 PNG 到 `../../figures/`。

**命名 / 慣例**：一張定案圖一支腳本；PNG 放 `../../figures/`。**新增繪圖腳本前須依 `../../../../CLAUDE.md` 繪圖腳本規則**：先確認功能組、一任務一腳本（原地改到定案、定案前不另開）、定案後才存圖、真實節點不內插。

**相關**：見 `../../README.md`、`../../../../CLAUDE.md`。
