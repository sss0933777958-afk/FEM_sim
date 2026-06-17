# magnetic_sim/ANSYS/backup/hung/analysis/variants — 替代 / 舊版擬合腳本（alternative scripts）（歸檔/archived）

> **本設計已歸檔（非活躍），保留供參考。** Hung hexapole（build-workflow design）。目前唯一活躍設計為 `../../../main/`。

**用途**：非當前主 pipeline 的替代 / 前一版擬合方法；路徑已調整、保持可重跑，供日後與主流程對照或探索替代法。

**內容**：
- `fit_KI_v1.m` — v1（用 method [A] 的單一 ell）；主流程已改用 `../fit/fit_ell_perlayer.m`（per-layer ell）。
- `fit_J_fittedKI.m` — 用**fitted** K_I；主流程 `../fit/fit_J.m` 用**ideal** K_I（0.46% err，已確認足夠）。

**資料來源 / 流向**：讀 `../../results/.../*.dat` → 輸出 `.mat` 至 `../../data/variants/`（見 `data/README.md`）。

**命名 / 慣例**：`variants` = 變體腳本，非主流程；對應 `.mat` 收在 `../../data/variants/`。

**相關**：見 `../README.md`（main vs variant 對照）、`../../data/README.md`、上層 `../../README.md`。
