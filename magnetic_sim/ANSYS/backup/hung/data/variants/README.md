# magnetic_sim/ANSYS/backup/hung/data/variants — 替代 / 舊版擬合資料（alternative .mat）（歸檔/archived）

> **本設計已歸檔（非活躍），保留供參考。** Hung hexapole（build-workflow design）。目前唯一活躍設計為 `../../../main/`。

**用途**：非當前主 pipeline 的擬合結果 `.mat`（由 `../../analysis/variants/` 產出或 legacy orphan）。

**內容**：
- `joint_6coil_fit.mat` — `../../analysis/variants/fit_J_fittedKI.m` 產出（[J] fitted K_I，已被 ideal K_I 的 `fit_J.m` 取代）。
- `B6x_hung_fit.mat` — **orphan**，`fit_B6x_hung.m` 拆成 1C/6C 前的殘留，已不產出 / 不被引用。
- `wp_fitting_data.mat` — **orphan**，repo 無腳本引用，來源不明，保留供安全。

**資料來源 / 流向**：上游 `../../analysis/variants/`（讀 `../../results/.../*.dat`）→ 本層 `.mat`。

**命名 / 慣例**：對應主 pipeline 的 `.mat` 在上一層 `../`（`data/` root）；本層只放變體 / legacy。

**相關**：見 `../README.md`（producer/consumer 表與 active pipeline）、`../../analysis/variants/README.md`、上層 `../../README.md`。
