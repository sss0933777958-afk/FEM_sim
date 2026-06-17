# magnetic_sim/ANSYS/backup/hung/analysis/fit — 電荷模型擬合 pipeline（charge-model fitting）（歸檔/archived）

> **本設計已歸檔（非活躍），保留供參考。** Hung hexapole（build-workflow design）。目前唯一活躍設計為 `../../../main/`。

**用途**：把 FEM 場配進點電荷模型的主擬合流程（[pre-J] ell → [J] joint → [B-6x]）。

**內容**（代表檔）：
- `fit_ell_percoil.m` / `fit_ell_perlayer.m` — 步驟 1/2，per-coil、per-layer 的 |ℓ| 與 sphere `pos` 初值。
- `fit_J.m`（主結果，err ~0.46%）/ `fit_J_50um.m` / `fit_J_l250*.m` — [J] joint 6-coil ideal K_I 擬合。
- `fit_B6x_1C.m` / `fit_B6x_6C.m` / `fit_B6x_allcoil.m` — [B-6x] 18 pos + C 係數版本。
- `fit_J.m`、`sweep_J_validity_hung.m` — charge-fit + validity sweep。

**資料來源 / 流向**：讀 `../../results/coilN/<變體>/*.dat`（透過 `../core/import_ansys_data.m`）→ 輸出 `.mat` 至 `../../data/`（見 `data/README.md` 的 producer/consumer 表）。

**命名 / 慣例**：依序執行（後段腳本載入前段 `.mat`）；`_l250` = ℓ=250 µm 變體。擬合電流須對齊 FEM 激發電流。

**相關**：見 `../README.md`（執行順序與檔名對照）、`../../data/README.md`、上層 `../../README.md`。
