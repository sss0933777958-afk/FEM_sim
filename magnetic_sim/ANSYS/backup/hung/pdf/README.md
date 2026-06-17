# magnetic_sim/ANSYS/backup/hung/pdf — 擬合報告 PDF（compiled reports）（歸檔/archived）

> **本設計已歸檔（非活躍），保留供參考。** Hung hexapole（build-workflow design）。目前唯一活躍設計為 `../../../main/`。

**用途**：擬合結果的編譯後報告 PDF（charge-model |ℓ| / Tip→Charge 表）。

**內容**：
- `fitting_l250.pdf` — ℓ=250 µm 擬合報告。
- `fitting_l500.pdf` — ℓ=500 µm 擬合報告。
- `scripts/` — 產生這些 PDF 的 LaTeX 原稿（見 `scripts/README.md`）。

**資料來源 / 流向**：`../analysis/fit/` 擬合 `.mat` → `scripts/*.tex`（手動填表）→ 編譯為本層 `.pdf`（pdf 為報告產物）。

**命名 / 慣例**：`l250` / `l500` = ℓ（pole-tip 至 WP 距離 / 電荷距離）變體。

**相關**：見 `scripts/README.md`、上層 `../README.md`。
