# magnetic_sim/ANSYS/backup/hung/pdf/scripts — 報告 LaTeX 原稿（report LaTeX source）（歸檔/archived）

> **本設計已歸檔（非活躍），保留供參考。** Hung hexapole（build-workflow design）。目前唯一活躍設計為 `../../../main/`。

**用途**：擬合報告的 LaTeX 原始檔（`.tex`），編譯後產出上層 `../*.pdf`。

**內容**：
- `fitting.tex` — 擬合結果表（booktabs + siunitx + colortbl）。
- `fitting_l250.tex` / `fitting_l250_cube20.tex` — ℓ=250 µm（±50 µm 擬合區、cube20 變體）報告表。

**資料來源 / 流向**：表中數值來自 `../../analysis/fit/`（讀 `../../data/*.mat`）的擬合結果，手動填入 `.tex` → 編譯為 `../*.pdf`（pdf 為報告）。

**命名 / 慣例**：`l250` = ℓ=250 µm；`cube20` = ±20 µm 擬合 cube 變體。

**相關**：見 `../README.md`（PDF 清單）、`../../analysis/fit/README.md`、上層 `../../README.md`。
