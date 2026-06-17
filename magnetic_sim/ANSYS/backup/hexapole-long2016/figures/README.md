# magnetic_sim/ANSYS/backup/hexapole-long2016/figures — 發表圖（歸檔/archived）

> **本設計已歸檔（非活躍），保留供參考。** Long Fei 2016 dissertation 六極 hexapole（full）。現役設計請見 `../../../main/`。

**用途**：本設計的發表／報告圖（publication figures）輸出根目錄。

**內容**：目前僅一個子夾 `long/`（見 `long/README.md`），收 dissertation 圖重製（fig2_3/2_4/2_6、fig26 系列）、P2 cone 場分析圖、charge-fit 驗證圖、PPT slide 圖。`.png` 圖檔本身放在 `long/`。

**資料來源 / 流向**：`analysis(*.m 繪圖腳本) → figures(本夾 .png)`。圖由 `../analysis/generate_*.m` / `plot_*.m` 讀 `../results/*.dat` 與 `../data/*.mat` 產出。

**命名 / 慣例**：場圖一律畫真實 FEM 節點原值（不用 scatteredInterpolant / 格點內插，除非明確標示內插）；WP 區用 mT；論文極名 P1–P6。

**相關**：見 ../README.md、figures/long/README.md、../analysis/。
