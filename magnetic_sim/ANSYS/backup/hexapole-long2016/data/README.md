# magnetic_sim/ANSYS/backup/hexapole-long2016/data — 擬合結果資料（歸檔/archived）

> **本設計已歸檔（非活躍），保留供參考。** Long Fei 2016 dissertation 六極 hexapole（full）。現役設計請見 `../../../main/`。

**用途**：存放 `../analysis/` 各擬合腳本（charge model fitting）的輸出成果，供繪圖與報告引用，免重跑擬合。

**內容**：MATLAB `.mat`：
- `charge_model_fit.mat` — [A] 單極點電荷模型擬合（fit `ell` ~835 µm）。
- `joint_6coil_19param_fit.mat` — [J] 6-coil 聯合擬合（19 參數）。
- `all6_bias_fit.mat` — [B-6x] 全激發 + bias 擬合（**最終方法**，誤差 ~0.07%）。

**資料來源 / 流向**：`analysis(*.m) → data(本夾 .mat) → figures(.png) / 報告`。由 `../analysis/fit_*.m` 寫入，由繪圖腳本（`plot_fig26_B6x.m` 等）讀回。

**命名 / 慣例**：檔名對應擬合方法（[A]/[J]/[B-6x]，見 `../docs/fitting-methods.md`）；內容沿用 dissertation notation（`ell`, `K_I`, `R_a`, `delta_i` 等）。

**相關**：見 ../README.md、../docs/fitting-methods.md、../docs/charge-model-fitting.md。
