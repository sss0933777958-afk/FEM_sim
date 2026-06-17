# magnetic_sim/ANSYS/backup/hexapole-long2016/analysis — MATLAB 後處理／擬合／繪圖（歸檔/archived）

> **本設計已歸檔（非活躍），保留供參考。** Long Fei 2016 dissertation 六極 hexapole（full）。現役設計請見 `../../../main/`。

**用途**：讀 FEM 輸出（`../results/coilN/*.dat`）做超疊加驗證、等效電荷模型擬合（charge model fitting）、論文圖（dissertation figures）重製。

**內容**：MATLAB `.m` 腳本：
- 共用工具（utility）：`import_ansys_data.m`（讀 `.dat`）、`mt_constants.m`（幾何／物理常數）、`point_charge_model.m`、`load_prpath.m`、`filter_iron_nodes.m`。
- 擬合（fitting）：`fit_charge_model.m`（[A] 單極點電荷）、`test_joint_6coil_fit.m` / `test_3d_charge_fit.m`（[J] 6-coil 聯合）、`fit_all6_with_bias.m`（[B-6x] 全激發 + 18 bias 參數，**最終方法**）。
- 驗證（verification）：`verify_superposition.m`、`verify_B_perpendicularity.m`、`analyze_P2_cone_Bfield.m`。
- 繪圖（figures）：`generate_figures_2_3.m` / `_2_4.m` / `_2_6.m`、`generate_fig26_coil1.m`、`plot_fig26_B6x.m`。

**資料來源 / 流向**：`apdl → results(.dat) → analysis(本夾) → data(.mat) → figures(.png)`。本夾讀 `../results/coilN/*.dat`，寫擬合結果到 `../data/*.mat`，圖存到 `../figures/long/`。

**命名 / 慣例**：擬合方法以 dissertation notation 命名（[A]→[J]→[B-6x]，見 `../docs/fitting-methods.md`）；極性以論文極名 P1–P6（非 APDL coil 索引）；符號慣例見 `../docs/notation-glossary.md`、`../docs/coil-winding-sign-convention.md`。

**相關**：見 ../README.md、../docs/（fitting-methods、notation-glossary、charge-model-fitting）。
