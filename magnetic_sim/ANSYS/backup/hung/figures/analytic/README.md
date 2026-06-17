# magnetic_sim/ANSYS/backup/hung/figures/analytic — 解析 / 分析圖（analysis figures）（歸檔/archived）

> **本設計已歸檔（非活躍），保留供參考。** Hung hexapole（build-workflow design）。目前唯一活躍設計為 `../../../main/`。

**用途**：非 coil-specific 的分析圖：電荷模型擬合誤差、方向向量 d̂、Hall sensor 靈敏度、幾何標註、mesh 收斂。

**內容**（代表檔）：
- 擬合：`charge_model_cost_landscape.png`、`charge_positions_J.png`、`fitting_J_idealKI_{quiver_error,RMSE}.png`、`fitting_B6x_quiver_error.png`。
- 方向向量：`dhat_vs_pole_axis.png`、`dhat_3D_Long.png`。
- Hall sensor：`sensor_sensitivity_*.png`、`sensor_ratio_*.png`、`sensor_Btop_*.png`。
- 幾何 / mesh：`hexapole_tip_annotated.png`、`mesh_steel_filleted_3D.png`、`WP_convergence_Hung.png`。

**資料來源 / 流向**：由 `../../analysis/plot/`、`../../analysis/fit/` 讀 `../../data/*.mat` / `../../results/.../*.dat` 產出（pipeline：apdl→results(.dat)→analysis→data(.mat)→figures）。

**命名 / 慣例**：coil-specific 場圖在同層 `../coil1`–`../coil6`；本層只放 analytic（分析 / sweep / 標註）圖。

**相關**：見 `../README.md`（逐檔說明與圖風格）、上層 `../../README.md`。
