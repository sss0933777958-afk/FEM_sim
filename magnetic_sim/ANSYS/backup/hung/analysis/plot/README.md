# magnetic_sim/ANSYS/backup/hung/analysis/plot — 繪圖腳本（figure generation）（歸檔/archived）

> **本設計已歸檔（非活躍），保留供參考。** Hung hexapole（build-workflow design）。目前唯一活躍設計為 `../../../main/`。

**用途**：產生 B-field 分布、擬合誤差（quiver / RMSE）、電荷投影、|B| contour 等圖。

**內容**（代表檔）：
- `plot_Bfield_2d.m` / `plot_Bfield_3d.m` — Long2016-style 向量 + contour（2D / 3D 鐵芯磁通）。
- `plot_Bcontour_xaza*.m` — x_a–z_a 平面 |B| contour（baseline / l_sensitivity 版本）。
- `plot_Bvector_P1_*.m` — P1 軸與投影向量圖。
- `plot_J_quiver.m` / `plot_J_rmse.m` / `plot_charge_proj.m` / `plot_J_*.m` — [J] FEM vs 模型誤差與電荷投影。

**資料來源 / 流向**：讀 `../../results/coilN/<變體>/*.dat`（場圖）或 `../../data/*.mat`（擬合圖）→ 輸出 `.png` 至 `../../figures/coilN/` 或 `../../figures/analytic/`。

**命名 / 慣例**：場圖一律畫真實 FEM 節點原值（除非明確要求內插並標示）。`coilN` = 各極自激；`analytic/` = 非 coil-specific 分析圖。

**相關**：見 `../README.md`、`../../figures/README.md`、上層 `../../README.md`。
