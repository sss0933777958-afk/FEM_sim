# …/fix_l/code/plot/ — fix-ℓ 主程式的繪圖

**用途**：fix-ℓ 校正主程式對應的繪圖腳本（收斂、電荷/磁路視覺化等）。
**內容**：`plot_fixl_convergence.m`（fit 收斂）、`plot_P1_circuit_charge.m` / `plot_P2_circuit_charge.m`（單極磁路+**在軸**電荷，吃 `Rum`；本批 R=150、'zoom' 裁到尖端/WP 強場區）、`plot_pole_circuit_charge_3d.m`（3D）、`make_R040_charge_figs.m`（R=40 µm 電荷圖批次）。

**資料來源 / 流向**：讀 FEM 場（`ANSYS_data/.../coil1|coil5/standard`）+ `MATLAB_data/.../charge_fit/fit_KI_ball/fit_KI_R<R>.mat`（取 ℓ）→ 圖存 `../../figures/`。場圖一律畫真實 FEM 節點原值，不內插。all-source：P1（下極）`B=−B_FEM`、P2（上極）keep `B=+B_FEM`，兩極尖端皆射出。

**命名 / 慣例**：新增繪圖腳本前須依 `CLAUDE.md` 的繪圖腳本規則：先確認屬哪個功能組 → 一任務一腳本、原地反覆改到使用者定案 → 定案後才存最終圖（定案前用 MCP preview，不落地）；使用者沒說「新增」就不開第二支。

**相關**：見上層 `../README.md`、`../../../../../CLAUDE.md`（繪圖腳本規則）。
