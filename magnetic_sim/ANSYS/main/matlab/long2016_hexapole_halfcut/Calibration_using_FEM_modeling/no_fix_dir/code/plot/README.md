# …/no_fix_l/code/plot/ — no-fix-ℓ 主程式的繪圖

**用途**：18-param bias 校正主程式對應的繪圖腳本。
**內容**：
- `plot_nofixl_convergence.m`（fit 收斂）。
- `plot_P1_circuit_charge.m` / `plot_P2_circuit_charge.m`（單極側視磁路 + **離軸** 等效電荷，吃 `Rum`；本批 R=150）。charge 由 `calib_bias.mat` 的 `ℓ̂·(R'·Pc_18(:,k))` 算（18-param bias，actuator→measure），有別於 fix_dir 的在軸 `ℓ·d̂`。P1 'zoom' 裁到尖端/WP 強場區。

**資料來源 / 流向**：讀 FEM 場（`ANSYS_data/.../coil1|coil5/standard`）+ `MATLAB_data/.../charge_fit/calib_bias.mat`（R、Pc_18、ell_hat）→ 圖存 `../../figures/`。場圖一律畫真實 FEM 節點原值，不內插。all-source：P1（下極）`B=−B_FEM`、P2（上極）keep `B=+B_FEM`，兩極尖端皆射出。

**命名 / 慣例**：新增繪圖腳本前須依 `CLAUDE.md` 的繪圖腳本規則：先確認屬哪個功能組 → 一任務一腳本、原地反覆改到使用者定案 → 定案後才存最終圖（定案前用 MCP preview，不落地）；使用者沒說「新增」就不開第二支。

**相關**：見上層 `../README.md`、`../../../../../CLAUDE.md`（繪圖腳本規則）。
