# …/no_fix_l/code/plot/ — no-fix-ℓ 主程式的繪圖

**用途**：18-param bias 校正主程式對應的繪圖腳本。
**內容**：`plot_nofixl_convergence.m`（fit 收斂）。

**資料來源 / 流向**：讀 FEM 場（`ANSYS_data/.dat`）/ fit 結果 → 圖存 `../../figures/`。場圖一律畫真實 FEM 節點原值，不內插。

**命名 / 慣例**：新增繪圖腳本前須依 `CLAUDE.md` 的繪圖腳本規則：先確認屬哪個功能組 → 一任務一腳本、原地反覆改到使用者定案 → 定案後才存最終圖（定案前用 MCP preview，不落地）；使用者沒說「新增」就不開第二支。

**相關**：見上層 `../README.md`、`../../../../../CLAUDE.md`（繪圖腳本規則）。
