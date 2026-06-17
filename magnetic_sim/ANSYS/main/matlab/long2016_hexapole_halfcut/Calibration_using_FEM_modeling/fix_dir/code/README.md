# …/fix_l/code/ — fix-ℓ 校正程式碼

**用途**：fix-ℓ 點電荷校正主程式的所有 MATLAB 碼，依角色分三組。
**內容**：
- `main/` — 主程式 `main.m`（驅動：load → select → fit → error → 出 `.tex`）。
- `function/` — 模型數學輔助函式（`load_coils`、`select_ball`、`fit_KI_fixl`、`charge_residual`、`region_field_err`、`unpack_params`、`write_KI_tex`）。
- `plot/` — 該主程式的繪圖腳本。

**資料來源 / 流向**：`main.m` 讀 `ANSYS_data/.dat` → 呼叫 `function/` 算 → `plot/` 畫圖、`function/write_KI_tex` 寫 `../results/fix_l/*.tex`。

**命名 / 慣例**：`code/{main,function,plot}/`；數學在 `function/`、driver 在 `main/`。

**相關**：見上層 `../README.md`。
