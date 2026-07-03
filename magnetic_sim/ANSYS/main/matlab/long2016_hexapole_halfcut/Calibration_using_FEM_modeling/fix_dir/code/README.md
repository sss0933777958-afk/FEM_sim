# …/fix_l/code/ — fix-ℓ 校正程式碼

**用途**：fix-ℓ 點電荷校正主程式的所有 MATLAB 碼，依角色分三組。
**內容**：
- `main/` — 主程式 `main.m`（驅動：load → select → fit → error → 出 `.tex`）。
- `main_function/` — **main.m 掛勾的 pipeline 函式**（`load_coils`、`select_ball`、`fit_KI_fixl`、`charge_residual`、`unpack_params`、`region_field_err`、`calc_range_metrics`）。
- `function/` — 非 main 的輔助（`emit_model_results`、`svd_SiHI_WP`、`write_KI_tex`）。
- `plot/` — 該主程式的繪圖腳本。

**資料來源 / 流向**：`main.m` 讀 `ANSYS_data/.dat` → 呼叫 `function/` 算 → `plot/` 畫圖、`function/write_KI_tex` 寫 `../results/fix_l/*.tex`。

**命名 / 慣例**：`code/{main,function,plot}/`；數學在 `function/`、driver 在 `main/`。

**相關**：見上層 `../README.md`。
