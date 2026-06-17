# …/no_fix_l/code/ — no-fix-ℓ（18-param bias）校正程式碼

**用途**：18-param bias 點電荷校正主程式的所有 MATLAB 碼，依角色分三組。
**內容**：
- `main/` — 主程式 `main.m`（驅動：load → select → fit_bias → gauge → error → 出 `.tex`）。
- `function/` — 模型數學輔助函式（`load_coils_actuator`、`select_ball`、`fit_bias`、`bias_resid`、`build_A`、`make_Pc`、`gauge_KI`、`region_field_err`、`write_KbarI_tex`）。
- `plot/` — 該主程式繪圖。

**資料來源 / 流向**：`main.m` 讀 `ANSYS_data/.dat`（經 `ansys_path`）→ `function/` 算 → `plot/` 畫圖、`write_KbarI_tex` 寫 `../results/no_fix_l/*.tex`。

**命名 / 慣例**：`code/{main,function,plot}/`；數學在 `function/`、driver 在 `main/`。

**相關**：見上層 `../README.md`。
