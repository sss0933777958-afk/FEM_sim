# …/no_fix_l/code/function/ — no-fix-ℓ 非 main 輔助函式

**用途**：非 main 輔助（legacy tex）。**main 校正流程函式（pipeline + 輸出 `emit_model_results`）已移到 `../main_function/`**（load_coils_actuator/select_ball/fit_bias/bias_resid/build_A/make_Pc/gauge_KI/region_field_err/calc_range_metrics/emit_model_results）。
**內容**：
- `write_KbarI_tex.m` — 輸出純結果 `.tex`。

**資料來源 / 流向**：`load_coils_actuator` 讀 `ANSYS_data/.dat` → fit/gauge/誤差 → `write_KbarI_tex` 寫 `../../results/no_fix_l/*.tex`。

**命名 / 慣例**：純函式（一檔一函式）；I=1 A 對齊 FEM；coil_sign 全 source 翻上極。

**相關**：見上層 `../README.md`。
