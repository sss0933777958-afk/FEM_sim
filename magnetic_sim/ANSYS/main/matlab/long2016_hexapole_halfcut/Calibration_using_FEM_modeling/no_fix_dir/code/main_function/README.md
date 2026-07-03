# …/no_fix_l/code/main_function/ — main.m 掛勾的 pipeline 函式

**用途**：`main.m` 的 18-param bias 校正 pipeline 依賴函式（由 `../function/` 移來，2026-07-03 重構）。
**內容**：`load_coils_actuator`（載 6-coil FEM、旋 actuator）、`select_ball`、`fit_bias`（lsqnonlin fit {ℓ̂,ê17}）、`bias_resid`（殘差）、`make_Pc`（ê+Pc_base→電荷）、`build_A`（設計矩陣）、`gauge_KI`（profile g_j → ^Bĝ_I/K̄）、`region_field_err`、`calc_range_metrics`（σ_tot/iso_tot）。

**慣例**：一檔一函式。**Hall_sensor_base_{fix,no_fix}_dir 的 main/plot 也 addpath 本夾**（借用 load_coils_actuator/select_ball/build_A/make_Pc/fit_bias）。非 main 輔助（emit/tex）在 `../function/`。

**相關**：見 `../README.md`。
