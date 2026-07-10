# …/no_fix_l/code/main_function/ — main 校正流程函式（pipeline + 輸出）

**用途**：`main.m` 的 18-param bias 校正流程函式（由 `../function/` 移來）：pipeline（main 直接呼叫）+ 輸出階段 `emit_model_results`（獨立跑、讀 main 存的 `.mat` 出 PDF）。
**內容**：`load_coils_actuator`（載 6-coil FEM、旋 actuator）、`select_ball`、`fit_bias`（lsqnonlin fit {ℓ̂,ê17}）、`bias_resid`（殘差）、`make_Pc`（ê+Pc_base→電荷）、`build_A`（設計矩陣）、`gauge_KI`（profile g_j → ^Bĝ_I/K̄）、`region_field_err`、`calc_range_metrics`（σ_tot/iso_tot）、**`emit_model_results`**（輸出階段：讀 fit `.mat` 出結果 PDF 到 `../results/`）、**`fit_singlepole_bias`**（單極有-bias fit `{ℓ̂,e_y,e_z,ĝ_I}`：電荷離軸橫向位移 `pc=ℓ·[1,e_y,e_z]`、lsqnonlin；被 `../main/main_singlepole_bias.m` 呼叫；單極路徑 reuse fix_dir 的 `load_singlepole`）。

**慣例**：一檔一函式。**Hall_sensor_base_{fix,no_fix}_dir 的 main/plot 也 addpath 本夾**（借用 load_coils_actuator/select_ball/build_A/make_Pc/fit_bias）。非 main 輔助（tex）在 `../function/`。

**相關**：見 `../README.md`。
