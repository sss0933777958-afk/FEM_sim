# …/fix_l/code/main_function/ — main.m 掛勾的 pipeline 函式

**用途**：`main.m` 的校正 pipeline 依賴函式（由 `../function/` 移來，2026-07-03 重構）。
**內容**：`load_coils`（載 6-coil FEM，旋 actuator 框）、`select_ball`（取 R 內節點）、`fit_KI_fixl`（lsqnonlin fit {K̄,ℓ̂,ĝ_I}）、`charge_residual`（殘差）、`unpack_params`（打包/解包）、`region_field_err`（區域 RMS 誤差）、`calc_range_metrics`（R≤150µm 球 σ_tot/iso_tot）。

**慣例**：一檔一函式；main.m 與 `../plot/` 的誤差圖都 addpath 本夾 + `../function/`。非 main 輔助（emit/svd/tex）留在 `../function/`。

**相關**：見 `../README.md`。
