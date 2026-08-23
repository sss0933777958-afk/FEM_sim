# doc/charge_model_fitting/Hall_snesor_base_fix_l/ — fix-ℓ 成本函數推導

**用途**：Hall-sensor 校正「**固定 ℓ̂、求每極 d**」(fix-ℓ) 流程的成本函數推導文件。
對應程式包：`matlab/long2016_hexapole_halfcut/Calibration_using_FEM_modeling/Hall_sensor_base_fix_dir/`
（模型 `b_ij = S_i·V_j·d`，no-gain；先 fit ℓ̂ 再閉式解對角 d）。

**內容**：
- `pdf/cost_function.pdf` — 成本函數 `J(d)` / `J(ℓ̂)` 推導。

**相關**：no-fix-ℓ（18-param bias）版見 `../Hall_sensor_base_no_fix_l/`；程式見上述 fix_dir 包
（`code/function/sensor_cost_lhat.m`、`solve_d.m`、`sensor_residual.m`）。

> 註：資料夾名 `snesor` 為既有拼字（應為 sensor）；改名屬結構變動，未經指示不擅動。
