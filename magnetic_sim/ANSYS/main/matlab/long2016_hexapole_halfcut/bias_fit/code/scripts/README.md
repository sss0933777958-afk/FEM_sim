# long2016_hexapole_halfcut/bias_fit/code/scripts/ — bias_fit 計算腳本

**用途**：no_fix_l 18-param bias 電荷模型的 fit / 掃描 / 比較 / 診斷腳本（不出圖）。

**內容**：
- `sweep_nofixl_vs_R.m` — 對 R=50:5:500µm 各做一次 bias fit，存 `ℓ̂(R)/gB(R)/K̄_I(R)/e_hat(R)` 趨勢。
- `calib_fem_bias.m` — bias 模型主校正骨架（actuator 框、profiled 殘差、lsqnonlin）。
- `compare_models_nrmse.m` — bias 模型 vs 舊全-K̂ 模型 NRMSE 比較。
- `check_charge_in_pole.m` — 檢查 fit 出的等效電荷是否落在極體內。

**資料來源 / 流向**：`ansys_path(model,'coilN','standard')` 讀 6 顆 coil FEM `.dat`（`'all'`）→ fit → `.mat` 寫 `matlab_path(model,'charge_fit',...)`；每半徑 `.tex` → `../../results/per_radius/`。模型電流 `I_actual=1`（對齊 FEM 1A）。

**命名 / 慣例**：純計算放此（`code/scripts`），繪圖放 `../plot/`；APDL→paper 索引 `[1,3,6,5,2,4]`。

**相關**：見 `../../README.md`、`../../../common/README.md`、`../../../../../CLAUDE.md`。
