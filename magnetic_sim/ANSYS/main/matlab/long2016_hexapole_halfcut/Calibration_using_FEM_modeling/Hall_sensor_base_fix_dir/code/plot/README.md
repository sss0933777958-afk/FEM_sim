# …/sensor_d/code/plot/ — d 模型驗證腳本

**用途**：驗證求出的 `d` 自洽（診斷用，印數字、非出圖）。讀 `calib_sensor_d.mat`（由 `../main/main.m` 產出）。

**內容**：
- `verify_d_vs_kmQS.m` — 驗 d-form `b=g_H S diag(d) V` 與 Q-form `b=(k_m/ℓ̂²) S Q̂`（`Q̂=(1/μ0)diag(d)V`）逐元素相同（機器精度）；含橋接常數 `g_H` vs `k_m/(ℓ̂²μ0)`。
- `verify_dform_Qform_nrmse.m` — 驗兩 form 在 R≤50 µm 的 NRMSE 完全相等（機器精度）。

**注意**：放在 `plot/` 是比照 fix_l 三層架構；這兩支實為**驗證/診斷**腳本（無圖）。

**相關**：見上層 `../README.md`。
