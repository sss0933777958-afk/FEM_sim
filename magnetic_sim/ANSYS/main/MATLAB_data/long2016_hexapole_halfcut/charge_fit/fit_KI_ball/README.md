# MATLAB_data/long2016_hexapole_halfcut/charge_fit/fit_KI_ball — 各取樣半徑的 K_I 球擬合

**用途**：存放在**不同取樣球半徑 R**（以 WP 為心、半徑 R 的球面 / 球內節點）上做點電荷 K_I 擬合的結果與收斂過程，供 `fitting_trend/` 選定最佳 R*。

**內容**（代表檔）：
- `fit_KI_R040.mat` / `fit_KI_R050.mat` / `fit_KI_R100.mat` / `fit_KI_R150.mat` / `fit_KI_R500.mat`：各半徑（40/50/100/150/500 µm）單點 K_I 擬合結果。
- `fit_KI_ball_sweep.mat`：跨半徑整批掃描彙整（個別 R200–R450 已刪，可由此重生）。
- `KI_convergence_gB50.mat`：固定 gB 的收斂歷程；`KI_perpoint_gB50_R040.mat`：R040 逐點殘差。

**資料來源 / 流向**：由 `matlab/long2016_hexapole_halfcut/charge_fit/`（scripts，如 `sweep_KI_radius*`/`plot_KI_convergence`）讀 6 顆 coil 的 1A FEM `.dat` 產生。resolver = `matlab_path('long2016_hexapole_halfcut','charge_fit')`，再進子夾 `fit_KI_ball/`。

**相關**：見 [../README.md](../README.md)（charge_fit 總覽）、[../../README.md](../../README.md)（本 model 功能總覽）與 [../../../../CLAUDE.md](../../../../CLAUDE.md)。
