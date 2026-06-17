# MATLAB_data/long2016_hexapole_halfcut/charge_fit/fitting_trend — K_I 擬合的取樣半徑 R 掃描趨勢與選取

**用途**：存放點電荷 K_I 擬合隨**取樣半徑 R 掃描**的趨勢、NRMSE 視窗、重現性（CoV）與最佳 R* 選取的決策資料，用來定案校正半徑（結論 R*=150 µm）。

**內容**（代表檔）：
- `sweep_alln_vs_R.mat`：主檔，全節點擬合對 R 的掃描（ℓ̂、gB、‖K̂‖、NRMSE vs R）。
- `KI_trend_sweep.mat` / `_N1000.mat` / `_dense.mat`：不同抽樣數 / 密度的趨勢掃描。
- `nrmse_window_maxmin.mat` / `objective_R_select.mat` / `Rlo_decision.mat` / `Rhi_decision.mat`：NRMSE 視窗 max−min 與 R 選取目標函數 / 上下界決策。
- `combo_nrmse.mat` / `combo_basis_nrmse.mat` / `testset_error_50_500.mat`：電流組合 / 基底 / 測試集誤差。
- `fixl_convergence.mat` / `nofixl_convergence.mat` / `sweep_nofixl_vs_R.mat`：fix_l vs no_fix_l 兩種 gauge 的收斂與趨勢。
- `Rlo_cov.mat` / `repro_cov_vs_R.mat` / `conv_cov_vs_nredraws.mat`：重現性 CoV 對 R / 重抽次數。

**資料來源 / 流向**：由 `matlab/long2016_hexapole_halfcut/charge_fit/`（scripts，如 `sweep_KI_trend*`/`nrmse_window_maxmin`/`objective_R_select`）讀 6 顆 coil 的 1A FEM `.dat` 產生。resolver = `matlab_path('long2016_hexapole_halfcut','charge_fit')`，再進子夾 `fitting_trend/`。

**相關**：見 [../README.md](../README.md)（charge_fit 總覽）、[../../README.md](../../README.md)（本 model 功能總覽）與 [../../../../CLAUDE.md](../../../../CLAUDE.md)。
