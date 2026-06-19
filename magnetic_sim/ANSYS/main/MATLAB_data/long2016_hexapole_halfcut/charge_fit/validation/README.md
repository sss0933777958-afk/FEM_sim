# MATLAB_data/long2016_hexapole_halfcut/charge_fit/validation — 電流組合閉環驗證

**用途**：存放點電荷模型的**閉環驗證**結果（與同層 `calibration/`、`fit_KI_ball/`、`fitting_trend/`、`fitting_d/` 並列）。

**內容**：
- `validate_combos_R150.mat`：以 R*=150 校正參數對隨機電流組合做閉環檢查的 NRMSE（combos、normI、NRMSE、N、R_star）。

**資料來源 / 流向**：由 `validation/code/scripts/eval_validate_combos.m` 產生（讀 `calibration/calibration_final.mat`）；
對應圖在 `validation/figures/validate_combos_nrmse_R150.png`。目前無下游 loader。

**相關**：見 [../README.md](../README.md)（charge_fit 總覽）與 [../../../../CLAUDE.md](../../../../CLAUDE.md)。
