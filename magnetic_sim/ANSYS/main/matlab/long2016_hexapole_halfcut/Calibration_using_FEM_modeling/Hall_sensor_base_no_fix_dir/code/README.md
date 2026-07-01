# …/Hall_sensor_base_no_fix_dir/code/ — 求 d 一條龍程式（18 參數 bias 版）

**用途**：Hall-sensor per-pole `d`（18 參數 bias 版）的全部程式碼，分三層（比照 `../../Hall_sensor_base_fix_dir/code/`）：
- `main/` — 主程式 `main.m`（**唯一統一 driver**：18-param bias，`fit_bias` ℓ̂+ê → G(=D^v) → `extract_Vmat_interp`（下極 −β、Ø0.3×0.1 圓柱、10k 內插）抽 V → 解 Ĥ_V/D̄/^Bĝ_V（論文 notation）；輸出 **V、D̄、^Bĝ_V、Ĥ_V** 到 `../results/D_gap200um_mueq.pdf` + `../data/calib_D_gap200um_mueq.mat`（.mat field 名沿用 Dmat/g_V）。已取代舊 `main_Dmatrix`）。
- `function/` — 模型數學 + IO 輔助函式；`main.m` 全部從這裡呼叫（一檔一函式）。
- `plot/` — d-form vs Q-form 的兩支驗證腳本（診斷用，印數字非出圖）。

**相關**：見上層 `../README.md`。
