# …/Hall_sensor_base_no_fix_dir/code/ — 求 d 一條龍程式（18 參數 bias 版）

**用途**：Hall-sensor per-pole `d`（18 參數 bias 版）的全部程式碼，分三層（比照 `../../Hall_sensor_base_fix_dir/code/`）：
- `main/` — 主程式 `main.m`（載 calib_bias.mat → actuator 框建 M,c → 求 d → LaTeX 的 driver）。
- `function/` — 模型數學 + IO 輔助函式；`main.m` 全部從這裡呼叫（一檔一函式）。
- `plot/` — d-form vs Q-form 的兩支驗證腳本（診斷用，印數字非出圖）。

**相關**：見上層 `../README.md`。
