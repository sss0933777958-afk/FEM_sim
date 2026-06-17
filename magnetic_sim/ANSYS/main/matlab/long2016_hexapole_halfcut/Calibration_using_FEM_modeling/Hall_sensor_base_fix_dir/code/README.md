# …/sensor_d/code/ — 求 d 一條龍程式

**用途**：Hall-sensor per-pole `d` 的全部程式碼，分三層（比照 `../../fix_l/code/`）：
- `main/` — 主程式 `main.m`（求 d → 後處理 → 驗證印出 → LaTeX 的 driver；config 在頂部）。
- `function/` — 模型數學 + IO 輔助函式；`main.m` 全部從這裡呼叫（一檔一函式）。
- `plot/` — d-form vs Q-form 的兩支驗證腳本（診斷用，印數字非出圖）。

**相關**：見上層 `../README.md`。
