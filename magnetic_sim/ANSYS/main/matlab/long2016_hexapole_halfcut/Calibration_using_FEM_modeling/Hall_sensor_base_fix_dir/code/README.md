# …/Hall_sensor_base_fix_dir/code/ — 求 d 程式

**用途**：Hall-sensor 每極 `d` 的全部程式碼，分兩層：
- `main/` — 主程式 `main.m`（指定流程 driver：載 ℓ̂ → 真實節點抽 V → 建模型/殘差 → minJ 解 d；config 在頂部）。
- `function/` — 模型數學 + IO 輔助函式（一檔一函式），`main.m` 全部從這裡呼叫。

**相關**：見上層 `../README.md`。
