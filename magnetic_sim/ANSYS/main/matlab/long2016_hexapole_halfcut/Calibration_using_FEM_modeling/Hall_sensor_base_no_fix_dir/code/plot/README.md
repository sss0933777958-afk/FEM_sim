# …/Hall_sensor_base_no_fix_dir/code/plot/ — 驗證 + 磁路視覺化腳本

**用途**：畫 sensor 處磁路箭頭圖，理解 B 方向 vs n+（為何某些耦合近零）。

**內容**：
- **磁路箭頭圖**（coil1 = P1 激發、all-source、**真實 FEM 節點不內插**、格點抽樣降密度、箭頭單位方向/色=|B|）：
  - `plot_sensorBcircuit_P1exc.m`（**通用**，`(pole_i)` 任一極）— 在該極「pole_axis × n+」局部切面、**world-up 顯示**（看得出真極傾角）畫 sensor 處磁路；含真實 cone 輪廓 + n+。出 `../../figures/<Pn>sensor_Braw_P1exc.png`。本批用 `pole_i=4`（P4：P1 的鄰極上極，強耦合 ~15%）。
  - `plot_P2sensor_Braw_P1exc.m`（**P2 專用**，y=0 側視）— 同概念但 P2 剛好在 y=0，直接 x-z 框。出 `P2sensor_Braw_P1exc.png`（P2：P1 的對極，近零耦合）。
    > 通用版 `plot_sensorBcircuit_P1exc(2)` 也能畫 P2；目前兩支並存（暫定案保留）。

**相關**：見上層 `../README.md`、`../../figures/README.md`。
