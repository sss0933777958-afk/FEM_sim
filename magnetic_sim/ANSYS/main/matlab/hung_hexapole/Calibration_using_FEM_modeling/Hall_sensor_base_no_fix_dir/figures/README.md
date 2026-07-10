# …/Hall_sensor_base_no_fix_dir/figures/ — 圖檔

**用途**：本功能組的圖檔輸出（**直接放此層，不再多包子夾**）。

**內容**：sensor 處磁路箭頭圖（coil1 = P1 激發、all-source、**真實 FEM 節點不內插**、格點抽樣、箭頭=單位方向/色=|B|、cone 真實輪廓 + n+）：
- `P2sensor_Braw_P1exc.png` — P2 sensor（P1 的**對極**）；B 幾乎沿錐面、⊥ n+ → **近零耦合**（由 `plot_P2sensor_Braw_P1exc.m`，y=0 側視）。
- `P4sensor_Braw_P1exc.png` — P4 sensor（P1 的**鄰極上極**）；B 與 n+ 有明顯分量 → **強耦合 ~15%**（由 `plot_sensorBcircuit_P1exc.m(4)`，P4 axis–n+ 局部切面、world-up 顯示真傾角）。

**相關**：見上層 `../README.md`、`../code/plot/README.md`。
