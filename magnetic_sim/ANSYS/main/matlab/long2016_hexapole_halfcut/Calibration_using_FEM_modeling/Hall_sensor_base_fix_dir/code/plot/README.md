# …/Hall_sensor_base_fix_dir/code/plot/ — 場視覺化 + sign 診斷

**用途**：P2 sensor / P2 磁極的場圖與 sign 診斷腳本。場圖一律畫**真實 FEM 節點原值**（不內插，
除非檔名標 interp）；圖存 `../../figures/`。

**內容**：
- `plot_P2sensor_Braw_P1exc.m` — P2-sensor 局部磁路箭頭圖（coil1 P1 激發、all-source）。模式
  `circuit`/`zoom`/`pole`/`why`/`nodes`；參數 `VARIANT`(`sensor_spheres`/`standard`)、`FLIP`(`flip`/`raw`)。
- `plot_P1P2_air_circuit_3d.m` — 3D 空氣磁路圖。`FOCUS` = `sensor`/`full`/`p2pole`；參數
  `VARIANT`/`DATASET`/`FLIP`（如 `'graded_p2','p2reg','none'` 畫 P2 整根極）。
- `plot_P2pole_circuit_2d.m` — **P2 整根磁路 2D 剖面**（WP→支撐座，graded 密網格、含鋼件輪廓）。
- `plot_P2sensor_tets_3d.m` — P2 sensor 取樣圓柱被真實 FEM tet 包住的 3D 圖。
- `plot_interp_tet_schematic.m` — 重心法內插示意（一顆四面體 + 內部點 + λ 權重，純示意）。
- `diag_P2_Bn_map.m` — 沿 P2 cone 兩側 flank × SOFF × gap 掃 `B·n+` sign，找匯進(負)/回流(正)區。
- `diag_P2P1_single.m` — 單 coil1（P1 激發）P2 sensor 單點內插 `B·n+` + sign；可帶 variant
  （`standard`/`gap*_mueq`/`mueq_s<i>`）比翻負臨界。
- `diag_Vmat_sign.m` / `diag_Vmat_sign_center.m` — 印 6×6 all-source Vmat sign 表（圓柱平均 / 底面中心
  單點內插兩法），看哪些 off-diagonal 為正。
- `compare_sensor_methods.m` — 完美接觸下三法 Vmat 比較（M1 單點內插＝參考 / M2 1000 點內插 / M3 加密
  real-node），算誤差 %（self/cross 分開）+ **Task4 每 sensor 取樣圓柱覆蓋幾個 tet**；輸出
  `results/sensor_methods_compare.txt` + `_data.mat`。
- `diag_sensor_node_spread.m` — 加密網格下，各 sensor 圓柱內 ~170 真實節點的 `B·n+` 散布（mean/std/CoV/
  全距 + 沿軸趨勢），看「節點值是否一致」。

**相關**：見上層 `../README.md`；資料來源 `ANSYS_data/.../data/`（讀前查 `RESULTS_MAP.md`）。
