# …/Hall_sensor_base_fix_dir/code/ — 求 d 程式

**用途**：Hall-sensor 每極 `d` 的全部程式碼，分層：
- `main/` — 主程式：`main.m`（真實節點抽 V：fit ℓ̂ → extract_Vmat → minJ 解對角 d）；`main_interp.m`（**內插版**：standard 粗網格 tet 重心內插、**圓柱內均勻 1000 點**抽 V → 對角 d；存 `calib_sensor_d_interp.mat`）。
- `function/` — 模型數學 + IO 輔助函式（一檔一函式）。含 `extract_Vmat`（真實節點圓柱平均）與 `extract_Vmat_interp`（standard 粗網格 tet 重心內插，**圓柱內均勻撒 n_uniform 點，預設 1000**；用 `data/mesh/<variant>/csv/sensor_local_{nodes,elems}.csv` 連接性 + .dat 場）。
- `decouple/` — 滿 6×6 解耦 `D_H`：`solve_DH_full.m`（真實節點 Vmat）、`solve_DH_interp.m`（**內插版 1000 點 Vmat → D_H 6×6 + V_j 對角 6×6×6**；cost_J = fix_dir 自由電荷下界，存 `calib_DH_interp.mat`）。
- `plot/` — 場視覺化 + sign 診斷（9 支：plot_P2sensor_Braw_P1exc / plot_P1P2_air_circuit_3d / plot_P2pole_circuit_2d / plot_P2sensor_tets_3d / plot_interp_tet_schematic / diag_P2_Bn_map / diag_P2P1_single / diag_Vmat_sign(_center)）。**逐支用途見 `plot/README.md`**。圖存 `../figures/`。

**相關**：見上層 `../README.md`。
