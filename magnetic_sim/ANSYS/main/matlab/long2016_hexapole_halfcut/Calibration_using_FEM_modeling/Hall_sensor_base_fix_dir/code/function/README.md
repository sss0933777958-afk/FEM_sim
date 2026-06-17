# …/sensor_d/code/function/ — 求 d 模型/IO 輔助函式

**用途**：Hall-sensor per-pole `d` 管線的數學/IO 輔助函式；`main.m` 全部從這裡呼叫。多數逐字移植自 `../../../fixl_fit/code/scripts/calib_fem.m`（PAGE 1/2 區域函式）與 sensor_d 舊 `gen_*_latex.m`。

**內容**：
- `build_S.m` — 單點 3×6 空間函數矩陣 `S(pbar; ℓ)`。
- `cost_J.m` — 剖面化最小二乘成本 `J(ℓ)`（擬合 ℓ̂ 用；呼叫 `build_S`）。
- `build_sensor_geometry.m` — 每極 sensor 位置/法線 n+ + Ø0.3 mm disc 取樣格（721 點）。
- `extract_Vmat.m` — 抽 sensor 電壓 `Vmat`（disc 面積平均 signed B·n+）+ all-source 翻下極欄。
- `solve_d.m` — 閉式解 `d`（含增益 g_H；重用 page-1 的 M、c）。
- `sensor_residual.m` — sensor 模型相對 RMSE（呼叫 `build_S`）。
- `compute_KH.m` — derivation.pdf 流程：由 d 求 `ĝ_F`、`K̂_H`、`R_a`。
- `write_d_tex.m` / `write_KH_tex.m` — 輸出純結果 `.tex`。

**資料來源 / 流向**：`main.m` 讀 `ANSYS_data/.dat` → 這些函式算 ℓ̂/d/K_H → `write_*_tex` 寫 `../../results/sensor_d/*.tex`。

**命名 / 慣例**：純函式（一檔一函式）；I=1 A 對齊 FEM；電荷全 source、all-source 翻下極。

**相關**：見上層 `../README.md`。
