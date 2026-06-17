# …/Hall_sensor_base_no_fix_dir/code/function/ — 求 d 模型/IO 輔助函式（18 參數 bias 版）

**用途**：18 參數 bias 版 Hall-sensor `d` 管線的數學/IO 輔助函式；`main.m` 全部從這裡呼叫。

**內容**：
- `build_A.m` — 堆疊空間函數矩陣 A(3Np×6)，actuator 框、電荷格點 Pc_18（移植自 `bias_fit/calib_fem_bias.m`）。
- `sensor_residual_bias.m` — sensor 模型 actuator 框相對 RMSE（`bm=g_H·A·(V_j.*d)` vs all-source FEM）。
- `build_sensor_geometry.m` — 每極 sensor 位置/法線 n+ + Ø0.3 mm disc 取樣格（frame-independent，與 fix-ℓ 版同）。
- `extract_Vmat.m` — 抽 sensor 電壓 `Vmat`（disc 面積平均 signed B·n+）+ all-source 翻下極欄（frame-independent）。
- `solve_d.m` — 閉式解 `d`（含增益 g_H；吃 bias 的 M、c）。
- `compute_KH.m` — derivation.pdf 流程：由 d 求 `ĝ_F`、`K̂_H`、`R_a`。
- `write_d_tex.m` / `write_KH_tex.m` — 輸出純結果 `.tex`。

**註**：`build_sensor_geometry / extract_Vmat / solve_d / compute_KH / write_*_tex` 與 `../../Hall_sensor_base_fix_dir/code/function/` 相同（frame-independent，直接複用）；bias 版差異只在 `build_A`（取代 build_S）與 `sensor_residual_bias`（actuator 堆疊版）。

**資料來源 / 流向**：`main.m` 載 `calib_bias.mat` + 讀 `ANSYS_data/.dat` → 這些函式建 M,c/求 d → `write_*_tex` 寫 `../../results/*.tex`。

**命名 / 慣例**：純函式（一檔一函式）；I=1 A 對齊 FEM；actuator 框、Pc_18 離軸；all-source 翻下極。

**相關**：見上層 `../README.md`。
