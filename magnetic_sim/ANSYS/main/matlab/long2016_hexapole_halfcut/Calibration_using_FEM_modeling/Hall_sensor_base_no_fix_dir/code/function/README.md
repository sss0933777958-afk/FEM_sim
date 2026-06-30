# …/Hall_sensor_base_no_fix_dir/code/function/ — 求 d 模型/IO 輔助函式（18 參數 bias 版）

**用途**：18 參數 bias 版 Hall-sensor `d` 管線的數學/IO 輔助函式；`main.m` 全部從這裡呼叫。

**內容**：
- `build_S.m` — 單點空間函數矩陣 `S(p̄;ℓ̂)`（3×6，`(p̄−Pc_k)/‖·‖³`，p̄=p/ℓ̂）；與 fix-ℓ 版**完全相同**（取代舊堆疊版 `build_A`）。
- `sensor_residual_bias.m` — sensor 模型 **cost J = Σ‖ε‖²**（`bm=S_i V_j d`，**build_S 逐點**，no-gain；無 RMSE）。
- `build_sensor_geometry.m` — 每極 sensor 位置/法線 n+（只回 `sensor_pos`/`sensor_n`；frame-independent，與 fix-ℓ 版相同）。
- `extract_Vmat.m` — 抽 sensor 電壓 `Vmat`（**沿 n+ 圓柱選真實 FEM 節點**平均 signed B·n+、**非內插**）+ all-source 翻下極欄（含 `variant` 讀 mesh 變體子夾；與 fix-ℓ 版相同）。
- `solve_d.m` — 閉式解 `d`（**no-gain，無 g_H**；吃 bias 的 M、c）。

**註**：`build_S / build_sensor_geometry / extract_Vmat` 與 `../../Hall_sensor_base_fix_dir/code/function/` **完全相同**；`sensor_residual_bias` = fix-ℓ 版 `sensor_residual` 的 build_S 逐點 cost（只差檔名），`solve_d` 為 no-gain（吃 `main` 以 build_S 建好的 M、c，與 fix-ℓ 版數學等價）。**bias 唯一差異已收斂到 `main.m`**（載 calib_bias、旋進 actuator 框、用 Pc_18）。已移除 `build_A / compute_KH / write_*_tex`（堆疊 / g_H / K_H / LaTeX）。

**資料來源 / 流向**：`main.m` 載 `calib_bias.mat`（legacy MATLAB_data）+ 讀 `ANSYS_data/.dat` → 這些函式建 M,c / 求 d → 存 `../../data/calib_sensor_d_no_fix_dir.mat`（規則#2；`d`、`Vmat`、cost `J`）。

**命名 / 慣例**：純函式（一檔一函式）；I=1 A 對齊 FEM；actuator 框、Pc_18 離軸；all-source 翻下極。

**相關**：見上層 `../README.md`。
