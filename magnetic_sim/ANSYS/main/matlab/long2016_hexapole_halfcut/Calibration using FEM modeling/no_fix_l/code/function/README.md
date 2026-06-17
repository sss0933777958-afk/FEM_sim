# …/no_fix_l/code/function/ — no-fix-ℓ 模型輔助函式

**用途**：18-param bias 點電荷模型的數學/IO 輔助函式；`main.m` 全部從這裡呼叫。
**內容**：
- `load_coils_actuator.m` — 載入 6-coil FEM 場並旋到 actuator frame（讀 `.dat`）。
- `select_ball.m` — 取 WP 半徑 R 內節點。
- `fit_bias.m` — `lsqnonlin` fit `{ℓ, ê(17)}`（呼叫 `bias_resid`、`make_Pc`、`build_A`）。
- `bias_resid.m` — 殘差；`make_Pc.m` — 由 ê + Pc_base 組電荷座標；`build_A.m` — 建設計矩陣。
- `gauge_KI.m` — profile g_j 並 gauge 出 gB、K̄_I。
- `region_field_err.m` — 區域相對 RMS 場誤差；`write_KbarI_tex.m` — 輸出純結果 `.tex`。

**資料來源 / 流向**：`load_coils_actuator` 讀 `ANSYS_data/.dat` → fit/gauge/誤差 → `write_KbarI_tex` 寫 `../../results/no_fix_l/*.tex`。

**命名 / 慣例**：純函式（一檔一函式）；I=1 A 對齊 FEM；coil_sign 全 source 翻上極。

**相關**：見上層 `../README.md`。
