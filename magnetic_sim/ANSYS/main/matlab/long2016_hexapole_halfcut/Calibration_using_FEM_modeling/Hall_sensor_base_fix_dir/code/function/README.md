# …/Hall_sensor_base_fix_dir/code/function/ — 求 d 模型/IO 輔助函式

**用途**：Hall-sensor 每極 `d` 流程的數學/IO 輔助函式（一檔一函式）；`main.m` 全部從這裡呼叫。

**內容（5 支）**：
- `build_S.m` — 單點 3×6 點電荷空間核 `S(p̄; ℓ̂)`（庫倫式核，column k = (p̄−p̄_ck)/‖·‖³）。
- `build_sensor_geometry.m` — 6 顆 sensor 中心 `sensor_pos` + 法線 `sensor_n`（下極磨平面 / 上極錐面，n+ 出鋼）。
- `extract_Vmat.m` — 真實節點抽 sensor 電壓 `Vmat`：沿法線 n 開 0.15 mm 圓柱選真實節點對 B·n 平均（**不內插**，圓柱內無點則取沿 n 最近的盤內節點）+ all-source 翻下極欄。
- `solve_d.m` — 閉式解 `d`（minJ，含增益 g_H；重用 M、c）。
- `sensor_residual.m` — 算殘差 ε、成本 `J = Σ‖ε‖²`（只回 cost J，不算相對 RMSE；呼叫 `build_S`）。

**資料來源 / 流向**：`main.m` 讀 `ANSYS_data/.../data/coil1..6` 的 `.dat` → 這些函式算 M/c、Vmat、d、J。

**命名 / 慣例**：純函式；I=1 A 對齊 FEM；電荷全 source、all-source 翻下極。

**相關**：見上層 `../README.md`。
