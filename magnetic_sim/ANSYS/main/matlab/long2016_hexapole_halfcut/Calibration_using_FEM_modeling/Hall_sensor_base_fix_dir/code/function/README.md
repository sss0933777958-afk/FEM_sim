# …/Hall_sensor_base_fix_dir/code/function/ — 求 d 模型/IO 輔助函式

**用途**：Hall-sensor 每極 `d` 流程的數學/IO 輔助函式（一檔一函式）；`main.m` 全部從這裡呼叫。

**內容**：
- `build_S.m` — 單點 3×6 點電荷空間核 `S(p̄; ℓ̂)`（庫倫式核，column k = (p̄−p̄_ck)/‖·‖³）。
- `build_sensor_geometry.m` — 6 顆 sensor 中心 `sensor_pos` + 法線 `sensor_n`，**全域(WP 中心)座標**：
  - **下極(P1/P3/P6)**：磨平面，`sensor_pos = tip + 4.572mm·[cosθ;sinθ;0] + 0.41mm·[0;0;1]`、`n+=[0;0;1]`。
  - **上極(P2/P4/P5)**：自然錐面、軸用 CAD 實際傾角 `inc_up=cnst.upper_incline≈36.59°`（非理想魔術角）；沿錐面 `dir(inc_up+β,θ)`、法線 `dir(inc_up+β+90°,θ)`（`dir(e,a)=[cos(e)cos(a);cos(e)sin(a);sin(e)]`、β≈11.31°、n+ 指向 sensor 出鋼）。
- `extract_Vmat.m` — 真實節點抽 sensor 電壓 `Vmat`：沿法線 n 開 0.15 mm 圓柱選真實節點對 B·n 平均（**不內插**，圓柱內無點則取沿 n 最近的盤內節點）+ all-source 翻下極欄。
- `solve_d.m` — 閉式解 `d`（**no-gain，無 g_H**）：內外層雙重加總都在此檔，`d = (Σ_j V_j (Σ_i S_iᵀS_i) V_j)⁻¹ (Σ_j V_j (Σ_i S_iᵀ b_ij))`；模型 `b = S·V·d`。
- `sensor_residual.m` — 算殘差 ε、成本 `J = Σ‖ε‖²`（模型 `b = S·V·d`，無 g_H；只回 cost J；呼叫 `build_S`）。
- `sensor_cost_lhat.m` — profiled cost `J(ℓ̂)`（消去 d、per-j；給 `main.m` 一維 fminbnd 找 ℓ̂）。
- `extract_Vmat_interp.m` — **內插版** Vmat：粗網格在 sensor 圓柱內均勻撒點（預設 1000）做真·FEM tet
  重心內插（連接性取 `data/mesh/<variant>/csv/sensor_local_*.csv`），對 B·n+ 平均 + all-source 翻號。
  末參 `variant`（選填，預設 `'standard'`；第 11 參數）：場來源子夾，可填 `gap100um_mueq` 等同網格變體
  （網格拓樸仍用 standard CSV，只換場）。`calib_gap100um.m` 用 n_uniform=100 + variant='gap100um_mueq'。
- `extract_Vmat_interp_center.m` — **單點內插版** Vmat：在「圓柱底面中心」(= sensor 中心) 重心內插單點 B·n+
  （非平均）；可帶 `variant`（standard / gap*_mueq）比 sign。

**資料來源 / 流向**：`main.m` 讀 `ANSYS_data/.../data/coil1..6` 的 `.dat` → 這些函式算 Vmat、d、J（M/c 內外層加總都在 `solve_d` 內部）。

**命名 / 慣例**：純函式；I=1 A 對齊 FEM；電荷全 source、all-source 翻下極。

**相關**：見上層 `../README.md`。
