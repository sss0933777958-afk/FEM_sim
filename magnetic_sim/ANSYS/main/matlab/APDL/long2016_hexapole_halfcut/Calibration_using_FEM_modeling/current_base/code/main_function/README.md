# …/current_base/code/main_function/ — main 校正流程函式（唯一函式夾）

**用途**：`main.m` 的六極校正流程函式，一步一支、命名清楚；`main.m` 只當 driver 逐個呼叫，末端也在此產 PDF。
（2026-07-16 重整：合併原 `function/` + `main_function/` 成單一 `main_function/`，並改名簡化。）

**六極 pipeline（main.m 依序呼叫）**：
- `load_coils_actuator` — 載 6-coil FEM、剔鐵、all-source、旋轉 actuator frame。
- `select_ball` — 取 R 球內真實節點 → P、Bstack。
- `build_S_matrix` — 每點 S_matrix kernel（原 `build_A`；fitting 與 solve 共用）。
- `fitting` — 變數投影疊代解 ℓ̄（+ê 若 USE_BIAS），回傳 `(ℓ̄, ê, Pc, J)`。
  **合併了原 `fit_varpro` + `bias_resid` + `make_Pc`**（殘差與 e→Pc 組裝都是其內部 nested function）。
- `solve_KI_bar_gain` — 由擬合解 K̄_I、G、^Bĝ_I（gauge K̄(1,1)=5/6；原 `gauge_KI`）。
- `calc_ellipsoid` — SVD 控制範圍性能量 𝒞/κ（原 `calc_range_metrics`）。
- `emit_model_results` — **main 最後一步呼叫**：讀剛存的 `.mat` 出結果 PDF（K̄_I/ℓ̂/G/F/^Bĝ_I，18-param 加印 ê）到 `../results/{single_param|eighteen_param}/`。
- `emit_tex` — 共用 LaTeX 排版 helper（`T.mat`/`T.scalar_unit`/`T.e`），供 `emit_model_results` 呼叫。

**單極次要 pipeline（`../main/main_singlepole*.m` 用；不屬六極統一流程）**：
`load_singlepole`、`fit_singlepole`、`fit_singlepole_bias`（`sweep_singlepole_tipcut` 無 main 掛鉤 → 移至 `../function/`）。

**已刪除**：`fit_bias`（被 `fitting(...,true)` 取代）、`emit_charge_Q`。**已合併/改名** 見上。

**慣例**：一檔一函式。voltage_base 的 main/plot 也 addpath 本夾（借用 load_coils_actuator/select_ball/build_S_matrix/fitting）。

**相關**：見 `../README.md`。
