# …/Hall_sensor_base_fix_dir/code/decouple/ — 滿 6×6 解耦矩陣 D_H（對照對角 d）

**用途**：獨立對照腳本，解「滿矩陣」解耦校正 `D_H`（6×6），把交叉耦合的 sensor 電壓反解成每極電荷：
`q_j = D_H·Vmat(:,j)`（對照 `code/main/main.m` 的對角版 `q_j = diag(d)·Vmat(:,j)`）。模型場 `b = S·D_H·V`（無 g_H、all-source）。

**內容**：
- `solve_DH_full.m` — 載 ℓ̂（fix_dir fit_KI_fixl，R=150）→ 載 6-coil FEM/取 R≤150 球 → `extract_Vmat` → 閉式解
  `D_H = (M\C)/Vmat`（`M=ΣSᵀS`、`C=ΣSᵀb`）→ 算 cost_J、印 D_H 與非對角量級 → 存 `calib_DH_full.mat`。

**關鍵結果（R=150 µm）**：
- `cost_J(滿 D_H) = 4.2611e-03 T²` **= fix_dir 自由電荷 cost（完全對上）**；對角 d 則 0.1419（33×）。
- `D_H` 非對角大（`‖off‖/‖diag‖≈0.83`）→ sensor 耦合強，這就是對角 d 殘差大的原因。
- 原理：`M\C` = 自由電荷最佳解（= fix_dir 找的），`D_H=(自由電荷)·Vmat⁻¹` = 物理解耦器 `G⁻¹`；
  `q_j=D_H·Vmat(:,j)` 還原自由電荷 → cost 落在點電荷模型下界。**非對角 = cross-talk 反解權重**。

**重用**（不複製）：`../function/{build_S, build_sensor_geometry, extract_Vmat}`；ℓ̂ 來源 `../../../fix_dir/data/fit_fixl_R150um.mat`（規則#2）。
**輸出**：`../../data/calib_DH_full.mat`（本組 `data/`，規則#2；`D_H`/`cost_J`/`ell_hat`/`Vmat`/`exc_sign`）。

**與 main 的關係**：對照用，**不動** main（對角 d = LAB406 per-pole 模型）。差別純粹是「校正矩陣對角 vs 滿」。

**相關**：見上層 `../../README.md`、`../function/README.md`。
