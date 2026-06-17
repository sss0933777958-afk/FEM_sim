# long2016_hexapole_halfcut/bs_matrix/code/scripts/ — bs_matrix 計算腳本

**用途**：6×6 `B̄_S` / `B_S` / `V_out_V_in` 轉移矩陣計算與 LaTeX 產生。

**內容**（代表檔）：
- `compute_B_bar_matrix.m` — `B̄(i,j)=(130e-3/I)·B(i,j)`，sensor 盤面積平均（per `dh_derivation.pdf`）。
- `compute_B_s_final.m` — 最終 `B_S` 矩陣。
- `compute_V_out_V_in_kA0p3.m` / `compute_delay.m` — 含放大器增益的 V/V 矩陣與延遲。
- `gen_Vout_Vin_<dist>.m` — 各 sensor 距離變體（`1p572`…`4p572`，含 `gap200um_mueq`）的 V/V 產生器。
- `gen_Bbar_S_*.m` / `gen_Bs_latex.m` / `gen_Vout_Vin_latex.m` — LaTeX 產生器。

**資料來源 / 流向**：`ANSYS_data/<model>/coilN` 讀 6 顆 coil FEM `.dat`（FEM 0.6A = 操作電流，免縮放）→ 算 → `.mat` 寫 `MATLAB_data/<model>/bs_matrix/`；`.tex` → `../../results/` 與 `doc/Solve B_matrix/.../scripts/`。常數：`S_hall=130 V/T`、`k_A=0.3614 A/V`、`I_in=0.6 A`、APDL→paper `[1,3,6,5,2,4]`。

**命名 / 慣例**：`compute_*` 算矩陣、`gen_*` 出 LaTeX；`<dist>`（如 `4p572`）= sensor-to-tip 距離 mm（`p`=小數點）。讀 coilN 前依 `result-read-safety` 核指紋。

**相關**：見 `../../README.md`、`../../../common/README.md`、`../../../../../CLAUDE.md`。
