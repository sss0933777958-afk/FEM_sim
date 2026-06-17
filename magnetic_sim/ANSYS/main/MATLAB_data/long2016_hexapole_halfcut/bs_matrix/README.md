# MATLAB_data/long2016_hexapole_halfcut/bs_matrix — B̄_S 感測矩陣與 V_out/V_in 轉移矩陣

**用途**：存放 6×6 **B̄_S（sensor 端磁通密度感測矩陣）** 與 **V_out/V_in（電壓轉移矩陣）** 的各版本成果，含 sensor 距極尖距離掃描與 gap / 放大器增益 k_A 變體。

**內容**（代表檔）：
- `Bbar_S_4p572.mat`（+ `_gap200um`）：sensor 在錐尖後 4.572 mm 的 6×6 B̄_S；`_gap200um` 為 200 µm 氣隙對照。
- `B_bar.mat` / `B_s_matrix.mat` / `B_s_final.mat` / `B_s_kA_final.mat` / `B_bar_s_corrected.mat` / `B_bar_S_toWP.mat` / `Bs.mat`：B_S 系列中間 / 修正 / 縮到 WP 版。
- `Vout_Vin_1p572.mat … Vout_Vin_4p572.mat`：V_out/V_in **sensor 距離掃描**（1.572–4.572 mm）；`_gap200um` / `_gap200um_mueq` 為氣隙 + μ_r 等效變體。
- `V_out_V_in_kA0p3.mat` / `Vout_Vin_kA0p36.mat`：不同放大器增益 k_A（0.3 / 0.36）版。

**資料來源 / 流向**：由 `matlab/long2016_hexapole_halfcut/bs_matrix/`（scripts）讀 6 顆 coil 的 1A FEM `.dat`，套 S_hall（130 V/T）+ I_in（0.6 A 操作點縮放）+ k_A 後產生。resolver = `matlab_path('long2016_hexapole_halfcut','bs_matrix')`。

**相關**：見 [../README.md](../README.md)（本 model 功能總覽）、[../../README.md](../../README.md)（MATLAB_data 總覽）與 [../../../CLAUDE.md](../../../CLAUDE.md)。
