# MATLAB_data/long2016_hexapole_halfcut/freq_response — Fig 4.4 頻率響應 / cross-act

**用途**：存放下極半切六極的**頻率響應（Fig 4.4）** 成果：P1 self-act（自激自測）與 P5→P1 cross-act（P5 激發、P1 感測）的頻掃，含 COMSOL WP 平面對照。

**內容**（代表檔）：
- `fig44_P1.mat` / `fig44_P1_30freq.mat` / `fig44_P1_v2.mat`：P1 self-act 頻響（不同頻點數 / 版本）。
- `P1field_P5act.mat`（+ `_1kHz`）/ `fig44_P5act_P1sensor.mat` / `fig44_P5act_freqsweep.mat` / `wp_P5act.mat`：P5 激發、P1 感測 cross-act 與 WP 場。
- `I_in_30freq.mat`：30 頻點輸入電流。
- `comsol_wp_plane_P1.mat`：COMSOL WP 平面對照（驗證 FEM 頻響）。

**資料來源 / 流向**：由 `matlab/long2016_hexapole_halfcut/freq_response/`（scripts）讀 ANSYS AC / 頻掃 `.dat`（及 COMSOL LiveLink 抽出的 WP 平面）產生。resolver = `matlab_path('long2016_hexapole_halfcut','freq_response')`。

**相關**：見 [../README.md](../README.md)（本 model 功能總覽）、[../../README.md](../../README.md)（MATLAB_data 總覽）與 [../../../CLAUDE.md](../../../CLAUDE.md)。
