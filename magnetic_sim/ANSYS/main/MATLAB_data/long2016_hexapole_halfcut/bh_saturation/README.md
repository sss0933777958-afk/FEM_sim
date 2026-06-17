# MATLAB_data/long2016_hexapole_halfcut/bh_saturation — 鐵材 B-H 飽和掃描成果

**用途**：存放下極半切六極在 AISI 1018 鐵芯 **B-H 非線性（飽和）** 模型下，掃不同激發電流的場結果，用來看極尖 / WP 區何時進入飽和。

**內容**：
- `bh_saturation_sweep.mat`：B-H 非線性飽和電流掃描（不同 I 下的峰值 |B| 等），對照線性 μ_r 估飽和拐點。

**資料來源 / 流向**：由 `matlab/long2016_hexapole_halfcut/bh_saturation/` 的腳本讀 ANSYS_data 的 B-H 版 `.dat`（含 `_DC_BH` 系列解）產生。resolver = `matlab_path('long2016_hexapole_halfcut','bh_saturation')`。

**相關**：見 [../README.md](../README.md)（本 model 功能總覽）、[../../README.md](../../README.md)（MATLAB_data 總覽）與 [../../../CLAUDE.md](../../../CLAUDE.md)。
