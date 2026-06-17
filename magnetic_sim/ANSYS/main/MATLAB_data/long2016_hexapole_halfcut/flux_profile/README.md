# MATLAB_data/long2016_hexapole_halfcut/flux_profile — 沿極軸磁通剖面 Φ(s)

**用途**：存放沿磁極軸向座標 s 的**磁通 Φ(s) = ∮ B·dA** 積分剖面，用來看磁通在錐形極內的衰減 / 洩漏（P1 為 D-shape 半切盤、P2 為含 36.59° tilt frame 的整盤）。

**內容**：
- `P1_flux_profile.mat` / `P2_flux_profile.mat`：P1 / P2 沿極軸磁通 Φ(s)（現行 smrt5 批次）。
- `P1_flux_profile_smrt4.mat` / `P2_flux_profile_smrt4.mat`：舊 smrt4 批次對照（smrt4 halfcut 不可靠，僅留參考）。

**資料來源 / 流向**：由 `matlab/long2016_hexapole_halfcut/flux_profile/`（plot 組）讀 P1 / P2 自激的 FEM `.dat`（極內各截面節點 B），對各截面積分得 Φ(s)。resolver = `matlab_path('long2016_hexapole_halfcut','flux_profile')`。

**相關**：見 [../README.md](../README.md)（本 model 功能總覽）、[../../README.md](../../README.md)（MATLAB_data 總覽）與 [../../../CLAUDE.md](../../../CLAUDE.md)。
