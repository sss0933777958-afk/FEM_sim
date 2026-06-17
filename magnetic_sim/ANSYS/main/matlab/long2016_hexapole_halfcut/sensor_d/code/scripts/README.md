# long2016_hexapole_halfcut/sensor_d/code/scripts/ — sensor_d 計算/驗證腳本

**用途**：Hall-sensor per-pole 常數 `d` 的求解、驗證與 LaTeX 產生。

**內容**：
- `gen_d_latex.m` — 出 `d_v2.tex`（無增益，d~1e-3）+ `d_final.tex`（含 `g_H`，d~1e-8），`d_v2 = d_final × g_H`。
- `gen_KH_latex.m` — 出 `K̂_H` LaTeX。
- `verify_d_vs_kmQS.m` — 驗 `d` 對 quasi-static km。
- `verify_dform_Qform_nrmse.m` — 比 d-form vs Q-form 殘差 NRMSE。

**資料來源 / 流向**：讀 `MATLAB_data/long2016_hexapole_halfcut/charge_fit/calib_sensor_d.mat` → 算/排欄（激發 coil → paper pole）→ LaTeX 寫 `doc/charge_model_fitting/long2016_hexapole_halfcut/scripts/`。

**命名 / 慣例**：`gen_*`=出 LaTeX、`verify_*`=驗證；`d_v2`=Calibration_V2 無增益版、`d_final`=hw.pdf 含增益版。

**相關**：見 `../../README.md`、`../../../common/README.md`、`../../../../../CLAUDE.md`。
