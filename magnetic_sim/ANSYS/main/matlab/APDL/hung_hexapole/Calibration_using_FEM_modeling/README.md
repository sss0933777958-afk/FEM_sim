# hung_hexapole/Calibration_using_FEM_modeling/ — current_base / voltage_base

**2026-07-15**：與 long2016 一致的 2 專案結構（原 4 子夾合成 `current_base` + `voltage_base`），
`USE_BIAS` + variable projection（Varpar）統一 fix/no_fix；PDF 產生集中在 `code/main_function/emit_*.m`。
- **`current_base/`** — 電荷/電流模型（`load_coils_actuator/select_ball/fit_varpro/make_Pc/gauge_KI`）。
- **`voltage_base/`** — Hall-sensor/電壓模型（依賴 current_base + `build_sensor_geometry/extract_Vmat_interp`）。
- 讀 `ANSYS_data/hung_hexapole/data/<variant>/coilN/`（`gap_200um`/`no_gap`）。

> ⚠ **hung pipeline 尚未 runtime 驗證**：hung 上/下極 FEM 網格節點數不同（401632 vs 446944），
> `load_coils_actuator` 的「6 coil 共用同一 air-node 格」假設在 hung 會於 coil2 size-mismatch
> （**既有限制、非本次重構造成**）。要跑 hung 需先讓 loader 支援逐 coil 各自 air 節點。
> 結構已比照 long2016 對齊；程式碼由 long2016 canonical 版帶入（model/path 已改 hung）。
