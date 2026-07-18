# NTU_hexapole/Calibration_using_FEM_modeling/ — current_base / voltage_base（骨架）

**2026-07-15**：與 long2016 一致的 2 專案結構（4 子夾合成 `current_base` + `voltage_base`）。
NTU 目前為**空骨架**（尚無校正程式碼）；待有 FEM 場後，比照 long2016 的
`current_base`（電荷/電流模型，`USE_BIAS` + variable projection）與 `voltage_base`（Hall-sensor）填入。

內層架構（同 long2016）：`code/{main,main_function,function,plot}/ + data/ + figures/ + results/`。
- `main.m` 只擬合 + 存 `.mat`；產 PDF 的 `emit_*.m` 放 `code/main_function/`。
