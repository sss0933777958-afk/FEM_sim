# code/function/ — 與 main 無掛鉤的功能函式（工具 / 診斷）

**用途**：放置**不被 `main.m` 呼叫**的功能函式（獨立工具、診斷、替代方法、實驗性算法）。
main 校正流程用到的函式一律在 `../main_function/`；本夾的函式只被 `../plot/` 腳本或使用者手動呼叫。

**內容**：
  - `compute_vmatrix.m`
  - `extract_Vmat.m`
  - `extract_Vmat_elemB.m`
  - `extract_Vmat_interp_center.m`
  - `make_scaled_coil1.m`

**相關**：`../main_function/README.md`（main pipeline 函式）、`../README.md`。
