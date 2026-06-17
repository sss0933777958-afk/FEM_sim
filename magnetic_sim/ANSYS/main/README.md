# magnetic_sim/ANSYS/main/ — 活躍設計工作根（4-pole MEMS Quadrupole，原 kuo/）

這是目前**唯一活躍**的設計工作根（ANSYS 求解器）。

**完整導覽與規則見本夾的 [`CLAUDE.md`](CLAUDE.md)** —— 含 9 夾架構地圖、資料流 pipeline（CAD_model→IGES→IGES_converted→apdl→ANSYS_data→matlab→MATLAB_data+figures）、resolver 用法（`common/{ansys_path,matlab_path}.m`），以及 **🎨 繪圖腳本規則** 與 **🔒 不擅自更動檔案架構鐵則**。

**9 個頂層夾**（各有自己的 README）：
`ANSYS_data/`（FEM .dat/.db）、`MATLAB_data/`（分析 .mat）、`matlab/`（分析碼）、`apdl/`（APDL 腳本）、`CAD_model/`（SolidWorks/STEP 原檔）、`IGES/` + `IGES_converted/`（幾何匯出）、`doc/`（LaTeX/PDF + workflows SOP）、`.claude/`（本地設定）。

**相關**：repo 級規則見 `../../../.claude/rules/`；repo 總覽見 `../../../CLAUDE.md`。
