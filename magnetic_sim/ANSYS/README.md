# magnetic_sim/ANSYS/ — ANSYS 求解器子層

**用途**：磁學模擬中以 ANSYS MAPDL（APDL magnetostatic FEM）求解的部分。分活躍設計與歸檔兩層。

**內容**：
- `main/` — ★ 目前唯一活躍設計：4-pole MEMS Quadrupole（原 `kuo/`）。完整 pipeline（CAD → IGES → APDL → ANSYS_data → matlab → MATLAB_data → doc）；有自己的 `CLAUDE.md`（資料夾地圖 + 資料流 + 繪圖規則）與各子夾 README。
- `backup/` — 歸檔（非活躍）設計：`hexapole-long2016/`、`hung/`。

**相關**：工作目錄全域規則見 `main/CLAUDE.md` 與 `../../.claude/rules/main-workspace.md`、`main-workflows.md`；歸檔說明見 `backup/README.md`；COMSOL 部分見 `../COMSOL/README.md`。
