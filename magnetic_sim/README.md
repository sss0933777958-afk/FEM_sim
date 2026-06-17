# magnetic_sim/ — 磁學模擬類別

**用途**：`FEM_sim/` 下的磁學（magnetostatic / eddy-current）模擬類別，收磁鑷 hexapole / quadrupole 相關工作。目前是 `FEM_sim/` 唯一的物理類別；未來其他 FEM 類別（靜電 / 結構 / 熱）會與本夾**並列**為兄弟資料夾。

**內容**：依**求解器**分子層：
- `ANSYS/` — ANSYS MAPDL 求解器子層（含 `main/` 活躍設計 + `backup/` 歸檔）。
- `COMSOL/` — COMSOL 求解器子層（與 ANSYS/ 並列），目前存 `.mph` 模型。

**相關**：repo 總覽見 `../CLAUDE.md` / `../README.md`；ANSYS 活躍設計導覽見 `ANSYS/main/CLAUDE.md`。
