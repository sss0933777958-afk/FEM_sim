# ANSYS_data/long2016_hexapole_halfcut/csv/ — sensor-local mesh CSV（gap_200um basegap 網格）

**用途**：Hall-sensor 局部 FEM 網格連接性 CSV，供 `extract_Vmat_interp.m` 對 `coilN/gap_200um` 場做真-FEM tet 重心內插（在 6 顆 sensor 取樣圓柱內撒點）。

**內容**：
- `sensor_local_nodes.csv` — `nodeID, x, y, z`[m, ANSYS 框]（15749 節點）。
- `sensor_local_elems.csv` — `elemID, n1..n8`（SOLID96 8 槽，77929 elems）。

**來源**：`apdl/long2016_hexapole_halfcut/postproc/MT_Export_SensorLocalMesh_basegap.txt` RESUME **`db/mesh_graded_basegap/mesh_graded.db`**（gap_200um 真實 200µm air-slab 網格，706k 節點）→ NSEL 6 sensor ±1mm box → CDWRITE。**node ID 與 `coilN/gap_200um/*.dat` 完全對齊**（座標殘差 ~1e-14 m）。

**⚠ 為何在這**：gap_200um 資料在 basegap 網格（706k），舊 `data/mesh/graded/csv`（plain-graded 636k）node ID 對到差 ~20mm 的錯位置 → V 矩陣全毀。故 gap_200um 專用此對齊 CSV（使用者 2026-07-10 指定放 `ANSYS_data/<model>/csv/`）。兩包 `main.m` 的 `mesh_csv_dir` 指到本夾。

**相關**：`data/mesh/{graded,standard}/csv/`（其他 variant 的 CSV，非 gap_200um 用）、`MT_Export_SensorLocalMesh_basegap.txt`、memory `long2016-hall-sensor-base`。
