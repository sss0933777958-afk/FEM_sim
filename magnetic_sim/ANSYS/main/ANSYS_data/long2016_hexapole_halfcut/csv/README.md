# ANSYS_data/long2016_hexapole_halfcut/csv/ — sensor-local mesh CSV（gap_200um basegap 網格）

**用途**：Hall-sensor 局部 FEM 網格連接性 CSV，供 `extract_Vmat_interp.m` 對 `coilN/gap_200um` 場做真-FEM tet 重心內插（在 6 顆 sensor 取樣圓柱內撒點）。

**內容**：
- `sensor_local_nodes.csv` — `nodeID, x, y, z`[m, ANSYS 框]（15749 節點）。
- `sensor_local_elems.csv` — `elemID, n1..n8`（SOLID96 8 槽，77929 elems）。

**來源**：`apdl/long2016_hexapole_halfcut/postproc/MT_Export_SensorLocalMesh_basegap.txt` 匯出時 RESUME **`db/mesh_graded_basegap/mesh_graded.db`**。**⚠ 本子夾（頂層 csv/）的 706k CSV 是 gap_200um 時代匯的**（那時 db=706k real 200µm slab）；**該 db 已於 2026-07-12 覆寫成 1.21M gap_300um 網格**——故此 deck 現在重跑會產 gap_300um CSV（輸出已改指 `csv/gap_300um/`，見下）。頂層 706k CSV 檔留著給 gap_200um 用（node ID 對齊 `coilN/gap_200um/*.dat`，殘差 ~1e-14 m）。

**⚠ 為何在這**：gap_200um 資料在 basegap 網格（706k），舊 `data/mesh/graded/csv`（plain-graded 636k）node ID 對到差 ~20mm 的錯位置 → V 矩陣全毀。故 gap_200um 專用此對齊 CSV（使用者 2026-07-10 指定放 `ANSYS_data/<model>/csv/`）。兩包 `main.m` 的 `mesh_csv_dir` 指到本夾（gap_300um 時自動改指 `csv/gap_300um/`）。

## `gap_300um/`（2026-07-12 新增）
- `sensor_local_{nodes,elems}.csv`（15507 節點 / 76384 SOLID96 elems）——從 **gap_300um 的 1.21M conformal real-slab basegap 網格** 匯出（同 deck、RESUME 現在的 `mesh_graded_basegap.db`）。**node ID 與 `coilN/gap_300um/*.dat` 對齊**（座標殘差 5e-14 m 已驗）。
- 供 `Hall_sensor_base_{fix,no_fix}_dir/code/main/main.m`（`VARIANT='gap_300um'`）的 `extract_Vmat_interp`；`mesh_csv_dir` 於 gap_300um 時自動指到本子夾。頂層 706k 與本 1.21M 兩套 **不可互換**（不同網格）。

**相關**：`data/mesh/{graded,standard}/csv/`（其他 variant 的 CSV，非 gap_200um 用）、`MT_Export_SensorLocalMesh_basegap.txt`、memory `long2016-hall-sensor-base`。
