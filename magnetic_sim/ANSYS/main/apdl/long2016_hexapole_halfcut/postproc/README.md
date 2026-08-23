# apdl/long2016_hexapole_halfcut/postproc/ — 後處理（抽場 / dump）

**用途**：對已求解、存好的 `.db` / 結果做 POST1 後處理，抽出 B 場到 `.dat`。不重跑求解，只 RESUME 既有解再 PRNSOL。

**內容**：
- `dump/` — `MT_Dump_*.txt`：RESUME 存好的 master `.db`，對指定區域（slab / WP / circuit）抽座標 + B 場到 `.dat`。
- `MT_Export_SensorLocalMesh{,_graded,_basegap}.txt` — RESUME 對應網格 `.db`（standard / plain-graded / **basegap**），dump 6 顆 Hall-sensor 附近的局部 SOLID96 air tet 連接性到 CSV（供 `extract_Vmat_interp.m` 做真-FEM tet 重心內插）。**`_basegap`（2026-07-10）= gap_200um 真實 200µm slab 網格（706k），CSV → `ANSYS_data/<model>/csv/`**（舊 `_graded` 的 636k CSV 與 gap_200um 706k 差 ~20mm 不可用）。純網格 export、不求解。

**資料來源 / 流向**：input = `ANSYS_data/long2016_hexapole_halfcut/coilN/<case>/3DMTmagneticfield.db`（已解）；output = 同 CWD 的 `coilN_coord_*.dat` + `coilN_bfield_*.dat`，供 MATLAB 讀。

**命名 / 慣例**：`MT_Dump_*` / `MT_Post_*`；後處理腳本不改幾何/不求解；座標選取用 `NSEL` 範圍框；輸出 .dat 命名含 coil index 與區域 tag（`_wp`/`_circuit`）；改動標 `[ADDED]`/`[MODIFIED]`。

**相關**：見 `../README.md`、`./dump/README.md`、`reference/workflows/apdl-postproc.md`、`.claude/rules/{apdl-editing,result-read-safety}.md`。
