# apdl/long2016_hexapole_halfcut/postproc/dump/ — 從存好的 .db 重抽場到 .dat

**用途**：RESUME 既有 master `.db`、SET,LAST 取既有解，對指定空間區域抽座標 + B 場 dump 成 `.dat`。用來補抽當初 sim 沒輸出的區域（例如 graded sim 只 dump 了 2mm WP sphere，需要更廣的磁路區域）。

**內容**：
- `MT_Dump_P1_circuit.txt` — 從 P1 graded master `.db` 重抽整個磁路 y≈0 側視 slab（|y|<1mm、x∈[-5,60]mm、z∈[-17,5]mm，APDL frame），給 P1 磁路側視圖（cone→arm→post→yoke→coil）。

**資料來源 / 流向**：input = `ANSYS_data/long2016_hexapole_halfcut/coil1/graded/3DMTmagneticfield.db`；output = 同 CWD 的 `coil1_coord_circuit.dat` + `coil1_bfield_circuit.dat`。

**命名 / 慣例**：`MT_Dump_<pole>_<region>.txt`；只 RESUME + POST1，不 /SOLU、不改幾何；`NLIST,...,COORD` 抽座標、`PRNSOL,B,COMP` 抽分量；改動標 `[ADDED]`/`[MODIFIED]`。

**相關**：見 `../README.md`、`doc/workflows/apdl-postproc.md`、`.claude/rules/result-read-safety.md`。
