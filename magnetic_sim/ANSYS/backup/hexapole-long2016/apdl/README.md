# magnetic_sim/ANSYS/backup/hexapole-long2016/apdl — APDL 建模／求解／後處理腳本（歸檔/archived）

> **本設計已歸檔（非活躍），保留供參考。** Long Fei 2016 dissertation 六極 hexapole（full）。現役設計請見 `../../../main/`。

**用途**：ANSYS MAPDL 輸入腳本——建幾何、佈 mesh、加 SOURC36 線圈激發、求解（DSP / `magsolv,3`），以及 POST1 抽場到 `.dat`。

**內容**：純文字 `.txt`：
- 6 顆線圈求解腳本 `MT_Modeling_Geometry_Meshing_Solving_Coil1..6.txt`（+ 變體 `Coil5_sph10.txt`、測試 `MT_TestA_ProtH41.txt`）。
- POST1 抽場腳本 `post_extract_coil1..6.txt`（+ `post_extract_coil5_sph10.txt`）。
- 路徑／電路追蹤輔助 `post_path_P2_cone.txt`、`post_trace_circuit.txt`、`post_trace_tips.txt`。

**資料來源 / 流向**：`apdl(本夾) → [ANSYS MAPDL 求解] → results/coilN/*.dat → analysis → data → figures`。求解腳本內 `/CWD` 指向 `results/coilN/`。

**命名 / 慣例**：6 顆求解腳本**僅 `CURR_ARRAY` 不同**（該顆線圈 = 1A、其餘 0A），其餘內容必須同步（線性超疊加 basis）；改動以 `[ADDED]` / `[MODIFIED]` 英文註解標記，保留 `!****` 註解碼；`D,ALL,MAG,0` 邊界為 DSP 求解必要條件。APDL coil 索引 1–6 對應論文極名見 ../README.md 的 Pole Naming Convention 表。

**相關**：見 ../README.md、../docs/（simulation-parameters、ansys-environment、troubleshooting、coil-winding-sign-convention）。
