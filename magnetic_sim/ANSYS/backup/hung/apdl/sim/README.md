# magnetic_sim/ANSYS/backup/hung/apdl/sim — 主求解腳本（6-coil simulation）（歸檔/archived）

> **本設計已歸檔（非活躍），保留供參考。** Hung hexapole（build-workflow design）。目前唯一活躍設計為 `../../../main/`。

**用途**：Hung hexapole 主設計（D-shape + 40 µm fillet）的 6 顆 coil 自激磁靜態求解。

**內容**（代表檔）：
- `MT_Hung_Simulate_Coil[1-6]_filleted.txt` — 6 顆極各自激發（P1–P6）。
- `MT_Hung_Simulate_Coil1_filleted_conv.txt` / `_l250*` / `_round_filleted_conv.txt` — NREFINE 收斂版、ℓ=250 µm 與 round-fillet 變體。

**資料來源 / 流向**：讀 `../geom/` 建的幾何 → ANSYS MAPDL 求解 → 輸出 `.db/.rmg` 至 `../../results/coilN/<變體>/`，再由 `../postproc/` 抽 `.dat`。

**命名 / 慣例**：6 腳本共用幾何 / 材料 / mesh / BC，**只差 `CURR_ARRAY` 與輸出 `/CWD`**（coil N 為 CURR_ARRAY(N)=1）；TURNS=70、I=1 A、鋼 μ_r=280 線性、solver magsolv,3（DSP，需 `D,ALL,MAG,0`）。`coilN` = 各極。

**相關**：見 `../README.md`（coil 對照、shared 設定、執行指令）、`../variants/README.md`、上層 `../../README.md`。
