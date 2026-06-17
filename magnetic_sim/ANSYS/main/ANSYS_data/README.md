# ANSYS_data/ — FEM 求解輸出（原始資料）

ANSYS MAPDL 求解後的**原始 FEM 資料**。與 MATLAB 分析成果（`../MATLAB_data/`，`.mat`）分開。

## 結構：`<model>/<case>/`

```
ANSYS_data/
├── long2016_hexapole_halfcut/   ← 主力（Long Fei 下極半切六極）
│   ├── coil1 … coil6/           ← 每顆 pole 自激解（baseline）
│   ├── mesh/                    ← mesh-only .db/.cdb
│   └── RESULTS_MAP.md           ← ★ 哪個 dir 是哪份結果的權威清單
├── kuo_quadrupole/              ← coil1..6 + mesh（4-pole 四極）
└── zhang_quadrupole/            ← coil1..6 + mesh
```

## 檔案類型
- `*.dat` — 抽出的場資料（節點座標 + Bx/By/Bz；MATLAB 讀這個）。
- `*.db` / `*.cdb` — ANSYS 模型 / mesh（GUI resume、重抽物理量用）。
- `RESULTS_MAP.md` — 每個 case dir 的物理意義 + 期望指紋（節點數、|B|max）。

## 規則（重要）
- **讀結果前先查 `<model>/RESULTS_MAP.md`**，並照 `.claude/rules/result-read-safety.md` 三層防呆（回報路徑+指紋→核對→不對就停）。`coilN` vs `coilN_gap200um_mueq` 節點數相同，只能靠 |B| 區分。
- **gitignore**：`.dat/.db/.cdb/.rmg/.esav/...` 等 FEM 重產物不進 git（見 repo `.gitignore`）。
- 讀取一律用 resolver `ansys_path('<model>','coilN',...)`（在 `../matlab/<model>/common/`），不要硬寫絕對路徑。
- 清理副產物前**必讀** `.claude/rules/sim-cleanup.md`（預設 half-clean）。
