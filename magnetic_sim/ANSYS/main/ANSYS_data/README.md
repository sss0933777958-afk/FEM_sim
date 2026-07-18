# ANSYS_data/ — FEM 求解輸出（原始資料）

ANSYS MAPDL 求解後的**原始 FEM 資料**。與 MATLAB 分析成果（`../MATLAB_data/`，`.mat`）分開。

## 結構：`<model>/data/<variant>/coil<N>/`（2026-07-14 全包統一為 `<variant>/coil<N>`）

FEM 場 `.dat` 一律 `<model>/data/<variant>/coil<N>/coil<N>_{coord,bfield}_{all,wp}.dat`
（**2026-07-14 使用者拍板：三包全翻成 `<variant>/coil<N>`，不再 `coil<N>/<variant>`**）。

```
ANSYS_data/
├── long2016_hexapole_halfcut/   ← 主力（Long Fei 下極半切六極）
│   ├── data/<variant>/coil1..6/ ← standard(coil2-6)/no_gap/gap_calibrate/gap_200um/gap_300um/graded/protgap/gap400
│   ├── data/{singlepole, coil1/singlepole, mesh}/  ← 非 coil-indexed（未翻轉，維持原位）
│   ├── db/ , csv/               ← .db 模型/mesh；csv=<variant>/*.csv（無 coil 層）
│   └── RESULTS_MAP.md           ← ★ 哪個 dir 是哪份結果的權威清單
├── hung_hexapole/               ← data/<variant>/coil1..6（variant∈{gap_200um,no_gap}）+ db/ + RESULTS_MAP.md
└── NTU_hexapole/                ← data/<類別>/coil1（類別∈{singlepole,upper_assembly}）+ db/
```

## 檔案類型
- `*.dat` — 抽出的場資料（節點座標 + Bx/By/Bz；MATLAB 讀這個）。
- `*.db` / `*.cdb` — ANSYS 模型 / mesh（GUI resume、重抽物理量用）。
- `RESULTS_MAP.md` — 每個 case dir 的物理意義 + 期望指紋（節點數、|B|max）。

## 規則（重要）
- **讀結果前先查 `<model>/RESULTS_MAP.md`**，並照 `.claude/rules/result-read-safety.md` 三層防呆（回報路徑+指紋→核對→不對就停）。`coilN` vs `coilN_gap_200um` 節點數相同，只能靠 |B| 區分。
- **gitignore**：`.dat/.db/.cdb/.rmg/.esav/...` 等 FEM 重產物不進 git（見 repo `.gitignore`）。
- 讀取一律用 resolver `ansys_path('<model>','coilN',...)`（在 `../matlab/<model>/common/`），不要硬寫絕對路徑。
- 清理副產物前**必讀** `.claude/rules/sim-cleanup.md`（預設 half-clean）。
