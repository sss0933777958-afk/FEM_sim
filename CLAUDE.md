# FEM Simulation Workspace

`FEM_sim/` 是**通用 FEM 模擬容器**。磁學相關的模擬全部收在 `magnetic_sim/` 這一類底下，再依**求解器**分子層：目前為 `magnetic_sim/ANSYS/`（未來可並列 `COMSOL/` 等），其下 `main/` 為活躍設計（4-pole MEMS Quadrupole，原 `kuo/`），`backup/` 歸檔非活躍設計（`hexapole-long2016/`、`hung/`）；未來其他模擬類別（如靜電 / 結構 / 熱）會與 `magnetic_sim/` **並列**為 `FEM_sim/` 的兄弟資料夾。

## Quick Triggers
- 當使用者說「**建 hexapole**」時，參照 `.claude/rules/hexapole-build.md` 執行建模流程
- 當工作涉及 `magnetic_sim/hung/` 目錄時，參照 `.claude/rules/hung-docs.md` 讀取必要文件
- 當工作涉及 `magnetic_sim/ANSYS/main/` 目錄（cwd 在 magnetic_sim/ANSYS/main/、討論 Quadrupole、編輯 magnetic_sim/ANSYS/main/* 檔案）時，參照 `.claude/rules/main-workspace.md`，所有新產物寫進 magnetic_sim/ANSYS/main/ 下對應子目錄，不得寫到其他設計目錄或外部路徑
- 當使用者要求 SolidWorks 出檔／解析 STEP／寫 MT_Geom／檢查模型 等 kuo 流程操作時（自然語觸發，例「出 STEP」「解析 STEP」「建 APDL 幾何」「檢查模型」），參照 `.claude/rules/main-workflows.md` 啟動對應 SOP（`magnetic_sim/ANSYS/main/doc/workflows/`）
- 當工作涉及「**清 sim 副產物 / 清 ANSYS results / 磁碟滿 / 整理 result dir / cleanup sim**」等清理動作時，**動手前必須先完整讀 `.claude/rules/sim-cleanup.md`**；該規則寫死「6 項不可影響工作 + 2 項不可失去能力」criteria、強制 dry-run、預設 half-clean（保 `.db` + 主 `.rmg`）、`--full` 要使用者明確批准；helper script `magnetic_sim/ANSYS/main/apdl/common/clean_sim_dir.sh`
- 當使用者貼一個 `IGES_converted/<...>.iges` 或 `IGES/<...>.iges` 路徑時，**Claude 必須從路徑識別這是哪個物理模擬模型**（Long Fei 半切六極 / Kuo Quadrupole / Zhang Quadrupole / ...），把該模型當作後續對話的「當前討論模型」，**不要問「這是哪個模型？」**；模型清單見 `.claude/rules/iges-model-id.md`
- 當工作涉及「**改 ANSYS 幾何 / 改 mt_constants / 寫新 APDL 幾何腳本 / 對齊 CAD / ANSYS 跟 CAD 不一致**」等動作時，**動手前必須先完整讀 `.claude/rules/ansys-cad-alignment.md`**；該規則寫死「CAD STEP/IGES 是 source of truth、ANSYS 數值必對齊 CAD、改前必量 CAD、不一致必通報使用者由其拍板、預設 Path A(改 ANSYS)」
- 當工作涉及「**跑 COMSOL / LiveLink / mphserver / 連 COMSOL server / 跑 .mph / COMSOL 頻掃**」等動作時，參照 `.claude/rules/comsol-livelink.md`（跨設計通用）；該規則寫死成功連線法 = 「獨立啟 `comsolmphserver.exe` + 另一個 `matlab.exe -batch` 內 `mphstart(2036)`」兩個 process，**不要用整合式 `comsolmphserver matlab`**（Win + R2025b 壞）；一鍵 launcher `magnetic_sim/ANSYS/main/comsol/kuo_quadrupole/run_matlab_with_comsol.ps1`
- 當工作涉及「**配 FEM 場的擬合 / 改 I_actual·I_in / 寫新 charge fit / 用 fit 參數預測不同電流**」時，**動手前必須先讀 `.claude/rules/fit-current-matches-sim.md`**；該規則寫死「模型電流 I 必須等於 FEM 激發電流（目前 1A），不可塞操作電流 0.6A（會把 1/0.6 假因子灌進 gB、預測需補正）；例外＝用 0.6A 把 1A 場縮到操作點的 V/V 矩陣」
- 當工作涉及「**讀 / 載入 ANSYS 結果（抽 .dat / import_ansys_data / 載入 coilN / postproc 或算矩陣·fit·畫場圖前載入 result）**」時，**動手前必須先讀 `.claude/rules/result-read-safety.md`**，照三層執行：①讀前回報絕對路徑+dataset+期望指紋，≥2 候選讓使用者選不自己猜 ②讀後核指紋（matched 節點數 + \|B\| max + case_tag，baseline vs `gap200um_mueq` 同節點數只能靠 \|B\| 低 ~30% 區分）對不上就停 ③查 `magnetic_sim/ANSYS/main/ANSYS_data/<topic>/RESULTS_MAP.md` 凌駕 memory（目前已建 long2016_hexapole_halfcut）

## Commands

### ANSYS 可用性

**執行 ANSYS 前必須先確認路徑存在**。本機實際安裝位置：

```
G:\ANSYS Inc\v252\ansys\bin\winx64\MAPDL.exe
```

若路徑不存在（例如換機器或版本），在標準位置（`C:\Program Files\ANSYS Inc\<version>\...`）或其他磁碟搜尋後再跑。

### 典型指令

```bash
ANSYS="G:\ANSYS Inc\v252\ansys\bin\winx64\MAPDL.exe"

# Run single coil (batch mode, no GUI) — run from magnetic_sim/hexapole-long2016/
cd magnetic_sim/hexapole-long2016
"$ANSYS" -b -np 4 -m 24000 \
  -dir "results/coil1" -j "coil1" \
  -i "$(pwd)/apdl/MT_Modeling_Geometry_Meshing_Solving_Coil1.txt" \
  -o "results/coil1/solve.out"

# Run all 6 coils sequentially
for i in 1 2 3 4 5 6; do
  "$ANSYS" -b -np 4 -m 24000 \
    -dir "results/coil${i}" -j "coil${i}" \
    -i "$(pwd)/apdl/MT_Modeling_Geometry_Meshing_Solving_Coil${i}.txt" \
    -o "results/coil${i}/solve.out"
done
```

## Architecture
```
FEM_sim/                   Git root — 通用 FEM 模擬容器
├── README.md                            Project overview
├── CLAUDE.md                            This file
├── .gitignore                           Excludes ANSYS outputs
├── .claude/rules/                       Path-scoped editing rules
├── magnetic_sim/                        磁學模擬類別（目前唯一類別）
│   └── ANSYS/                           ANSYS 求解器子層（未來可並列 COMSOL/ 等）
│       ├── main/                        ★ 活躍設計：4-pole MEMS Quadrupole (Harrison-style；原 kuo/)
│       │   ├── apdl/{geom,sim,postproc,sweep}/  APDL scripts
│       │   ├── matlab/                  MATLAB analysis (含 common/ resolver)
│       │   ├── ANSYS_data/<model>/<case>/  FEM .dat/.db (gitignored)
│       │   ├── MATLAB_data/<model>/<fn>/   MATLAB outputs (.mat/.csv)
│       │   ├── figures/                 All figures (incl. reports)
│       │   ├── IGES/ + IGES_converted/  Geometry exports (must sync)
│       │   ├── CAD/                     SolidWorks/STEP originals
│       │   ├── comsol/ + mph/          COMSOL LiveLink scripts + .mph models
│       │   ├── semulator/               SEMulator process flow
│       │   └── reference/               Paper PDFs (design-specific)
│       └── backup/                      歸檔（非活躍設計）
│           ├── hexapole-long2016/       Long 2016 dissertation hexapole design
│           └── hung/                    Hung hexapole design (build workflow)
└── (future: electric_sim/ 等其他 FEM 類別，與 magnetic_sim/ 並列)
```

## Hexapole Design Constraints (Mandatory)

These constraints apply to ALL hexapole designs in this repo. They are non-negotiable.

1. **Orthogonal pair axes**: 3 opposing pole pairs (P1-P2, P3-P4, P5-P6) must have mutually perpendicular connecting lines
2. **Tips on common sphere**: All 6 pole tips at distance R_norm from WP center (R_norm is adjustable)
3. **60-degree azimuthal offset**: Upper layer rotated 60 deg relative to Lower layer
4. **alpha = arctan(sqrt(2)) = 54.74 deg is FIXED**: derived from constraints 1-3, not a free parameter
   - `R_norm_xy = R_norm * sqrt(2/3)` and `R_norm_z = R_norm / sqrt(3)` — these formulas are locked
   - Lower poles at 0, 120, 240 deg; Upper poles at 60, 180, 300 deg

## Rules
- 6 Coil scripts are synchronized: only `CURR_ARRAY` values differ (one coil = 1, rest = 0)
- All code comments in English; explanations to user in Traditional Chinese
- Mark all APDL changes with `[ADDED]` or `[MODIFIED]` comments
- Always verify `D,ALL,MAG,0` boundary condition exists before `/SOLU`
- Preserve original commented-out code (prefixed `!****`) unless asked to remove
- Use tab indentation matching original style
- Use dissertation notation (B, Phi, q, K_I, rho, R_a, g_I, etc.) in all discussion and code comments
- Always refer to poles by paper name (P1-P6); mention APDL index only when editing APDL code

## Figure Production
- **Never generate figures without discussion first.** Before producing any figure:
  1. **Content**: Discuss what to show — which data, axes, normalization, range
  2. **Style**: Discuss visual details — font, title, legend, colors, line thickness
  3. **Preview**: Use MATLAB MCP to render a draft, review together, iterate
  4. **Finalize**: Only save to `figures/` after user confirmation
- Apply consistent figure style across the project
- **資料來源預設用真實模擬節點**：場圖一律直接畫 FEM 節點原值（節點實際位置 + 節點 Bx/By/Bz），**不要用 scatteredInterpolant / 規則格點內插，除非使用者明確要求內插**。因使用者要求而內插時，**必須在回覆/圖說明白標示「此圖為內插」**，絕不可把內插當 raw 呈現。可讀性需減量時，對節點做抽樣（每格挑最近 y=0 平面的節點）——那仍是節點原值、非內插。見 memory `plot-real-nodes`。

## Prohibitions
- NEVER commit ANSYS output files (*.rst, *.db, *.full, etc.)
- NEVER change geometry parameters without explicit user approval
- NEVER modify element types or material properties without approval
- NEVER remove boundary condition section (`[ADDED]` block near line 500)
- NEVER 跑任何 sim 清理（rm intermediates / rm result dir）前未先讀 `.claude/rules/sim-cleanup.md` 全文 — 該規則寫死「6 項不可影響工作 + 2 項不可失去能力」criteria；違反 = 違規
- NEVER 用 `--full` 模式清 sim 副產物 unless 使用者**明確同意**（預設一律 half-clean，保 `.db` + 主 `.rmg`）
- NEVER 繞過 helper `magnetic_sim/ANSYS/main/apdl/common/clean_sim_dir.sh` 直接手刻 `rm` ANSYS 檔（會跟規則的保留清單不同步）
- NEVER 改 ANSYS 幾何尺寸或 mt_constants 前未先量對應 CAD（SolidWorks STEP/IGES）並比對；發現不一致**不可自己選一個值**，必須通報使用者由其拍板（per `.claude/rules/ansys-cad-alignment.md`）
- NEVER change alpha (54.74 deg) or the R_norm_xy / R_norm_z formulas
- NEVER produce a pole configuration that violates pair-axis orthogonality
- NEVER 用 scatteredInterpolant / 格點內插畫場圖，**除非使用者明確要求內插**；也 NEVER 把內插圖當成 raw／節點原值呈現（預設一律真實模擬節點，per Figure Production）

## Notation Standard
All symbols and terms follow Fei Long's 2016 dissertation. See the full glossary:
- `magnetic_sim/hexapole-long2016/docs/notation-glossary.md` - **canonical** symbol/term mapping

Key conventions:
- Use **paper pole names** (P1-P6) in all user-facing text, figures, and discussion
- APDL coil indices (1-6) only in APDL code and raw data context
- Mapping: APDL {1,2,3,4,5,6} = Paper {P1,P3,P6,P5,P2,P4}
- Physical quantities use dissertation symbols: B, Phi, q, K_I, R_hat, L_i, rho, R_a, g_I, N_c
- Two meanings of rho: physical (500 um) vs fitted (900 um) — always clarify which
- Units: ANSYS outputs Tesla; figures use mT for WP region; dissertation Fig. 2.4 uses Gauss

## Detailed Docs
- `magnetic_sim/hexapole-long2016/docs/fitting-methods.md` - **[A]→[J]→[B-6x] fitting methods, [B-6x] is final**
- `magnetic_sim/hexapole-long2016/docs/model-validation.md` - APDL vs dissertation comparison
- `magnetic_sim/hexapole-long2016/docs/notation-glossary.md` - unified notation, dissertation alignment
- `magnetic_sim/hexapole-long2016/docs/coil-winding-sign-convention.md` - pole polarity & coil_sign correction
- `magnetic_sim/hexapole-long2016/docs/charge-model-fitting.md` - point-charge model derivation
- `magnetic_sim/hexapole-long2016/docs/ansys-environment.md` - ANSYS install, batch mode, hardware
- `magnetic_sim/hexapole-long2016/docs/simulation-parameters.md` - geometry, materials, mesh, solver
- `magnetic_sim/hexapole-long2016/docs/workflow.md` - 4-stage simulation-to-publication pipeline
- `magnetic_sim/hexapole-long2016/docs/troubleshooting.md` - known errors and fixes

## Compact Instructions
When context is compressed, preserve:
1. The 6 scripts differ ONLY in CURR_ARRAY (coil N has CURR_ARRAY(N)=1)
2. Boundary condition D,ALL,MAG,0 is mandatory for DSP solver
3. Results go to magnetic_sim/hexapole-long2016/results/coilN/ directories
4. User prefers Traditional Chinese explanations
5. Use paper pole names P1-P6 (not APDL indices) in all discussion
6. Notation follows Long 2016 dissertation — see `magnetic_sim/hexapole-long2016/docs/notation-glossary.md`
7. alpha = 54.74 deg is FIXED for all hexapole designs
