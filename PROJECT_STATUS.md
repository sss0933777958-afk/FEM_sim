# PROJECT_STATUS — FEM_sim 現況快照

> 本檔是**現況/交接快照**。詳細工作規則見 `magnetic_sim/ANSYS/main/CLAUDE.md` 與 `.claude/rules/`；
> 各層用途見每個資料夾的 `README.md`。最後更新：2026-06-17。

## 1. 身分 / GitHub
- **FEM_sim** = 通用 FEM 模擬容器（git repo 根 = `G:\my_workspace\code\FEM_sim\`）。
- 已推上 GitHub：**`https://github.com/sss0933777958-afk/FEM_sim`（Public）**，分支 `main`。
- remote：只有 **`origin` = `sss0933777958-afk/FEM_sim`**（`main` 追蹤它）。舊的 `kevinfan100/magnetic-tweezers-sim` 連結已移除。
- **歷史已重置為全新獨立專案**（2026-06-17）：砍掉繼承自 `magnetic-tweezers-sim` 的 53 個舊 commit，重開**單一 Initial commit**，作者只剩 sss。舊歷史備份在本機分支 `backup/old-history-magnetic-tweezers` + tag `pre-squash-2026-06-17`（未推遠端），且原樣仍存在 kevinfan100 的舊 repo。
- 大型 COMSOL `.mph`（約 2.5GB）**只留本地、被 gitignore**（單檔 >100MB，GitHub 不收）。

## 2. 現況結構
```
FEM_sim/                         通用 FEM 容器（CLAUDE.md / README / .claude/rules/）
└── magnetic_sim/                磁學模擬類別
    ├── ANSYS/                   ANSYS 求解器子層
    │   ├── main/                ★ 活躍設計：4-pole MEMS Quadrupole（原 kuo/）
    │   │   ├── CAD_model/ IGES/ IGES_converted/   幾何（CAD→IGES→轉換）
    │   │   ├── apdl/            APDL 腳本（geom/sim/postproc）
    │   │   ├── ANSYS_data/      FEM 輸出 .dat/.db（gitignore 重產物）
    │   │   ├── matlab/          MATLAB 分析（<model>/<功能組>/code，含 common/ resolver）
    │   │   ├── MATLAB_data/     分析成果 .mat/.csv
    │   │   └── doc/             LaTeX/PDF + workflows SOP
    │   └── backup/              歸檔（非活躍）：hexapole-long2016、hung
    └── COMSOL/mph/              COMSOL .mph 模型（gitignore，本地）
```
- 每一層資料夾都有 `README.md`（根到葉，共 ~168 份在 repo）。
- **resolver**（路徑解析，勿硬寫絕對路徑）：`magnetic_sim/ANSYS/main/matlab/<model>/common/{ansys_path,matlab_path}.m`。

## 3. 演進史（精要）
`magnetic-tweezers-sim` →(改名)→ `FEM_sim` →(加分層)→ `magnetic_sim/` →(加求解器層)→ `magnetic_sim/ANSYS/`、`magnetic_sim/COMSOL/`；
設計 `kuo` →(改名)→ `main`；`hung`、`hexapole-long2016` → `ANSYS/backup/`；刪除根層 `docs/`/`references/`/`studies/` 整套 scaffolding；全 repo 逐層補 README；
最後**重置 git 歷史**（2026-06-17）：丟掉繼承自 magnetic-tweezers-sim 的舊 commit 鏈，以單一 Initial commit 作為全新獨立專案起點。

## 4. Standing rules（長期遵守，詳見 main/CLAUDE.md + memory）
1. **不擅自更動檔案架構**（移動/改名/刪除/新建資料夾/重組）—— 要動先問；改檔內文不受限。
2. **改動同步「受影響的」README**（動到結構才連上層索引 + 架構圖）。
3. **繪圖腳本規則**：畫圖前先確認功能組 → 不屬於就先問再開新組 → 每個圖/任務只一支腳本原地改 → 定案後才存最終圖（前面只 MCP preview）→ 場圖一律真實 FEM 節點、不內插。
4. 回覆一律**繁體中文**（技術名詞用英文）。
5. 大型二進位（`.mph` 等）gitignore、本地保留；不上 GitHub。

## 4b. 本期工作（2026-06-17 後續：Hall-sensor 求 d + 視覺化）
- `magnetic_sim/ANSYS/main/matlab/long2016_hexapole_halfcut/` 下，使用者把 `Calibration using FEM modeling/` 改名 **`Calibration_using_FEM_modeling/`**，並 `fix_l→fix_dir`、`no_fix_l→no_fix_dir`、`sensor_d→Hall_sensor_base_fix_dir`。
- **新建 `Hall_sensor_base_no_fix_dir/`**：18-param bias 版求 Hall-sensor `d`（載 `calib_bias.mat`、actuator 框、`build_A`/`sensor_residual_bias`）；ℓ̂=0.857、sensor RMSE 8.39%、R_a 6.85e8；存 `calib_sensor_d_no_fix_dir.mat`。fix 版：ℓ̂=0.856、RMSE 15.52%、R_a 7.79e8。
- results/ 去冗餘：只留 `d_final.tex`+`KH_final.tex`（`d_v2` 無增益版＝d_final×g_H 已不輸出），且 results/figures 不再多包子夾。
- **sensor off-diag 對極微正號** 經查證＝對極近零耦合(~0.5%對角)非 bug（extract_Vmat≡既有 B_S 法 max|Δ|=2e-16）。
- **磁路視覺化** `code/plot/plot_sensorBcircuit_P1exc(pole_i)`（通用）+ `plot_P2sensor_Braw_P1exc`：真實節點、all-source、world-up、cone 真輪廓+n+；圖 `figures/{P2,P4}sensor_Braw_P1exc.png`。
- **reorg 失效路徑全 repo 補正**：`.m` 115 處/75 檔、`.md` 104 處/25 檔（`magnetic_sim\hexapole-long2016`→`\ANSYS\backup\...`、`\hung`→`\ANSYS\backup\hung`）。詳見 memory [[long2016-hall-sensor-base]]、[[fem-sim-github-status]]。

## 5. 未決 / 待辦（之後可處理）
- **命名兩套**：同一 Long Fei 模型在 `IGES_converted/` 叫 `long_fei`，但 `IGES/`、`ANSYS_data/`、`apdl/` 叫 `long2016_hexapole_halfcut`。
- 本機備份分支 `backup/old-history-magnetic-tweezers` + tag `pre-squash-2026-06-17` 保留舊歷史；確定不需要救回後可刪。
- `.txt` APDL 內舊路徑指標未補（待定是否需要）；`Hall_sensor_base_fix_dir` P2 專用繪圖腳本與通用版並存（待定整併）。

> 已解決：`settings.local.json` 殘留歷史、remote 命名（`gh`→`origin`）、`.m`/`.md` 失效 addpath —— 皆已處理。

## 6. 驗證指紋（沿用）
- fix_l 校正 R=150：`ell=0.856 mm / gB=8.43e-3 / NRMSE 1.23%`；no_fix_l：`ell=0.857 / gB=9.50e-3`。
- ANSYS：`G:\ANSYS Inc\v252\ansys\bin\winx64\MAPDL.exe`；MATLAB R2025b（`-batch` 收尾有已知 std::terminate teardown crash，output 先印出不影響）。
