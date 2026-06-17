# matlab/ — MATLAB 分析程式碼

讀 `../ANSYS_data/` 的 `.dat` 做擬合 / 算矩陣 / 校正 / 畫圖；成果 `.mat` 寫到 `../MATLAB_data/`，圖寫到各功能組 `figures/`。

## 結構：`<model>/<功能組>/`
```
matlab/
├── long2016_hexapole_halfcut/      ← 主力
│   ├── common/                     ← ★ resolver：ansys_path.m / matlab_path.m
│   ├── Calibration using FEM modeling/{fix_l,no_fix_l}/  ← 單一主程式：code/main/main.m
│   ├── fixl_fit/ bias_fit/ bs_matrix/ sensor_d/ validation/  ← 多腳本：code/scripts/
│   └── field_viz/ sensor_placement/                          ← 純繪圖：code/plot/
├── kuo_quadrupole/    ← 空 placeholder
└── zhang_quadrupole/  ← 空 placeholder
```
每個功能組視性質含 `code/{main|scripts|plot}/` + `figures/`（圖）+ `results/`（auto-gen `.tex`）。

## Resolver（不要硬寫絕對路徑）
`<model>/common/`：
- `ansys_path('<model>','coilN',...)` → 讀 FEM `../ANSYS_data/<model>/...`
- `matlab_path('<model>','<功能>',...)` → 讀寫 `../MATLAB_data/<model>/<功能>/...`

## 🎨 繪圖腳本規則（重點，完整見 `../CLAUDE.md`）
1. 畫圖前**先確認屬哪個功能組**；不屬於任何組 → 先問再開新組（`<model>/<新組>/code/plot/` + `figures/`）。
2. **每個圖/繪圖任務只一支腳本**，原地反覆改到定案；**定案前不另開新腳本**、不存最終圖（先用 MCP preview）。
3. 使用者沒說「新增」就不開第二支。
4. 場圖一律畫**真實 FEM 節點原值**，不內插（除非明確要求並標示）。

## 規則
- 擬合電流要對齊 FEM 激發電流（1A），見 `.claude/rules/fit-current-matches-sim.md`。
- 工作放置慣例（哪種產物去哪）見 `.claude/rules/main-workspace.md`。
