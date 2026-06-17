# doc/ — 文件（LaTeX 推導 / 報告）+ 流程 SOP

技術推導、報告 LaTeX 原稿與編譯 PDF，以及 Claude 可 follow 的跨主題流程 SOP。

## 結構
```
doc/
├── workflows/                ← ★ 流程 SOP（cad-export / apdl-fem-run / charge-model-fit / …）
│                               入口見 workflows/README.md
├── charge_model_fitting/     ← 電荷模型擬合推導（fitting_derivation / general / <model>/）
├── Solve_B_matrix/<model>/   ← B_S 矩陣推導
├── Br_analysis/<model>/      ← Br 分析
├── pole_WP_distance/<model>/ ← pole–WP 距離
├── fitting_trend/scripts/    ← R 掃描趨勢
├── error_definition/{pdf,scripts}
├── magenetic_flux_integral/{pdf,scripts}
└── ur_iteration/{pdf,scripts}
```

## 檔案類型
- `*.tex`（原稿）/ `*.pdf`（編譯產物）/ `*.m`（LaTeX 表格產生器）/ `*.md`（SOP、說明）。

## 規則
- **跑流程 / 操作前先看 `workflows/README.md`** 找對應 SOP；觸發詞 → SOP 對照見 `.claude/rules/main-workflows.md`。
- LaTeX 原稿放 `<analysis>/<topic>/scripts/`，編譯 PDF 放 `<analysis>/<topic>/pdf/`（analysis-first schema，見 `main-workspace.md`）。
- 中文報告需要圖時 **reference** `../matlab/<model>/<功能組>/figures/<file>.png`，不要把圖另存到 git root 外。
