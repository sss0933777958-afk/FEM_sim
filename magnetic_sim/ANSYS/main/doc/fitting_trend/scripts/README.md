# doc/fitting_trend/scripts/ — R 掃描趨勢 LaTeX 原稿 + 產生器

**用途**：存放 R 掃描趨勢分析的 **LaTeX 原稿（.tex）** 與 **LaTeX 表格 / 趨勢產生器（.m）**。
**內容**：`note.tex`（趨勢說明 TikZ 流程圖原稿）；`gen_KI_trend_per_radius.m`、`gen_KI_alln_per_radius.m`、`gen_nofixl_per_radius.m`（逐半徑提取 K̂_I / 參數）；`gen_KI_trend_latex.m`、`gen_KI_trend_N1000_latex.m`（輸出 .tex 趨勢表，N1000 = 每極 1000 取樣點變體）。
**命名 / 慣例**：葉夾同放 .tex 與 .m；`note.pdf` 編譯產物放上一層 ../note.pdf（本 topic 無 pdf/ 葉夾）；`per_radius`=逐半徑、`nofixl`=不固定 ℓ、`N1000`=取樣點數變體。
**相關**：見 ../README.md、../../../CLAUDE.md。
