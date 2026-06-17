# doc/fitting_trend/ — 擬合參數隨取樣半徑 R 的掃描趨勢

**用途**：分析點電荷擬合參數（ℓ̂、ĝ_B、K̂_I、NRMSE）如何隨取樣半徑 R 變化 —— R 從 40 µm 連續掃到 500 µm，逐半徑提取參數並找穩定收斂區間（校正半徑選擇）。
**內容**：頂層 `note.pdf`（趨勢說明，TikZ 流程圖）+ `scripts/`（.tex 原稿 + .m LaTeX 產生器，含 fix_l 與 no_fix_l 版、N1000 版）。代表檔：`note.tex`、`gen_KI_trend_per_radius.m`、`gen_KI_trend_latex.m`、`gen_KI_trend_N1000_latex.m`、`gen_KI_alln_per_radius.m`、`gen_nofixl_per_radius.m`。
**命名 / 慣例**：analysis-first schema：`<analysis=fitting_trend>/scripts(.tex/.m)`；屬跨單一模型的方法分析，無 model topic 層；`gen_*_per_radius.m` 逐半徑產生 trend、`*_latex.m` 輸出 .tex 表格、`nofixl`=不固定 ℓ 變體。
**相關**：見 ../README.md、../../CLAUDE.md；對應擬合結果 ../charge_model_fitting/long2016_hexapole_halfcut/；電流對齊見 .claude/rules/fit-current-matches-sim.md、符號慣例見 .claude/rules/charge-model-source-convention.md。
