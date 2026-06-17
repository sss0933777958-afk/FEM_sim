# doc/charge_model_fitting/long2016_hexapole_halfcut/ — 下極半切 hexapole 點電荷擬合

**用途**：Long Fei 2016 六極**下極半切** hexapole（V4 主軸）的點電荷擬合文件：擬合誤差定義，以及各組擬合結果 LaTeX 表格的 MATLAB 產生器（K̂_I、球殼 / cube40 取樣、R040 等）。
**內容**：`scripts/`（.tex 原稿 + .m LaTeX 產生器）。代表檔：`fit_error_def.tex`（誤差定義 ΔB=B_fit−B_FEM，仿 dissertation Fig 2.6c）、`gen_KI_latex.m`、`gen_KI_ball_latex.m`、`gen_KI_R040_latex.m`、`gen_fit_J_cube40_latex.m`（各取樣範圍的 K̂_I 表格產生器）。
**命名 / 慣例**：第二層 topic = 模型名 `long2016_hexapole_halfcut`；本 topic 只有 `scripts/`（無獨立 pdf/ 葉夾）；`.m` 產生器讀 MATLAB_data 的 fit 結果輸出 .tex 表格，後綴標取樣（`ball`/`R040`/`cube40`）。
**相關**：方法通法見 ../general/；誤差正規化見 ../../error_definition/；R 掃描趨勢見 ../../fitting_trend/；見 ../README.md、../../../CLAUDE.md。
