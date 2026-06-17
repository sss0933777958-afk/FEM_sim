# doc/charge_model_fitting/ — 點電荷模型擬合（推導 + 各模型結果）

**用途**：把 FEM 磁場配進 Long Fei 6-charge 點電荷模型 `B(p)=g_B·K̂·I_vec·kernel(p;ℓ,d̂)` 的所有文件：通用推導、方法詳解（[A]/[J]/[B-6x]）、初始猜測、k11 推導，以及各物理模型上的擬合結果與誤差。
**內容**：依用途分四層子夾：`fitting_derivation/`（共通推導 / reference）、`general/`（跨模型方法通法 .md）、`<model>/`（各模型結果：`kuo_quadrupole` / `long2016_hexapole_full` / `long2016_hexapole_halfcut`）。檔型 = .tex 原稿 / .pdf 編譯 / .m（LaTeX 表格產生器）/ .docx reference。
**命名 / 慣例**：analysis-first schema：`<analysis=charge_model_fitting>/<topic>/{scripts(.tex/.m), pdf, reference}`；第二層 topic = 物理模型名（特例 `fitting_derivation` 通用推導、`general` 跨模型方法）。
**相關**：見 ../README.md、../../CLAUDE.md；SOP 見 ../workflows/charge-model-fit.md（強制先跑 validity sweep，err<5% 才採信）；電流對齊見 .claude/rules/fit-current-matches-sim.md、符號慣例見 .claude/rules/charge-model-source-convention.md。
