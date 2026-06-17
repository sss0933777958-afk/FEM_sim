# doc/error_definition/ — 擬合誤差定義（NRMSE 正規化）

**用途**：定義擬合 / 場比對所用的誤差度量，特別是 NRMSE 的兩種正規化方式（除以最大場 vs 除以 RMS），供 charge fit、B 矩陣等各分析統一引用。
**內容**：`scripts/`（.tex 原稿）+ `pdf/`（編譯 PDF）。代表檔：`nrmse`（NRMSE_max 與 NRMSE_rms 兩式定義）。
**命名 / 慣例**：analysis-first schema：`<analysis=error_definition>/{scripts(.tex), pdf}`；屬跨模型通用定義，第二層直接是 scripts/pdf（無 model topic 層）。
**相關**：見 ../README.md、../../CLAUDE.md；擬合用法見 ../charge_model_fitting/、../fitting_trend/。
