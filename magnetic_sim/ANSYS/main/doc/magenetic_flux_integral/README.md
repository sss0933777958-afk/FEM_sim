# doc/magenetic_flux_integral/ — 磁通積分（網格數 / 格點值收斂）

**用途**：沿軸磁通 Φ=∫B·dA 數值積分的網格相關推導：網格數量與網格間距如何取（N=81×81…）、格點取值方式對積分結果的影響與收斂。
**內容**：`scripts/`（.tex 原稿）+ `pdf/`（編譯 PDF）。代表檔：`mesh_size`（網格數量與間距）、`grid_value`（格點取值）。
**命名 / 慣例**：analysis-first schema：`<analysis=magenetic_flux_integral>/{scripts(.tex), pdf}`；**注意原資料夾名拼字為 `magenetic`（少一個 t，沿用既有命名，勿改）**；屬跨模型通用方法，無 model topic 層。
**相關**：見 ../README.md、../../CLAUDE.md；後處理抽場 SOP 見 ../workflows/apdl-postproc.md；磁通相關 memory `project_long_fei_p1p2_flux_profile`。
