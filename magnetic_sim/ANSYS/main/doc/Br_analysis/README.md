# doc/Br_analysis/ — Br（尖端表面磁場）分析推導

**用途**：放磁極尖端錐面上磁場分量分析的推導文件 —— 尖端 FEM 量到的 B_FEM 如何拆成法向 B_F 與切向 B_n、兩者之間的 cosine（夾角）關係，供 sensor 法向投影與磁通推導使用。
**內容**：LaTeX 原稿（.tex）+ 編譯 PDF，依 model topic 分群。代表主題：`bf_cosine`（B_F / B_n 向量圖 + 夾角 cosine 投影）。
**命名 / 慣例**：analysis-first schema：`<analysis=Br_analysis>/<topic=model 名>/{scripts(.tex), pdf}`；第二層 topic = 物理模型名（目前 `long2016_hexapole_halfcut`）。
**相關**：見 ../README.md（doc 結構）、../../CLAUDE.md；流程 SOP 見 ../workflows/（如 apdl-postproc / field-plot）。
