# doc/Solve_B_matrix/ — B_S 轉移矩陣推導

**用途**：推導 sensor 輸出轉移矩陣 —— 6 顆極的自激 FEM 場如何映到 6 顆 Hall sensor 的 B_S / B̄_S，以及量測端 V_out/V_in 與 D_H（6×6 DC 轉移矩陣）的代數推導。
**內容**：LaTeX 原稿（.tex，含 auto-gen 的矩陣 body）+ 編譯 PDF，依 model topic 分群。代表主題：`Bs_derivation`（B_S 推導）、`Bbar_S`（含 baseline 與 gap200um 變體）、`Vout_Vin`、`D_matrix_hexapole` / `dh_derivation`（D_H 矩陣）。
**命名 / 慣例**：analysis-first schema：`<analysis=Solve_B_matrix>/<topic=model 名>/{scripts(.tex), pdf}`；第二層 topic = 物理模型名（目前 `long2016_hexapole_halfcut`）。
**相關**：見 ../README.md、../../CLAUDE.md；流程 SOP 見 ../workflows/bs-matrix-derive.md（強制先驗 coil1 資料夾陷阱）。
