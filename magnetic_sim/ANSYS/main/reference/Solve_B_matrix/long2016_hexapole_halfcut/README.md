# doc/Solve_B_matrix/long2016_hexapole_halfcut/ — 下極半切 hexapole 的 B_S 矩陣

**用途**：Long Fei 2016 六極**下極半切** hexapole 的 sensor 轉移矩陣推導：B_S / B̄_S（sensor disc 面積平均）、V_out/V_in、以及 D_H 六極轉移矩陣。
**內容**：`scripts/`（.tex 原稿，含 auto-gen 矩陣 body）+ `pdf/`（編譯 PDF）。代表檔：`Bs_derivation`、`Bbar_S_4p572_body`、`Bbar_S_4p572_gap200um_body`（μ_r 等效氣隙變體）、`Vout_Vin_4p572_body`、`Bs_matrix` / `D_matrix_hexapole` / `dh_derivation`。
**命名 / 慣例**：第二層 topic = 模型名；矩陣 body `.tex` 由 MATLAB（matlab/<model>/.../results）產生；sensor 變體後綴如 `4p572`（sensor 位置 mm）、`gap200um`（氣隙）。
**相關**：見 ../README.md、../../../CLAUDE.md；SOP 見 ../../workflows/bs-matrix-derive.md。
