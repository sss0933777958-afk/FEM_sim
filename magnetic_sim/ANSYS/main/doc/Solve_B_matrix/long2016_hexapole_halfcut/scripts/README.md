# doc/Solve_B_matrix/long2016_hexapole_halfcut/scripts/ — B_S 矩陣 LaTeX 原稿

**用途**：存放本 topic B_S / V_out/V_in 矩陣推導的 **LaTeX 原稿（.tex）**，含 MATLAB auto-gen 的矩陣 body 片段。
**內容**：`Bs_derivation.tex`（B_S 代數推導）、`Bbar_S_4p572_body.tex`（baseline 面積平均矩陣）、`Bbar_S_4p572_gap200um_body.tex`（μ_r 等效氣隙變體）、`Vout_Vin_4p572_body.tex`（量測端電壓比矩陣）。
**命名 / 慣例**：葉夾 = .tex；`*_body.tex` 是 `\input` 進主文件的矩陣表片段（由 gen_*_latex.m 產生）；後綴 `4p572` = sensor 位置 mm、`gap200um` = 氣隙變體。
**相關**：PDF 產物見 ../pdf/；產生器在 matlab/long2016_hexapole_halfcut/.../results；見 ../README.md。
