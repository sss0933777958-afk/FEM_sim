# doc/charge_model_fitting/kuo_quadrupole/ — Kuo 四極點電荷/多極擬合

**用途**：Kuo 4-pole MEMS Quadrupole 的場擬合文件：面內磁場的 multipole-expansion 擬合模型，以及非正交磁極軸的基底投影法。
**內容**：`scripts/`（.tex 原稿）+ `pdf/`（編譯 PDF）。代表檔：`multipole_model`（四極面內 B 的多極展開擬合）、`nonorthogonal_basis_projection`（磁極軸非正交投影 λ=E\v 統一法）；PDF 含不同取樣範圍變體 `fitting_R0500`、`fitting_R0500_cube70`、`fitting_F20_cube70`、`fitting_R0500_sphere`。
**命名 / 慣例**：第二層 topic = 模型名 `kuo_quadrupole`；PDF 後綴標取樣設定（`R0500`=R 範圍、`cube70`/`sphere`=取樣區形、`F20`=變體 tag）。
**相關**：方法通法見 ../general/；投影法 memory `feedback_nonorthogonal_basis_projection`；見 ../README.md、../../../CLAUDE.md。
