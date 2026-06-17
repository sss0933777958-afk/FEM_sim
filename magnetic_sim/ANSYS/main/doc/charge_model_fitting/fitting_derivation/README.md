# doc/charge_model_fitting/fitting_derivation/ — 點電荷擬合通用推導

**用途**：點電荷模型擬合的**共通推導與校正流程**（跨模型、與特定 model 結果無關）：FEM 校正流程（Calibration using FEM modeling）、初始猜測值、k11 簡併推導、no_fix_l（不固定 ℓ）變體，以及來源 reference 文件。
**內容**：`scripts/`（.tex 原稿）、`pdf/`（編譯 PDF + 來源報告 PDF）、`reference/`（外部來源 .docx）。代表檔：`initial_guess_cjk`（K̂ 初值對角 5/6、非對角 −1/6）、`k11_derivation`、`Calibration using FEM modeling_V2.pdf`、`no_fix_l.pdf`、`Magnetic hexapole model`。
**命名 / 慣例**：`fitting_derivation` 是 charge_model_fitting 下的通用推導 topic（非 model 名）；其下分 scripts/pdf/reference 三葉夾。
**相關**：方法通法見 ../general/；各模型結果見 ../long2016_hexapole_*、../kuo_quadrupole；見 ../README.md、../../../CLAUDE.md。
