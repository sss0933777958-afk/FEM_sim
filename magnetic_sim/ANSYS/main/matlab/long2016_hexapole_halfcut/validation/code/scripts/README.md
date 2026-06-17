# long2016_hexapole_halfcut/validation/code/scripts — 閉環驗證腳本

**用途**：電流組合閉環驗證的運算腳本（非繪圖）：把擬合模型對單極 / 組合電流的預測場跟 FEM 疊加場比 NRMSE。

**內容**：`eval_basis_nrmse.m`（6 單極 basis）、`eval_combo_nrmse.m`（5 組合 vs R）、`eval_testset_error.m`、`eval_validate_combos.m`、`validate_KI.m`、`error_per_point_R100.m`、`per_coil_err_R100.m`。

**資料來源 / 流向**：讀 ANSYS_data `coilN/*.dat`（`import_ansys_data` + `filter_iron_nodes`，WP frame，negate 為 source）+ 讀 fixl_fit 擬合 `.mat` → 算 NRMSE / 逐點誤差 → 結果供 `../plot/` 畫圖、誤差表視需要寫 fixl_fit `results/`。

**命名 / 慣例**：`code/scripts/` 放運算；`eval_*` = 評估、`validate_*` = 驗證、`*_R100` = 固定半徑。模型 1A 單位、電流直接代入不補正，見 `.claude/rules/fit-current-matches-sim.md`。

**相關**：見 `../../README.md`、`../../../../CLAUDE.md`。
