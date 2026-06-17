# long2016_hexapole_halfcut/validation — 電流組合閉環驗證功能組

**用途**：驗證 fixl_fit 擬合出的點電荷模型 {ℓ̂, ĝ_B, K̂_I}：用該模型對**多種電流組合**（單極 basis 與組合 combo）預測 WP 區場，跟 FEM 疊加場（superposition）比對 NRMSE，當作擬合品質的閉環驗證。含逐點誤差、per-coil / per-combo NRMSE 與疊圖。

**內容**：
- `code/scripts/` — `eval_basis_nrmse.m`（6 顆單極 basis NRMSE vs R）、`eval_combo_nrmse.m`（5 組電流組合 NRMSE vs 取樣半徑）、`eval_testset_error.m`、`eval_validate_combos.m`、`validate_KI.m`、`error_per_point_R100.m`、`per_coil_err_R100.m`。
- `code/plot/` — `plot_combo_nrmse_overlay.m`（組合 NRMSE 疊圖）。
- `figures/` — `basis_nrmse_*.png`、`combo_nrmse_C1..C5.png` + `_overlay.png`、`testset_error_vs_R.png`、`validate_combos_nrmse_R130/R150.png`。

**資料來源 / 流向**：讀 ANSYS_data `coilN/*.dat`（`import_ansys_data` + `filter_iron_nodes`，WP frame，source 符號 negate）做 FEM 組合場 + 讀 fixl_fit 擬合參數 → 算 NRMSE → PNG 寫本組 `figures/`（誤差表視需要寫 fixl_fit `results/`）。

**命名 / 慣例**：功能組 schema = `code/{scripts,plot}/` + `figures/`。NRMSE_max = sqrt(mean ‖B_model−B_FEM‖²)/max‖B_FEM‖×100。模型參數在 1A 單位，電流向量直接代入（不補正），見 `.claude/rules/fit-current-matches-sim.md`。

**相關**：見 `../README.md`、`../fixl_fit/README.md`（被驗證的擬合）、`../../CLAUDE.md`；讀結果防呆見 `.claude/rules/result-read-safety.md`。
