# long2016_hexapole_halfcut/fixl_fit — 固定-ℓ 點電荷 K_I 擬合功能組

**用途**：對 Long Fei 下極半切六極做**固定有效長度 ℓ（fixed-ℓ）的點電荷模型 K_I 擬合**。以 lsqnonlin 把 6 顆 coil 的 FEM 場（1A）配成模型 `B = gB·K̂_I·I_vec·kernel(p; ℓ, d̂)`，掃不同取樣半徑 R 找穩定區間，輸出擬合參數 {ℓ̂, ĝ_B, K̂_I} 與每半徑 LaTeX 結果表。這是**多腳本 + 繪圖 + LaTeX 結果**的完整功能組。

**內容**：
- `code/scripts/` — 擬合主力腳本：`sweep_KI_radius.m`（球取樣多半徑掃 fit）、`fit_KI_full.m`（單一全節點 fit）、`sweep_alln_vs_R.m`（NRMSE vs R）、`calib_fem.m`、`compare_KI_frob.m`、`KI_deviation_from_ideal.m`、`match_target_matrix.m`、`test_joint_6coil_fit_40um.m`。
- `code/plot/` — `plot_KI_convergence.m`、`plot_fixl_convergence.m`、`plot_fixl_params_vs_R.m`（收斂 / 參數 vs R 繪圖）。
- `figures/` — 輸出 PNG（`KI_cost_convergence_gB*.png`、`fixl_cost_convergence.png`、`sweep_alln_*_vs_R.png`）。
- `results/` — auto-gen LaTeX：`KI_trend_results.tex`、`testset_error_50_500.tex` + 子夾 `charge_fit/`、`per_radius/`。

**資料來源 / 流向**：讀 ANSYS_data `coilN/*.dat`（`import_ansys_data` + `filter_iron_nodes`，WP frame，source 符號 = negate FEM B）→ lsqnonlin 擬合 → `.mat` 寫 MATLAB_data（`charge_fit/...`，`matlab_path`）+ `.tex` 寫本組 `results/` + PNG 寫本組 `figures/`。

**命名 / 慣例**：功能組 schema = `code/{scripts,plot}/` + `figures/`（PNG）+ `results/`（auto-gen `.tex`）。⚠ 擬合電流必須 = FEM 激發電流（1A）；部分舊腳本仍寫 `I_actual=0.6`（同坑），改前須使用者確認，見 `.claude/rules/fit-current-matches-sim.md`。source 符號慣例見 `.claude/rules/charge-model-source-convention.md`。

**相關**：見 `../README.md`、`../../CLAUDE.md`；讀結果防呆見 `.claude/rules/result-read-safety.md`。
