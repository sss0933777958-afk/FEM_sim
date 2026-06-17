# long2016_hexapole_halfcut/fixl_fit/code/scripts — K_I 擬合腳本

**用途**：固定-ℓ 點電荷模型 K_I 擬合的運算腳本（非繪圖）。

**內容**：`sweep_KI_radius.m`（球取樣 R=50:50:500µm 逐半徑 lsqnonlin fit）、`fit_KI_full.m`（單一全節點 fit）、`sweep_alln_vs_R.m`（NRMSE/參數 vs R）、`calib_fem.m`、`compare_KI_frob.m`（Frobenius 比較）、`KI_deviation_from_ideal.m`、`match_target_matrix.m`、`test_joint_6coil_fit_40um.m`。

**資料來源 / 流向**：讀 ANSYS_data `coilN/*.dat`（`import_ansys_data` + `filter_iron_nodes`，WP frame，source 符號 negate）→ lsqnonlin 擬合 {K̂_I, ℓ, gB} → `.mat` 寫 MATLAB_data（`charge_fit/...`）+ `.tex` 寫 `../../results/`（含 `charge_fit/`、`per_radius/`）。

**命名 / 慣例**：`code/scripts/` 放運算；`sweep_*` = 掃參數、`fit_*` = 單次擬合、`compare_*` = 比較。⚠ 擬合電流對齊 FEM（1A）；部分腳本仍 `I_actual=0.6`，改前確認，見 `.claude/rules/fit-current-matches-sim.md`。

**相關**：見 `../../README.md`、`../../../../CLAUDE.md`。
