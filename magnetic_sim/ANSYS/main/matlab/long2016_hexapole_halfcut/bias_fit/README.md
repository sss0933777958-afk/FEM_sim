# long2016_hexapole_halfcut/bias_fit/ — no_fix_l 18-param bias 電荷模型功能組

**用途**：Long Fei 半切六極（halfcut hexapole）的 **no_fix_l「18 參數 bias 電荷模型」**腳本集。在 actuator 框下對每個取樣半徑 R 做 bias fit（理想格點 `Pc = Pc_base + E(e)`，`ℓ̂` + per-pole `e_hat` joint，g_j profile 掉後 gauge 出 `gB·K̄_I`），掃 R 趨勢、跟舊全-K̂ 模型比 NRMSE、檢查電荷是否落在極內，並產出收斂/趨勢圖與每半徑 LaTeX 表。

**內容**（功能組 root，下分標準 `code/{scripts,plot}` + `figures/` + `results/`）：
- `code/scripts/` — fit / 掃描 / 比較 / 診斷的計算腳本（`sweep_nofixl_vs_R`、`calib_fem_bias`、`compare_models_nrmse`、`check_charge_in_pole`）。
- `code/plot/` — 收斂圖、參數-R 趨勢圖、電荷側視動畫腳本（畫圖前須讀繪圖規則）。
- `figures/` — 已定案 `.png` / `.gif`（收斂、參數 vs R、P1/P2 電荷動畫幀）。
- `results/per_radius/` — auto-gen 的每半徑 `fit_KI_*.tex`（R050…R500）。

**資料來源 / 流向**：腳本經 `common/ansys_path('long2016_hexapole_halfcut','coilN','standard')` 讀 6 顆 coil 的 FEM `.dat`（`'all'` dataset，全 in-ball 節點）→ 旋進 actuator 框做 bias fit → 成果 `.mat` 寫 `matlab_path(model,'charge_fit',...)`（`MATLAB_data`）→ 圖進 `figures/` 與 `main/figures/long2016_hexapole_halfcut/` → 每半徑 `.tex` 進 `results/per_radius/`。
- 模型電流 `I_actual = 1`（對齊 FEM 1A 激發，per `fit-current-matches-sim` 規則）。
- 顯示符號全 source（翻上極 P2/P4/P5），per charge-model-source-convention。

**命名 / 慣例**：功能組 schema = `code/{scripts,plot}` + `figures/` + `results/`；計算腳本放 `scripts/`、繪圖腳本放 `plot/`；`results/` 只放 auto-gen `.tex`（不手改）。APDL coil j → paper pole 索引 `[1,3,6,5,2,4]`。

**相關**：見上層 `../README.md`（matlab/ 功能組 schema）、`../common/README.md`（resolver）、`../../../CLAUDE.md`。
