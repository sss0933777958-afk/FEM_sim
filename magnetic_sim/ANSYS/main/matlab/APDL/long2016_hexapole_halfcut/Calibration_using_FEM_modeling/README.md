# matlab/long2016_hexapole_halfcut/Calibration_using_FEM_modeling/ — FEM 點電荷校正（current_base / voltage_base）

**用途**：用 FEM 場校正 Long Fei 半切六極的「點電荷模型」與其 Hall-sensor 延伸。
**2026-07-15 重構**：原 4 子夾（`fix_dir`/`no_fix_dir`/`Hall_sensor_base_fix_dir`/`Hall_sensor_base_no_fix_dir`）
合成 **2 個專案**，各用**單一主程式 + `USE_BIAS` 參數**統一 fix / no_fix（差別只在有沒有加 bias）。

## 兩個專案
- **`current_base/`** — 電荷/電流模型（fit `{K̄, ℓ̂, ^Bĝ_I}`）。`code/main/main.m` 的 `USE_BIAS`：
  - `false` = **fix**（電荷在磁極軸 `ℓ̂·Pc_base`，ê=0，只擬合 ℓ̂）。
  - `true`  = **no_fix**（18-param bias：電荷離軸 `ℓ̂·(Pc_base+E(ê))`，擬合 ℓ̂+ê(17)）。
  - R150/gap_calibrate：BIAS ℓ̂≈**856 µm**, ^Bĝ_I≈**6.53 mT/A**；FIX(varpro) ℓ̂≈**868 µm**, ^Bĝ_I≈**5.75 mT/A**（region err 3.2%）。
- **`voltage_base/`** — Hall-sensor/電壓模型（輸出 `V, D̄, ^Bĝ_V, Ĥ_V`）。同 `USE_BIAS`；
  **依賴 current_base**（`load_coils_actuator/select_ball/build_A/make_Pc/fit_varpro`）+ 自身 sensor 程式
  （`build_sensor_geometry/extract_Vmat_interp`）。R150/gap_calibrate：BIAS ^Bĝ_V≈**7.14e-3 mT/mV**。

## 統一擬合 = 變數投影（Varpar / variable projection）
擬合函式 **`fitting`**（2026-07-16 合併原 `fit_varpro`+`bias_resid`+`make_Pc`）：逐激發把 6 個線性增益 `g_j`
用最小平方 **profile 掉**，優化器只看非線性 `θ={ℓ̂ (+ê(17) 當 USE_BIAS)}`，回傳 `(ℓ̄, ê, Pc, J)`。
`USE_BIAS=false` ⇒ ê 凍結 0 ⇒ 1-D 擬合 ℓ̂（= 舊 fix）；`USE_BIAS=true` ⇒ 18-param bias（= 舊 no_fix）。

## 內層架構（current_base 與 voltage_base 一致）
```
<project>/
├── code/
│   ├── main/          main.m（USE_BIAS；擬合→解矩陣→存 .mat→末端呼叫 emit 產 PDF）
│   ├── main_function/ 唯一函式夾：loader/fitting/solve/sensor + 產 PDF 的 emit_model_results(/gen_B_matrix)
│   └── plot/          繪圖/診斷（含搬自舊 function/ 的診斷 helper）
├── data/     (.mat 擬合成果)
├── figures/{single_param, eighteen_param, shared}/  (.png；依模型分：fix / 18-param bias / 比較·中性)
└── results/{single_param, eighteen_param}/          (.pdf；依 USE_BIAS 落對子夾，main 末端 emit 產)
```

## PDF/LaTeX（2026-07-16 更新：main 末端產 PDF）
- **`main.m` 最後一步就呼叫 `emit_model_results`（voltage 另加 `gen_B_matrix`）產 PDF** → `results/`（不再需手動另跑）。
- emit 函式仍在 **`code/main_function/`**、仍可單獨執行（自載 `data/*.mat`，寫 `.tex` + inline xelatex → `results/*.pdf`）。
- （2026-07-16 重整：`function/` 併入 `main_function/`；`build_A→build_S_matrix`、`fit_varpro+bias_resid+make_Pc→fitting`、
  `gauge_KI→solve_KI_bar_gain`、`extract_Vmat_interp→build_V_matrix`、新增 `solve_D_bar_gain`；刪 `fit_bias`/`emit_charge_Q` 等死檔。）

**資料來源**：讀 `ANSYS_data/long2016_hexapole_halfcut/data/<variant>/coilN/` 的 6-coil FEM 場（`.dat`，1 A）→
取 WP 半徑 R 內真實節點 → variable-projection 擬合 → 存 `.mat`。

**慣例**：單一主程式組（`code/main/main.m` + config）；model-first；電荷模型全 source（flip-sink 只翻下極）；
I_actual = 1 A 對齊 FEM 激發；`results/` 只放 `.pdf`。

**相關**：見上層 `../README.md`、`../../../CLAUDE.md`、`.claude/rules/{actuator-frame,charge-model-source-convention,fit-current-matches-sim,results-pdf-only}.md`。
