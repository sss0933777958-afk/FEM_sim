# 規則 #3：`current_base/` 與 `voltage_base/` 的 `results/` 只放 PDF

**使用者拍板（2026-06-26；2026-07-15 隨 Calibration 重構更新路徑）**：
`matlab/<pkg>/Calibration_using_FEM_modeling/{current_base,voltage_base}/results/`
**只留最終 `.pdf`**（矩陣 / 電荷 Q / d 向量等 xelatex 排版輸出）。
（原規則綁 `Hall_sensor_base_fix_dir/results/`；4 子夾已於 2026-07-15 合成 `current_base`/`voltage_base`。）

**2026-07-15 依模型分子夾**：`results/` 底下再分 **`single_param/`**（USE_BIAS=false、fix）與
**`eighteen_param/`**（USE_BIAS=true、18-param bias）——`emit_*.m` 依 `USE_BIAS`/`fit_{fixl,bias}` 自動落對子夾。
`results/` root 只留子夾 + README（no-mixed-files）；每子夾只放 `.pdf`。（`figures/` 另分
`single_param`/`eighteen_param`/`shared`，見 `figure-output.md`。）

## 🔒 檔名（Maxwell 分支，2026-08-17 使用者拍板縮短）

`matlab/Maxwell/results/<model>/{single,eighteen}/` 的 PDF 一律：

```
<base>_R<半徑µm>[_N<點數>][_<tag>][_soff<N>mm].pdf
```

- `<base>` = `current` | `voltage`
- `_R150` **一律標**（不再只在 ≠150 時標）
- `_N<點數>` 只在用等測度網格降取樣（`GRID_NRPT` 非空）時出現；全格點版沒有這段
- `_<tag>` = `rec.VARIANT` 去掉分支名 `maxwell` 與 `convN<n>` 後的殘餘（如 `mesh0p06`）
- 例：`current_R150.pdf`、`current_R150_N528.pdf`、`voltage_R150_N80_soff3mm.pdf`

模型名與 single/eighteen 由**資料夾**表達、Maxwell 分支由整棵樹表達 → 不再放進檔名
（舊格式 `model_results_current_maxwell_convN528.pdf` 38 字元 → 新 17 字元）。
實作在 `matlab/Maxwell/function/emit_results.m` 的 `stem` 組法。
⚠ **只改了 Maxwell 分支**；APDL 分支維持舊命名 `model_results_<base>_<variant>…`。

## 🔒 規則
- 只留 `*.pdf`；**`.mat / .tex / .aux / .log / .txt` 一律不留**。
  - `.mat` 分析成果 → 改放同專案 `data/`（見 `matlab-output-layout.md`）。
  - `.tex / .aux / .log` 排版中間檔 → 產 PDF 的 `code/main_function/emit_*.m`（`emit_model_results`/`gen_B_matrix`）xelatex 編完即刪。
- **（2026-07-16 更新）`main.m` 的最後一步就是產 PDF**：main 擬合 → 解矩陣 → 存 `.mat` → **呼叫 `emit_model_results`（voltage 另加 `gen_B_matrix`）寫 `.tex` + inline xelatex → `results/*.pdf`**。emit_* 仍是獨立函式（可單獨跑），但 main 會在末端自動呼叫、不再需要手動另跑。（舊規：main 只擬合不碰 PDF，已作廢。）
- ~~同步該專案 `results/README.md`~~ —— **[2026-08-17] `matlab/Maxwell/results/README.md` 已依使用者
  指示廢除刪除**（內容過期：舊檔名格式 + 舊 NMAE 定義 + 已作廢的 K̄(2,2) 揭露）。檔名規則改記在本檔上方。

## 觸發片語
- 動到 `current_base/results/` 或 `voltage_base/results/` 內容時
- 「results 只放 pdf」/「清 results 中間檔」

## 何時不適用
- 其他功能組的 `results/`（如 `bs_matrix/results`）目前仍放 `.tex` 排版原稿——本規則只綁 `current_base`/`voltage_base` 的 `results/`。

相關：`matlab-output-layout.md`、memory `reference_local_latex_compile`、`feedback_matlab_local_data_layout`。
