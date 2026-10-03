# 校正包結構凍結

適用四個校正包：`matlab/{Flux,Force}/APDL/Calibration_using_FEM_modeling/`、`matlab/{Flux,Force}/Maxwell/`。

## 結構

```
config/<model>/[<geom>/]   模型設定與幾何（R_act、Pc_base、sensor、方法旗標）
function/                  校正求解函式
main/main.m                driver（頂部旗標切換 model / variant / base / USE_BIAS）
utils/scripts/              額外計算的腳本：讀校正 .mat → 運算 → 存到 utils/data/
utils/data/                 額外計算的結果 .mat
plot/                       只做繪圖的腳本
common_path/                路徑 resolver
data/<model>/               校正結果 .mat，只留最終定案的版本（自描述：結果 + 設定條件）
results/<model>/{single,eighteen}/   PDF（single = USE_BIAS false、eighteen = true）
figures/<model>/{current,voltage}/   圖檔（不再分 single / eighteen）
```

資料流：`main + function`（校正 → `data/`）→ `utils/scripts/`（運算 → `utils/data/`）→ `plot/`（→ `figures/`）。

## 規則

1. **新增 / 改名 / 移動 / 刪除資料夾或檔案、重組檔案切分、新增校正模型，一律先問。**
2. 改既有檔的內文（修 bug、加參數或分支）不受限。
3. 模型設定一律由 `config/` 提供，`function/` 只消費，不在函式裡寫死幾何。Force 沒有自己的 `config/`，經 `load_flux_calib` 使用 `Flux/Maxwell/config/`。
4. 與校正無關的分析（從校正結果拿資料、額外計算、畫圖看）一律走 `utils/scripts/` → `utils/data/` → `plot/`，不另開暫存資料夾。
5. `data/` 只留最終定案的校正結果；被取代的舊版本刪除。
