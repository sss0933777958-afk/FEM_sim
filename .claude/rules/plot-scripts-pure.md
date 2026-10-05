# 繪圖腳本只畫圖

```
main/main.m + function/   校正 → 存 data/<model>/.mat
utils/scripts/<model>/    讀 .mat → 運算 → 存 utils/data/<model>/
plot/                     讀 .mat → 畫圖 → 存 figures/<model>/{current,voltage}/
```

`plot/` 底下的腳本（及任何 `plot_*.m`）只准：載入既有 `.mat`、取值、畫圖並輸出圖檔。

**不准**：解線性系統或擬合、讀 `.fld` / `.dat`、建取樣設計或跑收斂迴圈、產生新 `.mat`、呼叫 `model_config` 自己算幾何。

**判準**：把畫圖的行拿掉後還在算東西，那段就搬去 `utils/scripts/<model>/`，結果存到 `utils/data/<model>/`。

**可以留在繪圖腳本**：座標旋轉、單位換算、直方圖分箱、在既有資料上畫參考線這類呈現用的輕量運算。
