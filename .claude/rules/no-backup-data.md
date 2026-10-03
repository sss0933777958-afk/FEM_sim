# 不使用 `backup/` 的資料與程式

1. 不 `addpath` 到 `magnetic_sim/ANSYS/backup/`，不從那裡讀 `.m`、`.mat`、`.dat`、`.db`。
2. 模型設定一律走該校正包的 `model_config(model, geom)`，不呼叫 `mt_constants()`：
   ```matlab
   cfg = model_config('long2016_hexapole_halfcut', 'tip40um');
   ```
3. 每個模型用自己的 config，不可互抄。
4. 腳本裡的根目錄變數命名為 `CALROOT`（`CAL` 常與既有變數撞名）。

例外：刻意重現歷史結果時可讀 `backup/`，但必須在回覆與註解中標明「這是 backup 歷史資料」。
