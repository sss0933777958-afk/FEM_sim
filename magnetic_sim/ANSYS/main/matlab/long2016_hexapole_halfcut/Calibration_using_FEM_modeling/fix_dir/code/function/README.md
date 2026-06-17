# …/fix_l/code/function/ — fix-ℓ 模型輔助函式

**用途**：fix-ℓ 點電荷模型的數學/IO 輔助函式；`main.m` 全部從這裡呼叫。
**內容**：
- `load_coils.m` — 一次載入 6-coil FEM 場（讀 `.dat`）。
- `select_ball.m` — 取 WP 半徑 R 內的節點。
- `fit_KI_fixl.m` — `lsqnonlin` fit `{K̂, ℓ, gB}`（呼叫 `charge_residual`、`unpack_params`）。
- `charge_residual.m` — 殘差函式；`unpack_params.m` — 參數打包/解包。
- `region_field_err.m` — 區域相對 RMS 場誤差。
- `write_KI_tex.m` — 輸出純結果 `.tex`。

**資料來源 / 流向**：`load_coils` 讀 `ANSYS_data/.dat` → fit/誤差計算 → `write_KI_tex` 寫 `../../results/fix_l/*.tex`。

**命名 / 慣例**：純函式（一檔一函式）；I=1 A 對齊 FEM；電荷全 source。

**相關**：見上層 `../README.md`。
