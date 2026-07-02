# …/fix_l/code/function/ — fix-ℓ 模型輔助函式

**用途**：fix-ℓ 點電荷模型的數學/IO 輔助函式；`main.m` 全部從這裡呼叫。
**內容**：
- `load_coils.m` — 一次載入 6-coil FEM 場（讀 `.dat`）。
- `select_ball.m` — 取 WP 半徑 R 內的節點。
- `fit_KI_fixl.m` — `lsqnonlin` fit `{K̂, ℓ, gB}`（呼叫 `charge_residual`、`unpack_params`）。
- `charge_residual.m` — 殘差函式；`unpack_params.m` — 參數打包/解包。
- `region_field_err.m` — 區域相對 RMS 場誤差。
- `write_KI_tex.m` — 輸出純結果 `.tex`。
- `svd_SiHI_WP.m` —（獨立腳本）在 WP 中心組電流→場轉移 `T = S_i·Ĥ_I`（`Ĥ_I=gB·Khat`）並 SVD，回報 `U/Σ/Wᵀ`；WP 中心 `S_i=−dhat`（與 ℓ̂ 無關）。存 `../../data/svd_SiHI_WP_<variant>.mat`。
- `calc_range_metrics.m` — `m=calc_range_metrics(P,Ĥ_I,ell_m,dhat)`：對 R≤150µm 球內真實節點逐點 SVD → 回 `.sigma_tot`(mean ‖T‖_F)/`.iso_tot`(mean σmax/σmin)/`.sigma_min`/`.iso_worst`/`.Np`。`main.m` 呼叫、存進 fit mat、console 印。
- `emit_model_results.m` —（獨立）載 fit mat 出結果 PDF（K̄_I/ℓ̂/G/F/^Bĝ_I + **σ_tot/iso_tot** 控制範圍代表值）到 `../../results/`。

**資料來源 / 流向**：`load_coils` 讀 `ANSYS_data/.dat` → fit/誤差計算 → `write_KI_tex` 寫 `../../results/fix_l/*.tex`。

**命名 / 慣例**：純函式（一檔一函式）；I=1 A 對齊 FEM；電荷全 source。

**相關**：見上層 `../README.md`。
