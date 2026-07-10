# …/fix_l/code/function/ — fix-ℓ 非 main 輔助函式

**用途**：fix-ℓ 的非 main 輔助（SVD/legacy tex）。**main 校正流程函式（pipeline + 輸出 `emit_model_results`）已移到 `../main_function/`**（load_coils/select_ball/fit_KI_fixl/charge_residual/unpack_params/region_field_err/calc_range_metrics/emit_model_results）。
**內容**：
- （`write_KI_tex.m` 已於 2026-07-10 刪＝legacy；結果改由 `../main_function/emit_model_results.m` 產）
- `svd_SiHI_WP.m` —（獨立腳本）在 WP 中心組電流→場轉移 `T = S_i·Ĥ_I`（`Ĥ_I=gB·Khat`）並 SVD，回報 `U/Σ/Wᵀ`；WP 中心 `S_i=−dhat`（與 ℓ̂ 無關）。存 `../../data/svd_SiHI_WP_<variant>.mat`。

**資料來源 / 流向**：`load_coils` 讀 `ANSYS_data/.dat` → fit/誤差計算 → `../main_function/emit_model_results.m` 產 results PDF。

**命名 / 慣例**：純函式（一檔一函式）；I=1 A 對齊 FEM；電荷全 source。

**相關**：見上層 `../README.md`。
