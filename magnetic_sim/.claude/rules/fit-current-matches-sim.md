# 擬合電流必須對齊模擬激發電流(強制讀取)

**原則**:做任何「把 FEM 場配進模型」的擬合時,**模型裡的電流 `I` 必須等於該 FEM 解的激發電流**。不可把「實際操作電流」塞進擬合的 `I` —— 那會讓增益參數(gB)默默吃進一個假比例,且日後用參數預測不同電流時必須補正,極易出錯。

當任何工作涉及:
- 寫 / 改任何 `*_fit*.m`、`sweep_KI_*.m`、`fit_KI_*.m`、charge-model 擬合腳本
- 改 `I_actual` / `I_in` / 模型電流變數
- 用 fit 出的參數**預測不同電流**下的場(例如多電流組合、actuation 掃描)

**動手前先讀完此規則**。

對應 memory:`feedback_fit_current_matches_sim.md`
相關 memory:[[long-fei-b-bar-matrix-v4]]、[[long2016-halfcut-KI-fit]]、[[charge-model-source-convention]]

---

## 🔒 核心規則

點電荷模型:`B(p) = gB · K̂ · I_vec · kernel(p; ℓ, d̂)`。

- **FEM 目前在 1A 解**(`CURR_ARRAY=1` → 70 匝 × 1A;raw `.dat` 即 1A 場)。
- 擬合時餵的是 1A 的場 → **模型 `I` 必須 = 1**。
- 如此 fit 出的 `{ℓ̂, ĝ_B, K̂}` 是「**每安培單位模型**」:任意電流向量 `I_vec`(單位安培,= FEM 激發單位)可**直接代入**,`B = ĝ_B·K̂·I_vec·kernel`,**不需任何補正因子**。

## ❌ 反例(踩過的坑,2026-06-13)

`sweep_KI_trend*.m` / `eval_testset_error.m` 原本寫 `I_actual = 0.6`(誤把元件操作電流當模型 I),卻去配 1A 的 FEM 場:
```
解的是:  gB · K̂ · 0.6 · kernel ≈ B_FEM(1A)
```
→ `1A/0.6 ≈ 1.67` 的假因子被 **ĝ_B 吃掉**(ĝ_B 被灌大 1.67×)。
→ 拿這組參數預測 1A/2A 時,若直接代入會**高估 1.67×**,每張誤差圖浮一個假底。
**已修正為 `I_actual = 1`**(2026-06-13)。

數值影響(供查核):`I` 0.6→1 只讓 `gB ×0.6`;`ℓ̂`、`‖K̂‖_F`、所有 err **不變**(這些跟整體增益無關)。

## 例外:何時可用操作電流(非此坑)

`gen_Vout_Vin*.m`、`compute_B_bar*`、bs-matrix 等用 `I_in = 0.6 A` 當**真實操作電流做合法縮放**(把 1A 的 FEM 場線性縮到 0.6A 操作點輸出電壓)—— 那是正確的,不是這個坑。**判準**:0.6 是「配 FEM 用的 I」(❌ 錯)還是「把 1A 場縮到操作點」(✅ 對)。

## ⚠ 尚待對齊(同坑、本次未改)

以下 fit 腳本仍寫 `I_actual = 0.6` 且配 1A FEM(同坑),但會動到 document 版 K_I / R_a 既載結果,**改前須使用者確認**:
- `fit_KI_full.m` / `fit_KI_full_zh.m`(document 版電荷擬合,影響 R_a)
- `sweep_KI_radius.m` / `sweep_KI_radius_zh.m`
- `test_joint_6coil_fit_40um.m`
- `plot_KI_convergence.m`

發現 / 處理後在此清單更新。

## 觸發片語(任一即啟動)

- 「改 I_actual / I_in」/「擬合電流」/「模型電流」
- 「寫新 fit / charge fit」/「預測不同電流 / 電流組合」
- 「為什麼 gB 這麼大 / 預測高估」

## 何時不適用

- 純後處理 / 畫圖、不涉及擬合電流。
- 明確用操作電流做縮放的 V/V 矩陣(見上「例外」)。
