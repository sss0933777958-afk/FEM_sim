# Calibration 結果一律加印 gauge 前物理轉移矩陣 ᴮĤ（強制）

**使用者拍板（2026-07-20；2026-07-21 符號 bar→hat）**：`Calibration_using_FEM_modeling/` 的 **current_base 與 voltage_base**，結果 PDF **除了 gauge 後的無因次矩陣（K̄_I / D̄），一律再加印 gauge 前的完整物理轉移矩陣 ᴮĤ**：

| 專案 | 無因次（既有） | **新增：gauge 前物理矩陣** | 單位 |
|---|---|---|---|
| current_base | K̄_I（K̄(1,1)=5/6） | **ᴮĤ_I = ĝ_I · K̄_I** | **mT/A** |
| voltage_base | D̄（D̄(1,1)=5/6） | **ᴮĤ_V = ĝ_V · D̄** | **mT/mV** |

> **符號**：物理轉移矩陣一律寫 **hat（ᴮĤ）**，不是 bar（ᴮH̄）——與 solve 函式回傳的 `H_I`/`H_V`（code 註解即稱 Ĥ）一致。PDF LaTeX 標籤 = `{}^{B}\hat{H}_{I}` / `{}^{B}\hat{H}_{V}`。（K̄_I / D̄ / ĝ 保留各自 bar/hat accent，只有 H 是 hat。）

當工作涉及：
- 改 / 重寫 `current_base` 或 `voltage_base` 的 `emit_model_results.m`
- 新增 Calibration 變體 / 新 pipeline 的結果輸出
- 討論 K̄_I / D̄ / ĝ / ᴮĤ 的呈現

**動手前先讀完此規則。**

對應 memory：`feedback_calibration_transfer_matrix_output.md`
相關規則：`unit-reference.md`（ᴮĝ_I=mT/A、ᴮĝ_V=mT/mV）、`results-pdf-only.md`（results/ 只放 PDF）、`charge-model-source-convention.md`（號誌）、`pole-coil-numbering.md`（map / K̄_I 判準）。

---

## 🔒 核心定義（別記錯）

ᴮĤ = **gauge 前的完整物理轉移矩陣**，就是 solve 函式回傳、之前「算了但沒印」的那顆：
- current：`solve_KI_bar_gain` 回 `H_I = G·Fᵀ(FFᵀ)⁻¹`；ĝ_I=(6/5)H_I(1,1)、K̄_I=(5/(6·g11))H_I ⇒ **ᴮĤ_I = ĝ_I·K̄_I = H_I**。
- voltage：`solve_D_bar_gain` 回 `H_V = G·Vᵀ(VVᵀ)⁻¹`；ĝ_V=(6/5)H_V(1,1)、D̄=(5/(6·h11))H_V ⇒ **ᴮĤ_V = ĝ_V·D̄ = H_V**。

即 **K̄_I / D̄ 是 ᴮĤ 除掉 (1,1) 元素做 gauge 的無因次版**；ᴮĤ 把物理增益乘回去（帶單位）。三者關係：**ᴮĤ = ĝ · (K̄ or D̄)**，數值上 = solve 回傳的 H。

## 🔒 實作要點

1. **用「顯示版」的 gauge 矩陣乘 ĝ**，確保號誌一致（對角全正、all-source）：
   - current_base emit 的 `K_bar` 已含 `coil_sign=[1 -1 1 -1 -1 1]` 顯示翻號 → `ᴮĤ_I = ghat_I_B .* ? ` 寫成 `ghat_I_B * K_bar`（用翻號後的 K_bar）。
   - voltage_base 的 `D_bar` 已是 all-source（Bstack flip-sink、不再翻號）→ `ᴮĤ_V = ghat_V_B * D_bar`。
2. **位置**：印在對應的 K̄_I / D̄ 矩陣**正下方**（同一顆矩陣的兩種形式擺一起）。
3. **標籤 / 單位**（照 `unit-reference.md` + `figure-style.md`）：
   - `{}^{B}\hat{H}_{I}~[\mathrm{mT/A}]`、`{}^{B}\hat{H}_{V}~[\mathrm{mT/mV}]`（**hat，不是 bar**）。
   - 數值 auto-factor 照既有 `emit_tex` T.mat（|指數|≥2 才抽 ×10^n）——ᴮĤ_V ~10⁻³ 會抽 ×10⁻³、ᴮĤ_I ~10⁰ 不抽。
4. **不動 .mat / solve 函式**：H_I/H_V 早已算出（voltage 存為 `Dmat`、current 由 emit 現算 `ĝ·K̄`），只在 emit 端多印一顆矩陣，不改上游。

## 驗證
- ᴮĤ 對角必**全正**（= ĝ>0 × gauge 矩陣對角正）；**ᴮĤ(1,1) = h11 = (5/6)·ĝ**（因 gauge 矩陣(1,1)=5/6、ĝ=(6/5)h11 ⇒ ĝ·(5/6)=h11=ᴮĤ(1,1)）——不是 ĝ 本身。
- 兩邊 `model_results_<variant>.pdf` 應同時出現「無因次 K̄_I/D̄」+「物理 ᴮĤ」兩顆；標籤是 **hat**（`\hat{H}`）。

---

# 🔒 附則：bias 偏移一律「物理 e [µm]」+「無因次 e/ℓ̂」兩張都印（強制）

**使用者拍板（2026-07-31）**：`USE_BIAS=true`（eighteen）的結果 PDF，除了既有的物理偏移 `e [µm]`，
**一律再加印無因次偏移 `e/ℓ̂`**（＝求解器內部單位、`Pc = Pc_base + e` 直接用的那個 e）。

**理由**：`Pc = Pc_base + e/ℓ̂` 是以 **ℓ̂ 為長度單位**的無因次格點；PDF 只印 µm 版，讀者要自己除 ℓ̂ 才能重建 Pc、
算 `|Pc|`（電荷到 WP 的距離）或 `q/r²`（自場貢獻）——**印出來省掉這步、也避免除錯 / 乘錯方向**
（PDF 的 µm 是「乘過」ℓ̂，不是除過）。

## 實作要點
1. **兩張都印、上下相鄰**：先 `e [µm]`（物理），緊接 `e/\hat{\ell}`（無因次）。
2. **呈現形式不同**（使用者拍板 2026-07-31）：
   - `e [µm]` → **`T.e`**（帶 `P1..P6` 欄標 + `e_x/e_y/e_z` 列標的表格），給人讀。
   - **`e/ℓ̂` → `T.mat`**（**標準 `bmatrix`、無任何標籤**），照 `figure-style.md` 數值標註慣例 #4「矩陣用標準
     bmatrix、不用欄位標籤表格」——這張是要直接搬進算式 / 程式的矩陣。
   - ⚠ `emit_tex.m` 的 `emit_mat` 已改成**任意尺寸**（`[nr,nc]=size(Ms)`，原寫死 6×6）才能吃 3×6 的 E36；
     既有 6×6 呼叫輸出**逐字不變**。
3. **標籤**：`e~[\mu\mathrm{m}]` / `e/\hat{\ell}`。**無因次那張不標單位**（照 `figure-style.md`「無單位不標」）。
4. **不加 auto-factor**：兩張 fmode 都留空 `''`，維持 `%9.4f` 原值——`e/ℓ̂` 量級 ~10⁻¹~10⁻³，
   抽因子反而難跟上面的 µm 版逐格對照。
5. **不動 .mat / solve**：`rec.e` 存的本來就是無因次 e，µm 版是 emit 端乘 `ℓ̂·1e6` 得到；只在 emit 端多印一張。
6. 實作位置：`Calibration_using_FEM_modeling/function/emit_results.m` 的 `if rec.USE_BIAS` 區塊（current / voltage 共用）。

## 驗證
- 兩張的每個元素必滿足 **`e[µm] = (e/ℓ̂) × ℓ̂[µm]`**（ℓ̂ 印在同頁）。
  例（long2016 R150 eighteen current，ℓ̂=857.3 µm）：P1 `e_x` 19.6024 µm ÷ 857.3 = **0.0229** ✓
- `USE_BIAS=false`（single）**不印**這兩張（e ≡ 0）。

## 觸發片語
- 「改 calibration 結果輸出 / emit_results」
- 「加印 / 輸出 ᴮĤ / H_I / H_V / 物理轉移矩陣」
- 「加印 e/ℓ̂ / 無因次偏移 / 電荷格 Pc」
- 新增 Calibration variant 要出結果 PDF 時

## 何時不適用
- 非 Calibration 的其他分析 PDF（bs_matrix 等）——各自慣例。
- 純幾何 / mesh / 場圖。
