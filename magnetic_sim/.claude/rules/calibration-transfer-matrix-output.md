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

# 🔒 附則：bias 偏移**只印無因次 e/ℓ̂ 一張**（強制）

**使用者拍板（2026-08-04，取代 2026-07-31 的「兩張都印」）**：`USE_BIAS=true`（eighteen）的結果 PDF，
**e 矩陣只印無因次的 `e/ℓ̂`（＝已除以 ℓ̂ 的單位向量版）**，**不再印物理 `e [µm]`**。

**理由**：校正解出來的 `e` 本來就是**以 ℓ̂ 為長度單位**的無因次量（`rec.e` 存的即是），
`Pc = Pc_base + e/ℓ̂` 也直接吃這個值；µm 版是 emit 端額外乘 `ℓ̂` 生出來的衍生量，
既非求解器內部表示、又多一個乘錯方向的機會 → 拿掉，只留原生的無因次版。
（要 µm 值的人自己乘同頁印出的 ℓ̂ 即可。）

> 📌 **歷史**：2026-07-31 曾拍板「µm + 無因次兩張都印」（更早還一度用帶 `P1..P6` 欄標的 `T.e` 表格印 µm 版）。
> **兩者皆已作廢**，現行 = 只印 `e/ℓ̂` 一張標準 bmatrix。

## 實作要點
1. **只印一張**：`e/\hat{\ell}`，用 **`T.mat`**（標準 `bmatrix`、無任何欄/列標籤）。
   **不印** `e~[\mu\mathrm{m}]`；**不用** `T.e`（帶 `P1..P6` 欄標 + `e_x/e_y/e_z` 列標的表格）。
   照 `figure-style.md` 數值標註慣例 #4「矩陣用標準 bmatrix、不用欄位標籤表格」——它是要直接搬進算式 / 程式的矩陣。
   - ⚠ `emit_tex.m` 的 `emit_mat` 已是**任意尺寸**（`[nr,nc]=size(Ms)`，原寫死 6×6）才能吃 3×6 的 E36；
     既有 6×6 呼叫輸出**逐字不變**。**APDL 與 Maxwell 兩個分支都要改**（各有一份 `emit_tex.m`）。
2. **標籤**：`e/\hat{\ell}`。**不標單位**（無因次；照 `figure-style.md`「無單位不標」）。
3. **不加 auto-factor**：fmode 留空 `''`，維持 `%9.4f` 原值（`e/ℓ̂` 量級 ~10⁻¹~10⁻³，抽因子反而難讀）。
4. **不動 .mat / solve**：`rec.e` 存的本來就是無因次 e —— 這次改動只是**停止**在 emit 端乘 `ℓ̂·1e6` 那張。
5. 實作位置：**兩個分支各一份**——`matlab/APDL/Calibration_using_FEM_modeling/function/emit_results.m` 與
   `matlab/Maxwell/function/emit_results.m` 的 `if rec.USE_BIAS` 區塊（current / voltage 共用）。
   **改動要兩邊同步。**

## 驗證
- PDF 內**只有一張 e 矩陣**，標籤 `e/\hat{\ell}`、`\begin{bmatrix}`（3×6）、**無單位**。
- **不該再出現** `e~[\mu\mathrm{m}]` 那張，也不該出現 `P1..P6` 欄標的 e 表格。
- 數值可對回 `rec.e`（無因次）；要 µm 就乘同頁的 ℓ̂
  （例 long2016 R150 eighteen current，ℓ̂=857.3 µm：`e_x(P1)=0.0229` × 857.3 = 19.6 µm）。
- `USE_BIAS=false`（single）**不印** e（e ≡ 0）。

## 觸發片語
- 「改 calibration 結果輸出 / emit_results」
- 「加印 / 輸出 ᴮĤ / H_I / H_V / 物理轉移矩陣」
- 「加印 e/ℓ̂ / 無因次偏移 / 電荷格 Pc」
- 新增 Calibration variant 要出結果 PDF 時

## 何時不適用
- 非 Calibration 的其他分析 PDF（bs_matrix 等）——各自慣例。
- 純幾何 / mesh / 場圖。
