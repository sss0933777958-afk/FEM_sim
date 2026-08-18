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

---

# 🔒 附則三：擬合誤差一律印 **NMAE（向量範數版）**，不印 RMSPE（強制）

**使用者拍板（2026-08-15 改用 NMAE；2026-08-17 改為向量範數版）**。

| | 定義 | 出現位置 |
|---|---|---|
| ~~RMSPE~~（**停用**） | `sqrt(Σεᵢ² / Σbᵢ²)·100` | 只留 `.mat`（`rec.RMSPE`），**PDF 不印** |
| ~~NMAE-L1~~（**已作廢**） | `Σ\|εᵢ\| / Σ\|bᵢ\| · 100`（逐**分量**絕對值） | 只留 `.mat`（`rec.NMAE_L1`），**PDF 不印** |
| **NMAE**（現行） | **`[ Σ_j Σ_i ‖b_ij − S_i·ᴮĝ·M̄·u_j‖ / N_p ] / b̄ · 100`** | current 與 voltage 的 PDF 都印這個 |

其中 `M̄·u_j` = current 的 `K̄_I·F_j` / voltage 的 `D̄·V_j`；`b̄ = Σ_j Σ_i ‖b_ij‖ / N_p`（**與分子同結構**，
故 `N_p` 對消、比值無因次）。

- **關鍵差異**：先把每個「點 i × 激發 j」的 **3 維殘差取歐氏長度 ‖ε_ij‖** 再相加，
  不是逐分量取絕對值。同一份殘差下 **新值必 < 舊 L1 值**（實測比值穩定在 **0.73~0.76**）。
- `S_i·ĝ·M̄·u_j ≡ S_i·G_j`（因 `Ĥ·u = G`），故實作直接用 `S*G − Bstack` 當殘差，只改**聚合方式**：
  ```matlab
  e_ij = sqrt(sum(reshape(resid,  3, []).^2, 1));   % 每點每激發的 ‖ε‖
  b_ij = sqrt(sum(reshape(Bstack, 3, []).^2, 1));
  rm.NMAE = sum(e_ij) / sum(b_ij) * 100;
  ```
  ⚠ `reshape(·,3,[])` 正確的前提是 Bstack 以 `[bx;by;bz]` 逐點堆疊（現行 `main.m` 即如此）。
  ⚠ 唯一例外：`h11 < 0` 的退化區（ĝ 已取 abs）`ĝ·M̄·u = −G`，該區間本就判定不可用。
- 實作位置（**四支 solve + 兩支 emit，兩個分支都要同步**）：
  - `matlab/{Maxwell,APDL/Calibration_using_FEM_modeling}/function/solve_current.m`、`solve_voltage.m`
    → 算 `rm.NMAE` 與 `rm.NMAE_L1`（`rm.RMSPE` 保留不刪）
  - 兩支 `main/main.m` → `rec.NMAE = rm.NMAE`
  - 兩支 `function/emit_results.m` → LaTeX（current 版，voltage 把 `\hat{g}_{I}\bar{K}_{I}F_j` 換成 `\hat{g}_{V}\bar{D}V_j`）：
    `\[ \mathrm{NMAE} = \dfrac{\sum_j\sum_i \left\| b_{ij} - S_i\,{}^{B}\hat{g}_{I}\bar{K}_{I}F_j \right\| / N_p}{\bar{b}}\cdot 100 = %.2f\% \]`
- **三代指標的數字互不可並列比較**（RMSPE=L2 / NMAE-L1=分量 L1 / NMAE=向量範數）。要比就一起重跑。

## 換指標時的回填法（已驗證可重複，2026-08-15 建立、2026-08-17 再次沿用）

`.mat` 沒存 `P`/`Bstack`，但存了設定。做法＝**從 `.mat` 的
`model/GEOM/VARIANT/DATASET/R_select/GRID_NRPT` 重建取樣點，再用存著的 `l_hat/e/Fmap/V`
呼叫原 `solve_*` 取殘差** —— 只重算誤差、不重跑擬合。

- **驗證判準＝「RMSPE 與 `.mat` 逐位相同」**（2026-08-17 實測 12/12 皆 `0.0e+00`）。對不上就跳過不寫入。
- `rec.VARIANT` 存的是**輸出名**（含 `_convN…/_k22_…/_soff…`），要
  `regexprep(rec.VARIANT,'_convN\d+|_k22_[0-9p]+|_soff[0-9p]+mm','')` 剝掉才是資料 variant。
  ⚠ `maxwell_mesh0p06` 是**真的資料 variant**，不可剝。
- 回填後把 `rec.NMAE` 寫回 `.mat` 再 `emit_results`，否則 `.mat` 與 PDF 不同調。

## ⚠ 連帶：sensor 距非定案值時輸出檔名要帶 tag

`main.m` 原本只在 **`.mat`** 檔名加 `_soff<值>mm`，**PDF 沒有** → 換 sensor 位置重跑會覆蓋掉
4.572 mm 版的 PDF。2026-08-15 已改成把 soff tag 併進 `VAR_OUT`（PDF 檔名也帶），
`.mat` 則由 `VAR_OUT` 一次帶出，不再重複加後綴。

---

## 觸發片語
- 「改 calibration 結果輸出 / emit_results」
- 「RMSPE / NMAE / 擬合誤差指標」
- 「加印 / 輸出 ᴮĤ / H_I / H_V / 物理轉移矩陣」
- 「加印 e/ℓ̂ / 無因次偏移 / 電荷格 Pc」
- 新增 Calibration variant 要出結果 PDF 時

## 何時不適用
- 非 Calibration 的其他分析 PDF（bs_matrix 等）——各自慣例。
- 純幾何 / mesh / 場圖。
