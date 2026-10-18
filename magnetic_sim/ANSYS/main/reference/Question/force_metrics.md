# 問題：力模型能否定義類似磁通模型的 Efficiency index 與 isotropy（待討論）

狀態：**問題已提出、尚未分析**（2026-10-05）。

## 1. 磁通模型的現行定義

符號依 `reference/charge_model_fitting/fitting_derivation/pdf/Calibration using FEM modeling_V2.pdf`。
位置 p_i 處，電流到磁場的轉移是**線性**的：

b_ij = S(p_i/ℓ̂, ê) · ᴮĝ_I K̄_I · I_j

- p̄_i = p_i/ℓ̂，p̄_c = [p̄_c1 ⋯ p̄_c6]。
- S(p_i/ℓ̂, ê) = [s_i1 ⋯ s_i6]（3×6），s_ik = (p̄_i − p̄_ck) / ‖p̄_i − p̄_ck‖³。
- ᴮĤ_I = ᴮĝ_I · K̄_I：6×6 電流轉移矩陣 [mT/A]。
- 對 S(p_i/ℓ̂, ê)·ᴮĤ_I（3×6）做 SVD，得奇異值 σ₁ ≥ σ₂ ≥ σ₃（𝒞、κ 沿用 `solve_current.m` 與論文圖的記號）：
  - **Efficiency index** 𝒞 = σ₁σ₂σ₃（圖上用 𝒞^(1/3)，單位回到 mT/A）
  - **isotropy** κ = σ₃ / σ₁（0~1，1 為各方向等效）
- 實作：`matlab/Flux/Maxwell/function/hexapole/solve_current.m` 的 `control_metrics`（逐節點 svd，取平均 / 最差值）；
  直方圖 `figures/paper_fig_plot/plot/plot_gain_iso_hist.m`。

## 2. 力模型的形式

符號依 `reference/Force_model_fitting/Modeling_Calibration - Magnetic Force_CHM`（與 flux 重疊的符號以 flux 文件為準）：

m(p) = ᴹĝ_B · b(p)，f(p) = ½ ∇(m·b) = ½ ᴹĝ_B ∇(b·b)

p̄ = p/ℓ̂，r̄_n = p̄ − p̄_cn，r̄_n = ‖r̄_n‖

L(p̄; n, m, ⋮) = l_n,m(p̄) = (1/(r̄_n³ r̄_m³)) [ (1 − 3(r̄_n·r̄_m)/r̄_n²) r̄_n + (1 − 3(r̄_n·r̄_m)/r̄_m²) r̄_m ]

L(p/ℓ̂, ê) 是 6×6×3 無因次空間分布張量，L_k（k = 1, 2, 3 即 x, y, z）是它的第 k 個 6×6 對稱切片。

電流版：f_ijk = (ᴹĝ_B / 2ℓ̂) (ᴮĝ_I)² · I_jᵀ K̄_Iᵀ L_k(p_i/ℓ̂, ê) K̄_I I_j = I_jᵀ ᶠĤ_Iᵀ L_k(p_i/ℓ̂, ê) ᶠĤ_I I_j

電壓版：f_ijk = ᶠĝ_V · v_jᵀ M̄ᵀ L_k(p_i/ℓ̂, ê) M̄ v_j = v_jᵀ ᶠĤ_Vᵀ L_k(p_i/ℓ̂, ê) ᶠĤ_V v_j

即電流（電壓）到力是**二次型**（不可疊加），不是線性映射；每個力分量對應一個 6×6 對稱矩陣 ᶠĤ_Iᵀ L_k ᶠĤ_I。

## 3. 問題

1. 磁通模型的 Efficiency index 與 isotropy 的**概念與計算方式**能否直接套用到力模型？
   （S(p_i/ℓ̂, ê)·ᴮĤ_I 是線性映射才能用 SVD；力是電流的二次型，直接對 ᶠĤ_I 或某個 3×6 矩陣做 SVD 是否有物理意義？）
2. 若不能直接套用，有什麼數學方法能從力模型算出「這個 tweezer 產生力的效率」與「各方向產生力的均勻度」，
   作為與磁通 𝒞、κ 對應的指標？

## 4. 可考慮的方向（未驗證，僅供討論起點；文件未定義新記號，以下用文字描述）

- **工作點線性化**：在某一組偏置電流 I_j 附近，f_ij 對電流的 Jacobian 是 3×6，第 k 列 = 2 I_jᵀ ᶠĤ_Iᵀ L_k(p_i/ℓ̂, ê) ᶠĤ_I；
  對它做 SVD 得 𝒞、κ。結果依賴 I_j 的選擇。
- **單位電流下的最大可達力**：對每個力方向（單位向量，分量為權重），把三個切片 L_k 依權重加總成一個 6×6 對稱矩陣，
  在 ‖I_j‖ = 1 下該方向可達的最大力 = ᶠĤ_Iᵀ(加權後的 L) ᶠĤ_I 的最大特徵值；
  效率取各方向最大可達力的最小值（或可達力集合的體積），等向性取 min/max 比。
- 其他：可達力集合（力的 workspace / force ellipsoid）的幾何量、與電流限制（‖I_j‖∞ 或 ‖I_j‖₂）的關係。

## 5. 待決

- 指標要對應「電流範數」哪一種限制（‖I_j‖₂ 或各線圈 ‖I_j‖∞）。
- 是否需要和磁通的 𝒞、κ 在同一組節點、同一校正結果上比較。
