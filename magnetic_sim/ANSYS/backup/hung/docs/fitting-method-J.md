# Method [J]：Joint 6-Coil 磁荷擬合（Hung）

## 概述

同時使用 6 組 single-coil FEM 數據，聯合擬合 6 個磁荷位置。
每組 coil 激勵讓不同的 pole 成為主導，提供足夠的約束使位置唯一確定。

**結果：平均向量誤差 = 0.46%，位置 spread ≈ 0（唯一解）**

## 模型

Coil k 通 1A 時，WP 區域內位置 p 的磁場：

```
B_k(p) = C_k × k_m × Σ_i [ w_ki × (p - r_qi) / |p - r_qi|³ ]
```

- `k_m = mu_0 / (4*pi) = 1e-7`（磁 Coulomb 常數）
- `r_qi`：第 i 個磁荷位置（6 coils 共用同一組）
- `w_ki`：Coil k 激勵時，磁荷 i 的權重
- `C_k`：Coil k 的振幅因子（per-coil，解析求解），`C_k = N_c × I / (mu_0 × R_a_k)`

注意：Hung 使用正號慣例（無 -w），C_k > 0，所有 tip 為 SOURCE。

## 權重向量

```
w_k = K_I × I_vec_k
```

- `I_vec_k`：電流向量（第 k 個位置 = 1，其餘 = 0）
- `K_I = I_6 - ones(6)/6`：理想磁通分配矩陣

### K_I 的物理意義

K_I 描述磁通如何在 6 極之間分配。當 Coil k 激勵時：

```
w_k = K_I × [0,...,1,...,0]' = [−1/6, ..., 5/6, ..., −1/6]'
                                        ↑ 第 k 個
```

- 激勵極權重 = 5/6（主導磁荷）
- 其餘 5 極各 = −1/6（回流磁通均分）
- 每欄總和 = 0（磁通守恆）

K_I 來自 Hopkinson's law 的磁路分析（Long 2016, Eq. 2.4），假設回流磁通在 5 極間均分。

## 最佳化

### 參數

- **搜尋（fminsearch）**：18 個位置參數（6 poles × 3 絕對座標）
- **解析求解（VarPro）**：6 個 per-coil C_k

總有效參數：24 個

### 初始值

使用 per-coil single-charge ell sweep 的結果作為起始點：

1. 對每個 Coil k，只用 1 個磁荷，sweep ell 找最佳距離
2. 結果：Lower ell = 954 um, Upper ell = 908 um
3. 初始位置：`r_qi_init = ell_layer × d_hat_i`

d_hat_i 是原點到 tip i 的單位向量（magic angle 54.74° from Z）。

### Variable Projection（VarPro）

搜尋變數 x 是一個 18 × 1 的向量，包含 6 個磁荷的 3D 絕對座標：

```
x = [x_1, y_1, z_1, x_2, y_2, z_2, ..., x_6, y_6, z_6]'
```

reshape 成 3 × 6 矩陣後，每欄就是一個磁荷的位置 r_qi。
fminsearch 每次嘗試不同的 x，VarPro 對應算出最佳 C_k：

```
C_k = (b_unit_k' × b_fem_k) / (b_unit_k' × b_unit_k)
```

其中 b_unit_k 是 C=1 時的模型場。這把 optimizer 從 24 維降到 18 維。

### Cost Function

```
cost(x) = Σ_k ||C_k × b_unit_k - b_fem_k||²
```

對 6 個 coil 求殘差平方和。fminsearch（Nelder-Mead）最小化此值。

### 設定

- 擬合區域：+-100 um cube，隨機取樣 5000 節點
- TolX = 1e-10, TolFun = 1e-22, MaxFunEvals = 50000
- 2 組不同初始條件（per-layer ell, +50um noise）

## 結果

### 磁荷位置

```
Pole  類型     x [um]   y [um]   z [um]   |c| [um]
P1    Lower     ~797      ~1     ~-473      919.8
P2    Upper    ~-686      ~1     ~+514      916.5
P3    Lower    ~-410    ~+708    ~-454      919.8
P4    Upper    ~+343    ~-594    ~+518      916.5
P5    Upper    ~+344    ~+597    ~+514      916.5
P6    Lower    ~-407    ~-701    ~-460      919.8

Lower 平均 |c|: 919.8 um
Upper 平均 |c|: 916.5 um
```

### Per-Coil C_k 與 R_a

```
Coil    類型     C_k          R_a [A/Wb]
Coil1   Lower    7.171e-02    ~7.8e+08
Coil2   Upper    6.313e-02    ~8.8e+08
Coil3   Lower    7.168e-02    ~7.8e+08
Coil4   Upper    6.295e-02    ~8.8e+08
Coil5   Upper    6.314e-02    ~8.8e+08
Coil6   Lower    7.174e-02    ~7.8e+08
```

C_k 在同層之間高度一致（Lower ~0.072, Upper ~0.063），反映上下層磁路的差異。

### 擬合誤差

```
Coil    向量誤差 [%]
Coil1       0.37
Coil2       0.53
Coil3       0.39
Coil4       0.56
Coil5       0.55
Coil6       0.38
平均:       0.46
```

### 唯一性（Identifiability）

2 組 trial 收斂到幾乎相同的位置（cost 差異 < 1e-8）。

## 與 Long [J] 的比較

| 指標              | Hung [J]     | Long [J]     |
|-------------------|--------------|--------------|
| 平均誤差          | 0.46%        | 1.11%        |
| Lower 平均 |c|    | 920 um       | 816 um       |
| Upper 平均 |c|    | 917 um       | 767 um       |
| K_I               | 理想         | 理想         |

## 前置步驟

1. 3-pass NREFINE mesh 加密（WP +-50 um 有 12,321 節點）
2. 6 組 single-coil FEM 模擬，正確輸出 .rmg（-m 8000）
3. Per-coil ell sweep（`fit_single_charge.m`）→ 初始位置

## 檔案

- 擬合腳本：`magnetic_sim/ANSYS/backup/hung/analysis/fit/fit_J.m`
- 出圖腳本：`magnetic_sim/ANSYS/backup/hung/analysis/plot/plot_J_quiver.m`（% 誤差）、`magnetic_sim/ANSYS/backup/hung/analysis/plot/plot_J_rmse.m`（RMSE）
- 依賴：
  - `magnetic_sim/ANSYS/backup/hung/analysis/core/mt_constants.m` — 幾何常數
  - `magnetic_sim/ANSYS/backup/hung/analysis/core/import_ansys_data.m` — FEM 數據讀取
- FEM 數據：`magnetic_sim/ANSYS/backup/hung/results/coil[1-6]/filleted/coilN_{coord,bfield}_wp.dat`

## Hung 特有注意事項

1. **所有 tip 都是 SOURCE**（B 從 tip 射出）— 與 Long 相反（Long: lower=SINK, upper=SOURCE）。
   原因：SOURC36 的 N1/N2 swap。使用正號慣例，C_k 為正值。

2. **apdl_to_paper = 恆等映射** [1,2,3,4,5,6] — Long 使用 [1,3,6,5,2,4]。

3. **理想 K_I 即可使用** — 曾經以為需要擬合 K_I（因 Coil2-6 error >100%），
   但實際原因是 FEM 結果檔（.rmg）缺失。修正後理想 K_I 給出 0.46% error。
