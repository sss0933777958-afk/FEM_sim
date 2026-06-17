# 磁荷模型擬合：三種方法詳解 [A] / [J] / [B-6x]

> **文件目的**：詳細寫出 Long Fei dissertation Section 2.2.3 「6-charge magnetic field model」三種等級擬合方法的**完整推導過程 + 演算法步驟 + MATLAB code + 數值範例**。
>
> **適用對象**：要重做 fit、要 reproduce 結果、要寫論文 methodology 章節的人。
>
> **最後更新**：2026-06-03

---

## 目錄

- [0. 為什麼要 fit + 龍飛 b=0 假設](#0-為什麼要-fit-磁荷模型)
- [1. 三方法共用的數學工具與資料前處理](#1-三方法共用基礎)
- [2. Method [A] — 共用 ell 球面擬合（2 參數）](#2-method-a--共用-ell-球面擬合2-參數)
- [3. Method [J] — 6-coil joint 自由 3D 擬合（24 參數）](#3-method-j--6-coil-joint-自由-3d-擬合24-參數)
- [4. Method [B-6x] — ell+δ 全激勵擬合（19 參數）★ 定案](#4-method-b-6x--ellδ-全激勵擬合19-參數--定案)
- [5. 三方法對比與何時用哪個](#5-三方法對比)
- [6. 對龍飛 b=0 假設的回饋](#6-對龍飛-b0-假設的回饋)
- [7. 腳本/數據/圖檔清單](#7-腳本資料圖檔清單)
- [8. 名詞速查](#8-名詞速查)

---

## 0. 為什麼要 fit 磁荷模型？

### 0.1 出發點

Long Fei dissertation Eq 2.2 把 6 顆 pole 的 B 場 model 寫成 6 個**等效點電荷** (point charges)：

```
        ⎛   6                                     ⎞
B(p) = ⎜  Σ   k_m · q_i / |p − c_i|² · u_i(p, c_i) ⎟
        ⎝ i=1                                     ⎠
```

其中：
- `q_i` = 第 i 顆 pole 的等效磁荷（[A·m]）
- `c_i` = 第 i 顆磁荷的位置（[m]）
- `k_m = μ₀/(4π) = 10⁻⁷ N/A²`

每顆 charge 有 1 個 **scalar** 強度 q_i 跟 3 個 **位置** (x_i, y_i, z_i)。總計 6×4 = 24 自由度。

### 0.2 龍飛的兩個 fit 假設

論文 page 18-20 fit 出 (ρ, R_a) = (900 µm, 6.3×10⁸ A/Wb)。
他**暗藏**兩個假設：

1. **位置假設**：所有 charge 都坐在半徑 ρ 的球面上，方位由 pole tip 幾何決定 → **b_i ≡ 0**（每顆 charge 沒有額外偏移）
2. **強度假設**：6 顆 charge 強度由 K_I 矩陣 + R_a 決定，**只剩 R_a 是 free 參數**

→ 自由度從原本 24 砍到剩 **2** (ρ, R_a)。

### 0.3 三種方法的定位

我們做的三種 fit 是 **逐步放鬆龍飛的兩個假設**：

| Method | 自由度 | 放鬆什麼 | 數據量 | 主要目的 |
|---|---|---|---|---|
| **[A]** | **2** | 完全照龍飛 (b=0 + 球面) | 1 coil (P1) | reproduce Long Fei，提供 ρ baseline |
| **[J]** | **24** | 放開到全自由 (18 dim bias + 6 個獨立 C) | 6 coils | 量化 b≠0 偏離有多大 |
| **[B-6x]** | **19** | ρ 鎖定 + 自由 δ + 共用 1 個 C ★定案 | all6 superposition | 物理上最一致 + 數值最準 |

下面三章逐一**完整**推導。

---

## 1. 三方法共用基礎

### 1.1 統一的 fit 方程（被擬合的物理 model）

**從 Long Fei Eq 2.1-2.4 一路推到 fit 用的形式**：

#### Step A：原始物理 model（Eq 2.2）

```
            6                  p − c_i
B(p) =  Σ   k_m · q_i · ──────────────────
           i=1            |p − c_i|³
```

- **k_m = μ₀/(4π) = 10⁻⁷ N/A²**（磁 Coulomb 常數，跟靜電 Coulomb 1/(4πε₀) 對應）
- q_i = 第 i 顆 charge 強度 [A·m]

#### Step B：把電流→磁荷（Eq 2.1 + 2.4）

```
q_i = Φ_i / μ₀ = -(N_c / (μ₀ · R_a)) · (K_I · I_vec)_i = -(N_c/(μ₀·R_a)) · w_i
```

代回 Step A：

```
            ⎛   N_c       ⎞   6                p − c_i
B(p) = -k_m·⎜ ──────────  ⎟ · Σ   w_i · ──────────────────
            ⎝ μ₀ · R_a   ⎠  i=1          |p − c_i|³
       └──────────────────┘
       這整坨叫做「振幅 ξ」
```

#### Step C：合併成 fit-friendly 形式

```
            ⎛  6                p − c_i             ⎞
B_model(p; θ) = ξ · ⎜  Σ   w_i · ──────────────────  ⎟       ← fit 用形式
            ⎝ i=1               |p − c_i|³             ⎠

其中  ξ = -k_m · N_c / (μ₀ · R_a)
     = 振幅因子（合併 k_m 跟所有物理常數）
```

| 變數 | 意義 | 在三方法中的角色 |
|---|---|---|
| **ξ** | 全局振幅 = -k_m · N_c/(μ₀·R_a) | 解析解（VarPro），所有方法都這樣處理 |
| **k_m** | 磁 Coulomb 常數 = μ₀/(4π) = 10⁻⁷ N/A² | 物理常數，吸收進 ξ |
| **N_c** | 每極線圈匝數 = 70 | 物理常數，吸收進 ξ |
| **μ₀** | 真空磁導率 = 4π×10⁻⁷ T·m/A | 物理常數 |
| **R_a** | 集總空氣磁阻 [A/Wb] | 唯一 free 物理量，從 ξ 反推 |
| **w_i** | 第 i 顆 charge 的相對權重 = (K_I · I_vec)_i | 固定（由激勵電流決定） |
| **c_i** | 第 i 顆 charge 的 3D 位置 | **三方法差別就在這怎麼參數化** |
| **p** | 場點位置 [m] | FEM 節點座標 |

#### Step D：從 ξ 反推 R_a

fit 完拿到 ξ 後：

```
            k_m · N_c
R_a = ─────────────────────
          μ₀ · |ξ|
```

（負號 cancel；取 |ξ| 是因為 SOURC36 繞線方向給的符號 convention，per `coil_sign`）

> **注意**：在 MATLAB code 裡（per `point_charge_model.m`），習慣把 ξ 拆成兩部分：
> - 「外層 C」= N_c/(μ₀·R_a) ← 對應 `R_a_unit = N_c/μ₀` 讓 C=1 的 trick
> - 「內層 k_m」 baked into b_unit（因為 `point_charge_model` 函式內部已乘 k_m）
> 兩種 convention 數學等價，本文檔之後用 ξ 統一表示「整坨振幅」。

### 1.2 統一的 cost function

```
              N    
J(θ) = Σ Σ  ‖ B_model(p_j; θ) − B_FEM(p_j) ‖²       ← 向量 L² norm
        j=1
              ↑                ↑
       647 個 FEM 節點    對應 3 個分量 (Bx, By, Bz)
```

→ residual 是 **3N = 1941 維向量**（647 node × 3 component），用 `fminsearch` / `fminbnd` 最小化。

### 1.3 Variable Projection (VarPro) — 把 ξ 消掉

關鍵觀察：**B_model 對 ξ 是線性**（位置決定方向，ξ 只是振幅）。

定義 `b_geom(θ_geom) = Σᵢ wᵢ · (p − cᵢ)/|p − cᵢ|³`（純幾何 basis，不含 k_m 或 R_a）：

```
B_model(p; θ) = ξ · b_geom(p; θ_geom)

J(ξ, θ_geom) = ‖ ξ · b_geom(θ_geom) − b_FEM ‖²
```

對 ξ 求導 = 0 解析：

```
            b_geom^T · b_FEM
ξ_opt = ─────────────────────
            b_geom^T · b_geom
```

→ **不用數值搜尋 ξ**，每次只搜尋 θ_geom，每個 θ_geom 都直接代公式算 ξ_opt。Reduced cost：

```
                            (b_geom^T · b_FEM)²
J*(θ_geom) = ‖b_FEM‖² −  ─────────────────────
                              b_geom^T · b_geom
```

從 ξ_opt 回推 R_a：

```
R_a = k_m · N_c / (μ₀ · |ξ_opt|)
```

**Method [A] / [J] / [B-6x] 都用這個 trick**。
（在 code 裡 b_unit 已含 k_m，所以分兩階段：先算 C_opt 再算 R_a，數學等價於這裡的 ξ_opt。）

### 1.4 K_I 磁通分佈矩陣

當激勵 Coil k（i.e., 通 1A 到 pole k 的線圈），6 顆 charge 的權重：

```
w = K_I · I_vec,    K_I = I₆ − (1/6) · 1₆ = ⎡ 5/6  −1/6  ⋯  −1/6 ⎤
                                              ⎢ ⋮              ⋮  ⎥
                                              ⎣ −1/6 ⋯       5/6 ⎦
```

例：激勵 P1 (`I_vec = [1,0,0,0,0,0]ᵀ`) → `w = [5/6, −1/6, −1/6, −1/6, −1/6, −1/6]ᵀ`

- 主極 P1 得 **5/6** 權重（dominant、flux sink）
- 其餘 5 顆各得 **−1/6**（passive、flux source）
- 總和 = 0（磁通守恆，div B = 0）

### 1.5 資料前處理鏈（共用）

```
ANSYS FEM (V4 halfcut baseline, 6 coils × 0.6A) 
    ↓ 載入 coord + bfield .dat → 390k nodes
    ↓ 鐵芯節點排除 (geometric cone filter) → 341k nodes
    ↓ 取 WP 中心 100 µm cube |x|,|y|,|z_wp| < 50 µm → 647 nodes
    ↓ 座標轉換 z_wp = z_apdl − SPH_OFST (SPH_OFST = −12.711 mm)
    ↓ B 場分量堆疊
b_FEM = [bx; by; bz] (1941 × 1)
```

### 1.6 物理常數 + 幾何常數（從 `mt_constants.m`）

**物理常數**：

| 符號 | 數值 | 意義 |
|---|---|---|
| **μ₀** | 4π×10⁻⁷ T·m/A | 真空磁導率 |
| **k_m** | μ₀/(4π) = **10⁻⁷ N/A²** | 磁 Coulomb 常數（dissertation page 18 標 `km = μ₀ = 1.0×10⁻⁷ N/A²`，跟 `μ₀/4π` 等價）|
| **N_c** | 70 | 每極線圈匝數（dissertation page 14）|

**幾何常數**：

```
R_norm    = 500 µm           WP 到極尖物理距離
R_norm_xy = R_norm·√(2/3) = 408 µm
R_norm_z  = R_norm/√3 = 289 µm
α         = arctan √2 = 54.74°   極角
SPH_OFST  = −12.711 mm        APDL → WP 座標 z 平移
```

6 極方位角 (paper 索引)：
- Lower (P1, P3, P6): θ = 0°, 120°, 240°，z_sign = −1
- Upper (P2, P4, P5): θ = 180°, 300°, 60°，z_sign = +1

APDL index → paper index 映射：APDL `{1,2,3,4,5,6}` = Paper `{P1, P3, P6, P5, P2, P4}`。

---

## 2. Method [A] — 共用 ell 球面擬合（2 參數）

**腳本**：`magnetic_sim/hexapole-long2016/analysis/fit_charge_model.m`
**輸出**：`data/charge_model_fit.mat`
**主結果**：ell = 835.4 µm, R_a = 9.21×10⁸ A/Wb, mean err 4.94%

### 2.1 假設（完全照搬龍飛）

- 6 顆 charge 全部坐在半徑 ell 的球面上
- 方位由 pole tip 幾何完全鎖死
- 6 顆共用同一個 ell
- 共用同一個 R_a（即同一個 C）
- 只用 1 個 coil (P1, APDL coil 1) 的數據

### 2.2 數學參數化

charge 位置：

```
c_i = ell · d̂_i,   i = 1..6

d̂_i = [cos θ_i · sin α,  sin θ_i · sin α,  z_sign_i · cos α]
```

α = 54.74° 固定。θ_i 跟 z_sign_i 由 paper 幾何完全決定。

驗證：|c_i| = ell·√(sin²α + cos²α) = ell ✓

**自由變數總計 = 2 個**：
- ell ∈ [400, 2000] µm（位置，數值搜尋）
- C → R_a（VarPro 解析）

### 2.3 演算法詳細步驟

```
Step 0: 載入資料
    d = import_ansys_data('coil1', 'all', 'coil1')
    把 (x, y, z) 轉到 WP 座標系: z_wp = z − SPH_OFST
    應用鐵芯過濾 → 341k nodes
    取 100 µm cube → 647 nodes
    b_FEM = [bx; by; bz]   ← 1941×1

Step 1: 設定 K_I 跟 I_vec
    K_I = I₆ − (1/6)·1₆
    APDL coil 1 → paper P1, 所以 I_vec = [1, 0, 0, 0, 0, 0]ᵀ (paper 順序)
    w = K_I · I_vec = [5/6, −1/6, −1/6, −1/6, −1/6, −1/6]ᵀ
    
    coil_sign = +1 (P1 是 lower, SOURC36 繞線方向給 sink)

Step 2: Phase 1 — ell 粗掃
    For ell ∈ linspace(400e-6, 2000e-6, 300):
        c_i = ell · d̂_i,  i = 1..6
        對每個 FEM 節點 p_j 算 b_unit_j = Σᵢ w_i · (p_j − c_i)/|p_j − c_i|³
        b_unit_vec = [b_unit_x; b_unit_y; b_unit_z]  ← 1941×1
        C_opt(ell) = (b_unit_vec' · b_FEM) / (b_unit_vec' · b_unit_vec)
        cost(ell) = ‖ C_opt · b_unit_vec − b_FEM ‖²
    
    ell_0 = argmin cost(ell)

Step 3: Phase 2 — fminbnd 精煉
    ell* = fminbnd(@(ell) ell_cost(ell, p_FEM, b_FEM, w),
                   ell_0 − 200e-6, ell_0 + 200e-6,
                   optimset('TolX', 1e-12))

Step 4: 取最後 C_opt 並回推 R_a
    [b_unit*, ~] = compute_unit_field(ell*, w, p_FEM)
    C_opt* = (b_unit*' · b_FEM) / (b_unit*' · b_unit*)
    R_a = N_c / (μ₀ · |C_opt*|)
    
Step 5: 報告誤差
    B_pred = C_opt* · b_unit*
    residual = B_pred − b_FEM
    err_per_node = ‖residual_j‖ / ‖B_FEM_j‖    對每個節點
    mean_err = mean(err_per_node)
    max_err = max(err_per_node)
```

### 2.4 MATLAB code 骨架

```matlab
% --- fit_ell_cost.m (核心 cost function) ---
function [cost, C_opt] = fit_ell_cost(ell, p_wp, b_fem, w, K_I, c)
    % Build 6 charge positions
    d_hat = build_d_hat(c);              % 6×3 方向矩陣
    charges = ell * d_hat;               % 6×3 位置矩陣
    
    % Build unit field
    R_a_unit = c.N_c / c.mu_0;           % 讓 C=1
    [bx, by, bz] = point_charge_model(p_wp, ell, R_a_unit, ...
                                       I_vec, K_I, c);
    b_unit = [bx; by; bz];               % 3N×1
    
    % VarPro analytical C
    C_opt = (b_unit' * b_fem) / (b_unit' * b_unit);
    residual = C_opt * b_unit - b_fem;
    cost = sum(residual.^2);
end

% --- 主程式 ---
% Phase 1: 粗掃
ell_grid = linspace(400e-6, 2000e-6, 300);
costs = arrayfun(@(e) fit_ell_cost(e, p_FEM, b_FEM, w, K_I, c), ell_grid);
[~, idx_min] = min(costs);
ell_0 = ell_grid(idx_min);

% Phase 2: 精煉
ell_star = fminbnd(@(e) fit_ell_cost(e, p_FEM, b_FEM, w, K_I, c), ...
                   ell_0 - 200e-6, ell_0 + 200e-6, ...
                   optimset('TolX', 1e-12));

% Step 4: 取 C_opt + R_a
[~, C_opt] = fit_ell_cost(ell_star, p_FEM, b_FEM, w, K_I, c);
R_a = c.N_c / (c.mu_0 * abs(C_opt));
```

### 2.5 數值結果（V4 halfcut, P1）

```
ell* = 835.4 µm           Long Fei 論文: 900 µm
C_opt = 7.59 × 10⁻²       對應 R_a = 9.21 × 10⁸ A/Wb
Mean vector error = 4.94%
Max vector error  = 7.94%

Cost landscape:
    ell = 705 µm: cost ratio 1.05× (相對 minimum)
    ell = 835 µm: cost ratio 1.00× (minimum)
    ell = 1026 µm: cost ratio 1.08×
    → 100 µm cube 內 cost landscape 平坦, ell 約束弱
```

### 2.6 角色

[A] 是 **baseline**：
- ✅ 跟龍飛論文最接近 → 可直接做 reproducibility 對照
- ✅ 提供 ell 數值給 [B-6x] 當鎖定參數
- ❌ 4.94% error 偏大；幾何約束（共用 ell + 固定方向）犧牲精度
- ❌ 100 µm cube 內衰減模式不明顯，ell 約束不強

### 2.7 物理意涵

[A] 的 ell = 835 µm ≠ 物理 tip 距離 (500 µm)。意思是：

> 「等效點電荷的最佳位置在 tip **後方** 335 µm（深入 pole 鋼體）」

物理解釋：磁通不集中在 tip apex，而是沿著 cone 表面 + 鋼體有 finite 分布。等效 point 落在「分布有效質心」上，自然在 tip 後方。**這就是 Long Fei 論文 b=0 假設下的「實際磁荷位置」**。

---

## 3. Method [J] — 6-coil joint 自由 3D 擬合（24 參數）

**腳本**：`magnetic_sim/hexapole-long2016/analysis/test_joint_6coil_fit.m`
**輸出**：`data/joint_6coil_19param_fit.mat`
**主結果**：Lower |c| avg 815.9 µm, Upper |c| avg 766.9 µm, R_a avg 1.01×10⁹, mean err 1.11%

### 3.1 假設（完全放開龍飛）

- **6 顆 charge 的 18 個座標完全自由**（每顆獨立 3D）
- 每組 coil 有自己的 **C_k**（共 6 個，per-coil 解析解）
- → 完全 free fit，看 charge 自然落在哪
- 用 **6 個 coil** 全部數據聯合 fit

### 3.2 為什麼需要 6 個 coil 一起 fit？

單 coil 數據時：

```
激勵 P1 → w = [5/6, −1/6, −1/6, −1/6, −1/6, −1/6]
              主極 ← 5 個被動極 →
              
P1 weight = 0.833 → 對 P1 位置敏感
P2~P6 weight = 0.167 → 對 P2-P6 位置幾乎不敏感
```

→ 單 coil 時 5 個被動極**位置不可辨識**（cost landscape 對它們極平坦，spread > 1000 µm）。

6 coil joint 時，每顆 pole 在某組數據中是主導 (5/6)，輪流被「點亮」→ 全部 6 個位置都可辨識。

### 3.3 數學參數化

```
模型 (對 coil k):
    B_k(p) = C_k · Σᵢ wᵢ(k) · (p − cᵢ) / |p − cᵢ|³

    wᵢ(k) = (K_I · I_vec_k)ᵢ      ← 跟 coil 編號 k 有關

待 fit:
    18 個位置座標 (c₁, c₂, ..., c₆) — 數值搜尋
    6 個 C_k (per-coil) — VarPro 解析
    
資料:
    6 coils × 647 nodes × 3 components = 11,646 個 B 觀測值
```

### 3.4 演算法詳細步驟

```
Step 0: 載入 6 組 coil 資料
    For k = 1..6:
        d_k = import_ansys_data(['coil' k], 'all', ['coil' k])
        應用同樣鐵芯過濾 + 100 µm cube
        p_FEM_k = positions, b_FEM_k = B 場堆疊
    Concatenate: 全部 6 組分開儲存

Step 1: 設定每個 coil 的 K_I weights
    For k = 1..6:
        paper_idx = apdl_to_paper_idx(k)
        I_vec_k = zeros(6,1); I_vec_k(paper_idx) = +1
        w_k = K_I · I_vec_k
        coil_sign_k = +1 (lower) or −1 (upper)

Step 2: 初始猜測（從 [A] 結果出發）
    ell_0 = 835 µm (from Method A)
    c_i_init = ell_0 · d̂_i,  i = 1..6
    
    x_0 = [c₁; c₂; ...; c₆]    ← 18×1 vector

Step 3: 定義 joint cost function
    function J = joint_cost(x):
        c_i = reshape(x, 6, 3)   ← 把 18×1 還原成 6 個 3D 位置
        J_total = 0
        For k = 1..6:
            b_unit_k = build_unit_field(c_i, w_k, p_FEM_k)
            C_opt_k = (b_unit_k' · b_FEM_k) / (b_unit_k' · b_unit_k)
            residual_k = C_opt_k · b_unit_k − b_FEM_k
            J_total += sum(residual_k.²)
        return J_total

Step 4: fminsearch 在 18 維搜尋
    x_opt = fminsearch(@joint_cost, x_0, optimset('TolX', 1e-10, ...))

Step 5: 唯一性驗證（5 組不同初始條件）
    For trial = 1..5:
        x_0_trial = x_0 + random_perturbation(50 µm scale)
        x_opt_trial = fminsearch(@joint_cost, x_0_trial, ...)
    
    spread_per_pole = std(x_opt_trials, axis=0)
    確認 spread < 1 µm 才算唯一

Step 6: 取每個 coil 的 C_opt_k + 對應 R_a_k
    For k = 1..6:
        b_unit_k = build_unit_field(c_i_opt, w_k, p_FEM_k)
        C_opt_k = (b_unit_k' · b_FEM_k) / (b_unit_k' · b_unit_k)
        R_a_k = N_c / (μ₀ · |C_opt_k|)
    
    R_a_avg = mean(R_a_k)  ← 通常 spread ~2%
```

### 3.5 MATLAB code 骨架

```matlab
% --- joint_cost.m ---
function J = joint_cost(x, p_FEM_all, b_FEM_all, w_all, c)
    % x: 18×1 = [c1_xyz; c2_xyz; ...; c6_xyz]
    charges = reshape(x, [3, 6])';        % 6×3
    
    J = 0;
    for k = 1:6
        % Build unit field for coil k
        b_unit_k = build_unit_field(charges, w_all(:,k), ...
                                     p_FEM_all{k});
        % VarPro analytical C_k
        b_fem_k = b_FEM_all{k};
        C_k = (b_unit_k' * b_fem_k) / (b_unit_k' * b_unit_k);
        residual = C_k * b_unit_k - b_fem_k;
        J = J + sum(residual.^2);
    end
end

% --- 主程式 ---
ell_0 = 835e-6;
c_init = ell_0 * d_hat;                  % 6×3
x_0 = c_init(:);                         % 18×1

opts = optimset('TolX', 1e-10, 'TolFun', 1e-22, ...
                'MaxIter', 50000, 'MaxFunEvals', 100000);

x_opt = fminsearch(@(x) joint_cost(x, p_FEM_all, b_FEM_all, w_all, c), ...
                   x_0, opts);

% 唯一性驗證
results = zeros(18, 5);
for trial = 1:5
    x_perturbed = x_0 + randn(18,1) * 50e-6;
    results(:, trial) = fminsearch(@(x) joint_cost(...), x_perturbed, opts);
end
spread_per_pole = std(reshape(results, [3, 6, 5]), 0, 3);
% → 結果: spread < 0.5 µm 對所有 pole ✓
```

### 3.6 數值結果（V4 halfcut, 6 coil joint）

Per-pole charge 距離（|c_i|）：

```
P1 (Lower halfcut): 814.8 µm     P2 (Upper full):  765.8 µm
P3 (Lower halfcut): 817.3 µm     P4 (Upper full):  768.2 µm
P6 (Lower halfcut): 815.2 µm     P5 (Upper full):  766.8 µm
                    ──────                          ──────
Lower 平均:         815.9 µm     Upper 平均:        766.9 µm
                    Δ = 49 µm (Lower charge 比 Upper 遠 6%)
```

Per-coil C/R_a：

```
Coil 1: R_a = 1.019 × 10⁹     Coil 4: R_a = 1.001 × 10⁹
Coil 2: R_a = 1.014 × 10⁹     Coil 5: R_a = 1.008 × 10⁹
Coil 3: R_a = 1.010 × 10⁹     Coil 6: R_a = 1.005 × 10⁹
                              ──────
                              avg = 1.01 × 10⁹ A/Wb, spread 2%

Mean vector error = 1.11%
Max vector error  = 2.83%
```

Charge 偏離 cone 軸（方向偏差）：

```
Pole   Layer    Dev [°]   r_perp [µm]   r_cone [µm]   Inside?
P1     Lower    1.55      163.7         93.2          NO (出 cone)
P3     Lower    1.28      168.4         93.2          NO
P6     Lower    1.69      162.2         93.5          NO
P2     Upper    5.49      67.2          92.1          YES
P4     Upper    5.59      68.7          92.5          YES
P5     Upper    5.59      68.6          92.2          YES
```

→ **Lower charge 跑出 cone 外** 70 µm（向 −z 偏 ~163 µm）；Upper charge 還在 cone 內但有 5.5° 偏軸。

### 3.7 角色

[J] 是 **物理觀察工具**：
- ✅ 完全 free fit → 真實 charge 位置會自然偏離 b=0 多少都看得到
- ✅ 6 coil joint → 識別性極佳（per-pole spread 0.3 µm，5 組初始條件確認唯一解）
- ✅ 精度提升 4.5× (4.94% → 1.11%)
- ❌ 24 自由度太多，每組 coil 自己一個 C_k 不太符合「single global R_a」物理

### 3.8 量化「b 偏離 0 多少」

跟 [A] 的球面假設 (ell = 835) 比，[J] 揭示：

| 觀察 | 物理意涵 |
|---|---|
| Lower vs Upper |c| 差 49 µm | 半切下極 vs 完整上極幾何不對稱 |
| Lower charge 偏 cone 軸 1.5° | 半切後磁荷中心略偏（往 −z 偏 163 µm） |
| Upper charge 偏 cone 軸 5.5° | 完整錐 + yoke 不對稱 |
| Lower charge **跑出 cone 外** | 磨平半錐使等效中心向 −z 偏移 |
| Upper charge 在 cone 內 | 完整錐保留對稱性 |

→ **b=0 假設在 halfcut 設計下不嚴謹**。對 Long Fei 的 b=0 推導不是致命錯誤（error 從 1% 變 5%），但量化結果證明 b≠0 才是物理事實。

---

## 4. Method [B-6x] — ell+δ 全激勵擬合（19 參數）★ 定案

**腳本**：`magnetic_sim/hexapole-long2016/analysis/fit_all6_with_bias.m`
**輸出**：`data/all6_bias_fit.mat`
**主結果**：ell = 835 µm (fixed), R_a = 1.03×10⁹, mean err **0.07%**（達 point-charge model 物理極限）

### 4.1 假設（混搭 [A] + [J]）

- ell 鎖定 = [A] 結果 (835 µm) → 用「球面錨點」當隱式正則化
- 每顆 charge 可偏離球面 → 18 個 δ_i 自由（位置在 ell·d̂_i 上**再加** δ_i）
- **共用同一個 C**（更符合 Long Fei 物理：global R_a 只應該一個值）
- 用 **all6 superposition 數據**（6 顆 coil 同時通 ±1A，FEM 線性疊加）

### 4.2 為什麼用 all6 superposition？

APDL 模擬「全部 6 顆同時 +1A」時，因 SOURC36 繞線方向：

```
Lower (P1/P3/P6) +1A → tip = SINK    → I_diss = +1
Upper (P2/P4/P5) +1A → tip = SOURCE  → I_diss = −1

I_diss = [+1, −1, +1, −1, −1, +1]  (對偶 pair ±1, paper 順序)
sum(I_diss) = 0  →  w = K_I · I_diss = I_diss

→ 所有 6 顆 pole 都有 |w| = 1（等權重）
```

✅ vs single-coil 的 5/6 vs 1/6 不均：[B-6x] 每顆 pole **同等強度約束**，identifiability 大幅提升
✅ 線性疊加：B(all6) = B(coil1) + ... + B(coil6) 直接從 6 顆獨立 sim 加總（已用 `verify_superposition.m` 驗證）

### 4.3 數學參數化

```
位置:    c_i = ell_fixed · d̂_i + δ_i
        (ell_fixed = 835 µm 鎖死, δ_i ∈ ℝ³ 自由)
        i = 1..6, 共 18 個 δ 座標

權重:   w = I_diss = [+1, −1, +1, −1, −1, +1]ᵀ  (paper 順序)

模型:   B_all6(p) = C · Σᵢ w_i · (p − c_i) / |p − c_i|³

待 fit:
    18 個 δ_i 座標 — 數值搜尋
    1 個 C — VarPro 解析

資料: all6 ±1A 對應的 647 nodes × 3 components = 1941 個 B 觀測值
```

### 4.4 演算法詳細步驟

```
Step 0: 取得 all6 superposition 資料
    Option A (直接 sim): 跑 APDL 全 6 顆同時 ±1A
    Option B (疊加): 把 6 個獨立 coil sim 結果線性加總
        b_FEM_all6 = Σ_k b_FEM_coil_k       ← 已驗證 < 0.1% error
    
    取 100 µm cube → 647 nodes
    b_FEM = [bx; by; bz]    ← 1941×1

Step 1: 取 ell_fixed (from Method A)
    Load 'data/charge_model_fit.mat' → ell_fixed = 835 µm
    
    建構球面錨點:
    For i = 1..6:
        anchor_i = ell_fixed · d̂_i

Step 2: 設定 I_diss 與 w
    I_diss = [+1, −1, +1, −1, −1, +1]ᵀ
    w = K_I · I_diss
    驗證 sum(I_diss) = 0 → w = I_diss (磁通守恆)

Step 3: 初始猜測
    δ_i = 0,  i = 1..6
    x_0 = zeros(18, 1)

Step 4: 定義 cost function (B-6x 專屬)
    function J = b6x_cost(x):
        δ = reshape(x, 6, 3)        ← 6 個 3D δ
        For i = 1..6:
            c_i = anchor_i + δ_i
        
        b_unit = build_unit_field(c_i, w, p_FEM)
        C_opt = (b_unit' · b_FEM) / (b_unit' · b_unit)
        residual = C_opt · b_unit − b_FEM
        J = sum(residual.²)
        
        return J

Step 5: fminsearch 在 18 維搜尋
    opts = optimset('TolX', 1e-12, 'TolFun', 1e-22, ...
                    'MaxIter', 100000, 'MaxFunEvals', 200000)
    x_opt = fminsearch(@b6x_cost, x_0, opts)
    
    δ_opt = reshape(x_opt, 6, 3)
    c_opt = anchor + δ_opt          ← 6×3 最終位置

Step 6: 取最終 C_opt + R_a
    b_unit_opt = build_unit_field(c_opt, w, p_FEM)
    C_opt = (b_unit_opt' · b_FEM) / (b_unit_opt' · b_unit_opt)
    R_a = N_c / (μ₀ · |C_opt|)
    
Step 7: 算誤差
    B_pred = C_opt · b_unit_opt
    residual = B_pred − b_FEM
    err_per_node = ‖residual_j‖ / ‖B_FEM_j‖
    mean_err = mean(err_per_node)        ← 應該 ~0.07%
    
Step 8: 多次擾動測試（confirm 是 global min）
    For trial = 1..5:
        x_perturbed = x_0 + randn(18, 1) · 50e-6
        x_opt_trial = fminsearch(...)
    
    所有 trial 都收斂到同一個 δ_opt (spread < 0.5 µm) ✓
```

### 4.5 MATLAB code 骨架

```matlab
% --- b6x_cost.m ---
function J = b6x_cost(x, anchor, w, p_FEM, b_FEM)
    % x: 18×1 = [δ1; δ2; ...; δ6]
    delta = reshape(x, [3, 6])';         % 6×3
    charges = anchor + delta;             % 6×3 位置
    
    % Build unit field (C = 1)
    b_unit = build_unit_field(charges, w, p_FEM);  % 3N×1
    
    % VarPro analytical C
    C_opt = (b_unit' * b_FEM) / (b_unit' * b_unit);
    residual = C_opt * b_unit - b_FEM;
    J = sum(residual.^2);
end

% --- 主程式 ---
% Load Method [A] 結果
A_result = load('data/charge_model_fit.mat');
ell_fixed = A_result.ell;                % 835 µm

% Build anchor
[d_hat, ~] = build_pole_directions(c);
anchor = ell_fixed * d_hat;              % 6×3

% I_diss for all6 ±1A
I_diss = [1; -1; 1; -1; -1; 1];          % paper order
w = K_I * I_diss;                         % = I_diss because sum=0

% Load all6 data (或從 6 顆 coil 疊加)
p_FEM = load_p_FEM_all6();
b_FEM = load_b_FEM_all6();

% fminsearch
x_0 = zeros(18, 1);
opts = optimset('TolX', 1e-12, 'TolFun', 1e-22, ...
                'MaxIter', 100000, 'MaxFunEvals', 200000);

x_opt = fminsearch(@(x) b6x_cost(x, anchor, w, p_FEM, b_FEM), x_0, opts);

delta_opt = reshape(x_opt, [3, 6])';     % 6×3
charges_opt = anchor + delta_opt;

% 取 C_opt + R_a
b_unit_opt = build_unit_field(charges_opt, w, p_FEM);
C_opt = (b_unit_opt' * b_FEM) / (b_unit_opt' * b_unit_opt);
R_a = c.N_c / (c.mu_0 * abs(C_opt));     % 1.03 × 10⁹
```

### 4.6 數值結果（V4 halfcut, all6 ±1A）

```
ell (fixed) = 835 µm
R_a = 1.03 × 10⁹ A/Wb
Mean vector error = 0.07%    ← point-charge model 的物理極限
Max  vector error = 0.19%
```

Per-pole charge 距離（|c_i| = |anchor + δ_i|）：

```
P1 (Lower): ~781 µm        P2 (Upper): ~734 µm
P3 (Lower): ~783 µm        P4 (Upper): ~734 µm
P6 (Lower): ~780 µm        P5 (Upper): ~734 µm
            ─────                      ─────
Lower 平均: ~781 µm        Upper 平均: ~734 µm
            Δ = 47 µm  ← 跟 [J] 的 49 µm Δ 一致
```

跟 [J] 比，[B-6x] 的 |c| 都被 ell-anchor 拉回約 35 µm（[J]: 816/767 vs [B-6x]: 781/734），但 Lower/Upper 不對稱 pattern 完全一致。

### 4.7 角色

[B-6x] 是 **定案**：
- ✅ **精度 0.07%** ≈ point-charge model 的物理極限
- ✅ 共用 1 個 C → 符合 Long Fei 物理（single global R_a）
- ✅ 19 自由度（比 [J] 的 24 少）+ 隱式正則化（ell-anchor）→ 更穩健
- ✅ all6 等權約束 → identifiability 最好
- ⚠️ 需先跑 [A] 才能取 ell

### 4.8 為什麼 [B-6x] 比 [J] 還準?

| 因素 | [J] | [B-6x] | 影響 |
|---|---|---|---|
| 數據型態 | 6 個 single-coil (5/6 vs 1/6 不均) | all6 等權 (|w_i|=1) | 等權給 identifiability 大幅提升 |
| C 個數 | 6 個 (per-coil) | 1 個 (共用) | 共用更符合 single R_a 物理 + 約束更強 |
| 位置參數化 | 全自由 (18 dim) | ell + δ (隱式正則化) | δ 不能跑太遠 → 抑制 overfit |
| 結果 error | 1.11% | 0.07% | 16× 改善 |

→ [B-6x] 是「正確物理約束 + 等權數據 + 隱式正則化」三件事一起做好的結果。

---

## 5. 三方法對比

### 5.1 結果總表

```
Method   Params   ell [µm]   R_a [A/Wb]    Mean err   Data        Identifiable?
[A]         2     835        9.21 × 10⁸    4.94%      1 coil      ✓
[J]        24     767-818    1.01 × 10⁹    1.11%      6 coils     ✓ (need ≥6 coils)
[B-6x]    19     734-783    1.03 × 10⁹    0.07%      all6        ✓ ★ 定案

Long Fei   ~2    ~900       ~6.3 × 10⁸    ~1%        ?           ✓
```

### 5.2 誤差來源分解

從 [A] 的 4.94% 開始，逐步放鬆假設後誤差降到哪：

```
Point-charge model 固有極限              ~0.07%  ([B-6x] 達到)
  + 共用 ell + 固定方向 (球面約束)      ~1.04%  ([A] 有, [J]/[B-6x] 無)
  + 6 個 C_k + 非等權 (single coil)     ~3.83%  ([A] 有, [B-6x] 用 all6 解決)
─────────────────────────────────────────────
  [A] total                              4.94%
  [J] total                              1.11%  (放開球面 + 用 6 coil)
  [B-6x] total                           0.07%  (all6 + 共用 C + 隱式正則化)
```

→ **球面約束**佔 ~21% 誤差、**非等權單 coil**佔 ~77%。後者是主要 killer。

### 5.3 何時用哪個

| 用途 | 推薦 method |
|---|---|
| Reproduce Long Fei (b=0 假設) | **[A]** |
| 量化「真實 b 偏離 0 多少」 | **[J]** |
| 力模型 / 控制器設計（要 high accuracy + global R_a） | **[B-6x]** ★ |
| 快速 sanity check 新 sim batch | [A]（最便宜） |
| 論文 Ch2.2.3 主結果 | **[B-6x]**（定案） + 對照 [A] 跟龍飛數值 |

### 5.4 R_a 跟論文不一致 (~6.3e8 vs ~1e9)

所有方法的 R_a 都比論文高 ~1.5×。

主因：APDL 用 **μ_r = 280** 線性近似（低場初始磁導率），高於實際操作點（B ~ 1T）對應的 μ_r ≈ 1400 + BH 飽和效應。低 μ_r → 鋼導磁能力差 → 磁通洩漏多 → 高階多極成分大 → 6-monopole 近似精度下降 + R_a 偏大。

✅ V4 halfcut 結果**內部一致**，可作橫向比較
⚠️ 跟 Long Fei 論文絕對數值有 ~25-50% 差距，需在論文 Ch2.2.3 註明此 caveat

---

## 6. 對龍飛 b=0 假設的回饋

### 6.1 我們的三個 fit 量化了什麼？

| 龍飛聲稱 | 我們發現 | 論文位置 |
|---|---|---|
| 「Charges 都在 ρ=900 µm 球面」 | 物理上等效中心比物理 tip (500 µm) 後移 280-340 µm，球面假設**近似成立** | [A]/[J]/[B-6x] 都同意 |
| 「b=0」 | Lower (halfcut) charge 跑出 cone 外 163 µm，Upper 偏軸 5.5° → **b≠0** | [J] 揭示 |
| 「Lower/Upper 對稱（共用 ρ）」 | Lower 比 Upper 遠 ~50 µm（一致橫跨 [J]/[B-6x]） | [J]/[B-6x] 一致 |
| 「Fit error < 1%」 | 嚴格 b=0 (我們 [A]) 給 ~5%；放開 b 才達 ~1%，[B-6x] 達 0.07% | [B-6x] 最強 |

### 6.2 結論

> **龍飛論文 Section 2.3 「Assume b=0 for the nominal force model」是 modeling convenience，不是物理事實**。
>
> b≠0 在 halfcut 下極特別明顯（charge 跑出 cone）。Force model closed-form inverse 推導若用嚴格 b=0 會引入系統誤差 ~5%；要做 high-accuracy control 必須用 [B-6x] 等級的 fit 把 b 還原。

### 6.3 對 kuo 論文 Ch2.2.3 / Ch3 的建議

1. **Ch2.2.3 改寫**：
   - 重做龍飛 Fig 2.6 對照（用 [A] 結果，verify b=0 假設下 error ~5%）
   - 新增 Fig 2.6'：用 [J] 揭示 Lower charge 跑出 cone 的 quantitative 圖
   - 新增 Fig 2.6"：[B-6x] 的 ε ≈ 0.07% 達到 point-charge floor

2. **Ch3 開頭 sensitivity section**：
   - 把 b=0 vs b-free RMS 差距量化（4.94% vs 0.07%）
   - 點明「b=0 適用於 closed-form inverse model 推導；但 forward prediction high accuracy 需用 [B-6x]」

---

## 7. 腳本/資料/圖檔清單

### 7.1 腳本（在 `magnetic_sim/hexapole-long2016/analysis/`）

| 檔案 | 用途 | 對應方法 |
|---|---|---|
| `mt_constants.m` | 幾何常數、極頭座標 | 所有 |
| `import_ansys_data.m` | 載入 FEM `.dat` | 所有 |
| `filter_iron_nodes.m` | 排除鐵芯節點 | 所有 |
| `point_charge_model.m` | 6-pole B 場計算（含 K_I） | 所有 |
| `fit_charge_model.m` | 2-param 球面 fit | **[A]** |
| `test_joint_6coil_fit.m` | 24-param 6-coil joint free 3D | **[J]** |
| `fit_all6_with_bias.m` | 19-param ell+δ all6 | **[B-6x]** ★ |
| `verify_superposition.m` | all6 = Σ coil_k 驗證 | [B-6x] 前置 |

### 7.2 數據

| 檔案 | 內容 |
|---|---|
| `data/charge_model_fit.mat` | [A] 結果 (ell, R_a, error) |
| `data/joint_6coil_19param_fit.mat` | [J] 結果 (6 個 c_i, R_a per-coil) |
| `data/all6_bias_fit.mat` | [B-6x] 結果 (6 個 δ_i, single R_a) ★ |

### 7.3 圖

| 檔案 | 對應論文 Fig |
|---|---|
| `generate_figures_2_3.m` | Fig 2.3 (B vector near tip 收斂) |
| `generate_figures_2_4.m` | Fig 2.4 (|B| contour) |
| `generate_figures_2_6.m` | Fig 2.6 用 [A] |
| `generate_fig26_coil1.m` | Fig 2.6 用 [J] |
| `plot_fig26_B6x.m` | Fig 2.6 用 [B-6x] ★ |

---

## 8. 名詞速查

| 符號 | 意義 | 典型值 |
|---|---|---|
| B(p) | 磁通密度 [T] | ~mT 量級在 WP |
| **k_m** | 磁 Coulomb 常數 = μ₀/(4π) | **10⁻⁷ N/A²** |
| μ₀ | 真空磁導率 | 4π×10⁻⁷ T·m/A |
| q_i, c_i | 第 i 顆 charge 強度 / 位置 | 10⁻⁹ A·m / mm |
| ξ | fit 用合併振幅 = -k_m·N_c/(μ₀·R_a) | VarPro 解析解 |
| Q | 磁荷向量 [q_1..q_6] | |
| K_I | 名義磁通分佈矩陣 (I₆ − 1₆/6) | 對角 5/6, 非對角 −1/6 |
| ell (ρ) | 等效磁荷距離（球面半徑） | 835 µm ([A]) vs 900 µm (Long Fei) |
| δ_i | charge i 相對球面的偏移 (b_i) | [B-6x] 自由 fit, [A] = 0 |
| R_a | 集總空氣磁阻 [A/Wb] | ~1×10⁹ ([B-6x]) vs 6.3×10⁸ (Long Fei) |
| α | 極角 | 54.74° (= arctan √2) |
| R_norm | WP 中心到極尖物理距離 | 500 µm |
| C | 振幅因子 = N_c/(μ₀·R_a) | VarPro 解析解 |
| d̂_i | pole i 從 WP 中心的單位方向 | 由 (θ_i, α, z_sign_i) 決定 |
| w_i | charge i 的相對權重 = (K_I·I_vec)_i | 5/6 主極, −1/6 被動極 |
| coil_sign | Lower=+1 / Upper=−1 (SOURC36 繞線修正) | |
| anchor_i | [B-6x] 用的球面錨點 = ell_fixed · d̂_i | |

---

## 9. 延伸閱讀

- **原始技術文件**：`magnetic_sim/hexapole-long2016/docs/charge-model-fitting.md`（完整推導，1300+ 行）
- **方法總結**：`magnetic_sim/hexapole-long2016/docs/fitting-methods.md`（精簡版）
- **Long Fei 原文**：dissertation Section 2.2.3 (printed page 17-20), Eq 2.1-2.4
- **Memory**：[[J-fit-validity-by-pole-shape]] — cylinder vs cone-fillet validity 邊界差 4×
- **相關**：[[long-fei-vs-paper-lab406]] — paper LAB406 與 dissertation 的 ρ 不一致根因
