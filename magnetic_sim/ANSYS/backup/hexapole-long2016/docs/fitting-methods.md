# Point-Charge Fitting Methods: [A] → [J] → [B-6x]

> 最終擬合方法文檔。定案方法為 [B-6x]。
> 最後更新：2026-03-21

---

## 0. 共同基礎

### 0.1 物理模型

6 根 pole 各有一個等效磁荷（point charge）。
通電時，每個磁荷在 workspace (WP) 中心產生 Coulomb 型磁場，疊加得到總場。

### 0.2 基本公式

單個磁荷 i 在位置 c_i，觀察點 p 的磁場：

    B_i(p) = k_m × w_i × (p - c_i) / |p - c_i|^3

其中：
- k_m = mu_0 / (4*pi) = 10^-7  (磁 Coulomb 常數)
- w_i = 該磁荷的權重（由 K_I 矩陣和電流決定）
- c_i = 磁荷位置 [m]
- p = 場點位置 [m]

6 極疊加：

    B(p) = C × k_m × SUM_i [ w_i × (p - c_i) / |p - c_i|^3 ]

C = N_c / (mu_0 × R_a) 是整體振幅因子，包含：
- N_c = 70（線圈匝數）
- R_a = 空氣磁阻 [A/Wb]（待擬合）

### 0.3 K_I 矩陣與權重

K_I = I_6 - ones(6)/6（名義磁通分佈矩陣）

當 Coil k 通 1A 電流時，6 個磁荷的權重向量：

    w = K_I × I_vec

例如激勵 P1：I_vec = [1,0,0,0,0,0]

    w = [5/6, -1/6, -1/6, -1/6, -1/6, -1/6]

→ P1 權重 5/6（主導），其餘各 -1/6（被動回流）

注意：程式中用 -w_i（負號），因為 APDL 繞線使受激勵極為 flux sink。

### 0.4 共同數據

來源：ANSYS FEM 模擬（magnetostatic，murx=280）

每個 coil 的數據處理流程：
1. 載入 ANSYS 節點座標和 B 場 → 390,579 nodes
2. 轉換到 WP 中心座標：z_wp = z_apdl - SPH_OFST
3. 排除鐵心節點（geometric cone filter）→ 341,428 nodes
4. 取 100 um 立方體（|x|,|y|,|z_wp| < 50 um）→ 647 nodes

擬合數據：
- 位置：px, py, pz（647 × 3）
- B 場：bx, by, bz（647 × 3）
- 堆疊成向量：b_fem = [bx; by; bz]（1941 × 1）

### 0.5 幾何參數

- Pole tip 到 WP 中心：500 um（全部一樣）
- 擬合區域半寬：50 um
- alpha = 54.74°（pole tip 在球面上的極角）
- 6 極的 azimuthal angle：P1=0°, P2=180°, P3=120°, P4=300°, P5=60°, P6=240°
- Lower (P1,P3,P6)：磨平半錐體；Upper (P2,P4,P5)：完整錐體

### 0.6 APDL 到 Paper 的索引映射

APDL coil 1,2,3,4,5,6 = Paper P1,P3,P6,P5,P2,P4

K_I 和磁荷位置都用 paper 索引 (P1-P6)。載入 APDL coil k 的數據時，
I_vec 必須在 paper 空間中正確設定：

    paper_idx = apdl_to_paper_idx(k)
    I_vec(paper_idx) = 1

### 0.7 Variable Projection (VarPro)

所有方法共用的核心技巧：模型 B = C × b_unit，C 的最佳解是線性最小二乘：

    C_opt = (b_unit' × b_fem) / (b_unit' × b_unit)

這將 optimizer 的維度減少 1（或 6，if per-coil），只需搜尋位置參數。

---

## 1. Method [A]：共用 ell 球面擬合（2 參數）

**腳本：** `fit_charge_model.m` → `data/charge_model_fit.mat`

### 1.1 場景

論文的標準方法。6 個磁荷全部放在半徑 ell 的球面上，方向由 pole tip 幾何決定。
只用 Coil 1 (P1) 的 FEM 數據。

### 1.2 設定

磁荷位置：

    c_i = ell × d_hat_i

其中 d_hat_i 是第 i 極的單位方向向量（固定不變）：

    d_hat_i = [cos(theta_i)*sin(alpha), sin(theta_i)*sin(alpha), z_sign_i*cos(alpha)]

約束：6 個磁荷共用同一個 ell，方向完全由幾何決定。

### 1.3 推導

1. 建構 unit field（令 C=1）：b_unit(ell) = f(c_i = ell × d_hat_i)
2. VarPro 解析求解 C_opt(ell)
3. 一維搜尋 ell：粗掃 400~2000 um + fminbnd 精細搜尋
4. 回推 R_a = N_c / (mu_0 × C_opt)

### 1.4 自由度

- 優化：1（ell）+ 解析：1（C → R_a）= 總計 2

### 1.5 結果

    ell = 835.4 um
    R_a = 9.21 × 10^8 A/Wb
    向量誤差：mean = 4.94%, max = 7.94%

### 1.6 角色

[A] 提供 ell baseline，被 [B-6x] 用作固定的球面基線距離。

---

## 2. Method [J]：6-coil joint 自由 3D（24 參數）

**腳本：** `test_joint_6coil_fit.m` → `data/joint_6coil_19param_fit.mat`

### 2.1 場景

6 個磁荷的位置完全自由（18 個座標），同時使用全部 6 組 coil 的 FEM 數據聯合擬合。
每組 coil 有自己的 C_k（解析求解）。

### 2.2 為什麼需要 6 coil

單 coil 時，只有受激勵極有大權重（5/6），其餘 5 極各只有 1/6。
6 coil 時，每個 coil 使不同的極成為主導 → 所有位置可辨識。

### 2.3 設定

磁荷位置：c_i = (x_i, y_i, z_i)，完全自由，共 18 個座標。
資料：6 × 647 nodes × 3 components = 11,646 個數據點。
Per-coil C_k 讓每組數據各自調整強度。

### 2.4 推導

1. 對每個 coil k，建構 w_k = K_I × I_vec_k
2. 建構 unit field：b_unit_k = f(c_i, w_k)
3. VarPro 解析求解 C_k
4. 總 cost = SUM_k cost_k
5. fminsearch 在 18 維空間搜尋，5 組初始條件驗證唯一性

### 2.5 自由度

- 優化：18 + 解析：6（C_1..C_6）= 總計 24
- 數據/優化參數比：11646 / 18 = 647:1

### 2.6 可辨識性

5 組初始條件 → spread 0.3~0.4 um（極佳）。

### 2.7 結果

Per-pole |c| [um]:

    P1 (Lower): 814.8     P2 (Upper): 765.8
    P3 (Lower): 817.3     P4 (Upper): 768.2
    P6 (Lower): 815.2     P5 (Upper): 766.8

    Lower avg: 815.9 um | Upper avg: 766.9 um

    R_a avg: ~1.01 × 10^9 A/Wb（per-coil spread ~2%）
    Mean error: 1.11%

方向偏差：Lower ~1.5°, Upper ~5.5°
Lower charges 在 cone 外部（向 -z 偏移 ~163 um，磨平效應）。

---

## 3. Method [B-6x]：全激勵 b!=0 擬合（19 參數）★ 定案

**腳本：** `fit_all6_with_bias.m` → `data/all6_bias_fit.mat`

### 3.1 場景

用 all6 數據（所有 6 coil 同時 +1A），ell+delta 參數化，1 個共用 C。

### 3.2 為什麼 all6 = 每對 ±1

APDL 全 +1A 時，由於 SOURC36 winding 方向：
- Lower pole (P1,P3,P6) +1A → tip = SINK → model I_diss = +1
- Upper pole (P2,P4,P5) +1A → tip = SOURCE → model I_diss = -1

    I_diss = [+1, -1, +1, -1, -1, +1]
    sum(I_diss) = 0 → w = K_I × I_diss = I_diss
    所有 6 個 pole 都有 |w| = 1

### 3.3 優勢

- vs [J]：共用 C（更符合論文模型），自由度少（19 vs 24）
- vs single-coil：所有極等權約束（|w|=1 vs 5/6 和 1/6）→ 更好的 identifiability
- ell+delta 參數化提供隱含正則化（球面錨點），比直接 3D 位置更穩定

### 3.4 設定

模型：c_i = ell × d_hat_i + delta_i（18 delta params）
ell 固定為 [A] 的結果（835 um）
C 共用 1 個（VarPro 解析求解）

資料：all6 疊加場（線性材料 → B(all6) = B(coil1) + ... + B(coil6)，已驗證）

### 3.5 推導

1. 固定 ell = 835 um，建構 d_hat_i
2. 對每組 delta_i，計算 c_i = ell × d_hat_i + delta_i
3. 建構 unit field，VarPro 解析求解 C_opt
4. fminsearch 在 18 維空間搜尋 delta
5. 回推 R_a = N_c / (mu_0 × |C_opt|)

### 3.6 自由度

- 優化：18（delta）+ 解析：1（C → R_a）= 總計 19

### 3.7 結果

    ell (fixed): 835 um
    R_a = 1.03 × 10^9 A/Wb
    Mean fitting error = 0.07%

Per-pole |c| [um]:

    P1 (Lower): ~781     P2 (Upper): ~734
    P3 (Lower): ~783     P4 (Upper): ~734
    P6 (Lower): ~780     P5 (Upper): ~734

    Lower avg: ~781 um | Upper avg: ~734 um

位置和 [J] 一致（差異 ~45 um 以內）。

---

## 4. 方法對比

### 4.1 結果總表

    Method   Params  ell [um]     |R_a|       Mean err  Data
    ──────────────────────────────────────────────────────────
    [A]         2    835 (all)    9.21e8       4.94%     1 coil
    [J]        24    766~818      1.01e9       1.11%     6 coils
    [B-6x]     19    734~783      1.03e9       0.07%     all6   ★
    ──────────────────────────────────────────────────────────
    Diss.      ~2    ~900         ~6.3e8         ?         ?

### 4.2 誤差來源分解

    Point-charge 模型 floor     ~0.07%   (固有極限，[B-6x] 接近)
    + 方向固定 + 共用 ell       ~1.04%   ([A] 有、[J] 無)
    + 6 C_k vs 1 C + 非等權     ~3.83%   ([A] 有、[B-6x] 無)
    ────────────────────────────────────
    [A] total                   ~4.94%
    [J] total                   ~1.11%
    [B-6x] total                ~0.07%

### 4.3 Lower vs Upper 不對稱

所有方法一致顯示：Lower avg |c| > Upper avg |c|（差異 ~45-50 um）。
原因：Lower poles 的半錐體（磨平）使等效磁荷中心外移。

### 4.4 R_a 偏差

所有方法的 R_a ~ 0.8-1.0 × 10^9，高於論文的 ~6.3 × 10^8。
主因：murx=280（低場線性近似），使 WP 中心場比論文低 ~25%。

---

## 5. 腳本索引

    Method   Script                    Purpose                          Output
    [A]      fit_charge_model.m        b=0 sphere (2 params)            charge_model_fit.mat
    [J]      test_joint_6coil_fit.m    free 3D, 6 coils (24 params)     joint_6coil_19param_fit.mat
    [B-6x]   fit_all6_with_bias.m      ell+delta, all6 (19 params) ★    all6_bias_fit.mat

共同依賴：
- mt_constants.m — 幾何參數
- import_ansys_data.m — 載入 FEM 數據
- filter_iron_nodes.m — 排除鐵心節點
- point_charge_model.m — [A] 的 B 場計算

執行順序：
1. fit_charge_model.m → charge_model_fit.mat（提供 ell baseline）
2. test_joint_6coil_fit.m → joint_6coil_19param_fit.mat（中間參考）
3. fit_all6_with_bias.m → all6_bias_fit.mat（定案）

圖表腳本：
- generate_figures_2_3.m — 論文 Fig 2.3（純 FEM）
- generate_figures_2_4.m — 論文 Fig 2.4（純 FEM）
- generate_figures_2_6.m — 論文 Fig 2.6（用 [A]）
- generate_fig26_coil1.m — Fig 2.6 variant（用 [J]）
- plot_fig26_B6x.m — Fig 2.6 variant（用 [B-6x]）
- verify_superposition.m — 線性疊加驗證
