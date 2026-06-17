# Point-Charge Model Fitting：完整推導與記錄

> 來源：Long 2016 dissertation Section 2.2.3, Eq. 2.1-2.4
> 最後更新：2026-03-09

---

## 目錄

### Part I：論文完整 6-charge 模型（含 K_I 耦合）
1. [目標與動機](#1-目標與動機)
2. [物理背景：六極系統](#2-物理背景)
3. [數學模型：Eq 2.1-2.4](#3-數學模型)
4. [K_I 磁通分佈矩陣](#4-ki-磁通分佈矩陣)
5. [擬合策略：從 2D 降為 1D](#5-擬合策略)
6. [符號約定：Q 的負號](#6-符號約定)
7. [程式碼實現](#7-程式碼實現)
8. [6-charge 擬合結果](#8-6-charge-擬合結果)
9. [與論文的對比](#9-與論文的對比)
10. [已知限制](#10-已知限制)

### Part II：簡化單磁荷 3D 向量擬合（不含 K_I）
11. [動機與定位](#11-動機與定位)
12. [數學模型定義](#12-數學模型定義)
13. [3D 向量擬合完整推導](#13-3d-向量擬合完整推導)
14. [Variable Projection 等價性](#14-variable-projection-等價性)
15. [從 Q 反推 ell 和 R_a](#15-從-q-反推-ell-和-r_a)
16. [全 6 coil 結果](#16-全-6-coil-結果)
17. [Charge 幾何位置分析](#17-charge-幾何位置分析)

### Part III：兩種方法對比
18. [結構差異總表](#18-結構差異總表)
19. [為什麼 ell 和 R_a 不同](#19-為什麼結果不同)
20. [各自的意義與適用場景](#20-各自的意義)

### Part IV：參數空間探索與可辨識性實驗
21. [動機](#21-動機)
22. [測試的方法](#22-測試的方法)
23. [主結果表](#23-主結果表)
24. [ell 對比（每極）](#24-ell-對比每極)
25. [R_a 對比](#25-r_a-對比)
26. [可辨識性分析](#26-可辨識性分析)
27. [[J] 電荷位置分析](#27-j-電荷位置分析)
28. [誤差分解](#28-誤差分解)
29. [關鍵結論](#29-關鍵結論)

### 附錄
- [A. 程式碼與圖表索引](#a-程式碼與圖表索引)
- [B. 符號總表](#b-符號總表)

---

# Part I：論文完整 6-charge 模型（含 K_I 耦合）

## 1. 目標與動機

### 1.1 我們在做什麼

磁鑷（Magnetic Tweezers）有 6 根極頭（poles），每根纏繞 N_c = 70 匝線圈。通電後在 workspace（WP）中心產生磁場，用來操控磁珠。

我們已經用 ANSYS FEM 模擬了每根線圈單獨通 1A 電流時的完整 3D 磁場。現在要建立一個**解析模型**——用 6 個磁單極（point charges）來近似 FEM 的場。

### 1.2 為什麼需要解析模型

FEM 數據是離散的（~39 萬個節點），無法直接用於：
- 即時力計算（需要 B 場的解析梯度）
- 控制器設計（需要電流 → 力的閉式關係）
- 快速掃描不同電流組合

解析模型只有 **2 個待定參數**（`ell` 和 `R_a`），加上 K_I 矩陣所定義的 6 極耦合結構，擬合完成後可以在任意位置瞬間計算 B 場。

### 1.3 兩個待定參數

```
ell (等效磁荷距離)：6 個磁荷到 WP 中心的距離    預期 ~900 um（論文）
R_a (空氣磁阻)：    極頭 → WP 的等效空氣磁阻     預期 ~6.3x10^8 A/Wb（論文）
```

---

## 2. 物理背景

### 2.1 為什麼點電荷模型可行

極頭尖端半徑 = 40 um，WP 中心到極頭 = 500 um。
距離/尺寸比 = 500/40 = 12.5 >> 1。

在這個距離上，有限大小的極頭「看起來」就像一個數學上的點源——就像從遠處看一個帶電球體，場等效於所有電荷集中在球心。

### 2.2 六極系統的結構

6 根極頭排列在 WP 中心周圍，分上下兩層：

```
下層（Lower）：P1 (0°), P3 (120°), P6 (240°)    — z_wp < 0
上層（Upper）：P2 (180°), P4 (300°), P5 (60°)   — z_wp > 0
```

幾何參數：
```
R_norm    = 500 um    （WP 中心到極頭的幾何距離）
R_norm_xy = R_norm × sqrt(2/3) = 408 um    （xy 平面投影）
R_norm_z  = R_norm / sqrt(3) = 289 um       （z 方向投影）
alpha     = atan2(R_norm_xy, R_norm_z) = 54.74°  （極角）
```

所有極頭到 WP 中心的距離相等 = R_norm（因為 sqrt(R_norm_xy^2 + R_norm_z^2) = R_norm）。

### 2.3 座標系

- **APDL 座標**：原點在磁軛底面中心
- **WP 座標**：原點在 workspace 中心
- 轉換：`z_wp = z_apdl - SPH_OFST`，x 和 y 不變
- `SPH_OFST = -12.711 mm`
- B 場分量在兩個座標系中相同（純 z 平移）

---

## 3. 數學模型

### 3.1 Eq 2.1：磁通量 → 等效磁荷

```
q_i = Phi_i / mu_0
```

- q_i：第 i 根極頭的等效磁荷 [A*m]
- Phi_i：通過第 i 根極頭的磁通量 [Wb]
- mu_0 = 4*pi*10^-7 [T*m/A]

### 3.2 Eq 2.2：單磁荷的 Coulomb 場

```
B_i(p) = k_m × q_i × (p - c_i) / |p - c_i|^3
```

- p：場點位置 [m]
- c_i：第 i 個磁荷的位置 [m]
- k_m = mu_0 / (4*pi) = 10^-7 [N/A^2]

分母是 |p - c_i|^3（不是平方），因為：

```
(p - c_i) / |p - c_i|^3 = (1 / |p - c_i|^2) × (p - c_i) / |p - c_i|
                        = 距離平方反比衰減 × 單位方向向量
```

### 3.3 磁荷的位置

6 個磁荷放在以 WP 中心為原點、半徑為 ell 的球面上，方向與實際極頭相同：

```
c_i = ell × d_hat_i
```

其中 d_hat_i 是從 WP 中心指向第 i 根極頭的單位方向：

```
d_hat_i = [ cos(theta_i) × sin(alpha),
            sin(theta_i) × sin(alpha),
            z_sign_i × cos(alpha) ]
```

- theta_i：方位角（P1=0°, P2=180°, P3=120°, P4=300°, P5=60°, P6=240°）
- alpha = 54.74°：極角
- z_sign_i：下層=-1, 上層=+1

驗證：|c_i| = ell × sqrt(sin^2(alpha) + cos^2(alpha)) = ell。

**重要：ell ≠ 物理極頭距離 (500 um)。** ell 是讓點電荷模型最佳近似 FEM 的等效距離。擬合值 ~900 um > 500 um，因為真實極頭不是點，磁通從整個尖端表面擴散，等效源在更深處。

### 3.4 Eq 2.3：六極疊加

總場 = 6 個磁荷的場之和：

```
B(p) = sum_{i=1}^{6} k_m × q_i × (p - c_i) / |p - c_i|^3
```

### 3.5 Eq 2.4：電流 → 磁荷（Hopkinson 定律）

```
Q = (N_c / (mu_0 × R_a)) × K_I × I
```

推導：

**Step 1：磁動力（MMF）**
```
F_mmf_i = N_c × I_i
```
（N_c = 70，I_i = 第 i 個線圈的電流）

**Step 2：Hopkinson's Law（磁路 Ohm 定律）**
```
單一磁路：Phi = F_mmf / R_a
六極耦合：Phi = (N_c / R_a) × K_I × I
```

**Step 3：磁通 → 磁荷（Eq 2.1）**
```
Q = Phi / mu_0 = (N_c / (mu_0 × R_a)) × K_I × I
```

### 3.6 完整模型：電流 → B 場

合併 Eq 2.2-2.4：

```
B(p) = k_m × sum_{i=1}^{6} Q_i × (p - c_i) / |p - c_i|^3

其中 Q_i = (N_c / (mu_0 × R_a)) × (K_I × I)_i
     c_i = ell × d_hat_i
```

**已知量**：N_c, mu_0, k_m, theta_i, alpha, K_I, I
**待定參數**：ell（磁荷距離）, R_a（空氣磁阻）

---

## 4. K_I 磁通分佈矩陣

### 4.1 名義 K_I（Eq 2.8）

```
K_I = I_6x6 - (1/6) × 1_6x6

    = [ 5/6  -1/6  -1/6  -1/6  -1/6  -1/6 ]
      [-1/6   5/6  -1/6  -1/6  -1/6  -1/6 ]
      [-1/6  -1/6   5/6  -1/6  -1/6  -1/6 ]
      [-1/6  -1/6  -1/6   5/6  -1/6  -1/6 ]
      [-1/6  -1/6  -1/6  -1/6   5/6  -1/6 ]
      [-1/6  -1/6  -1/6  -1/6  -1/6   5/6 ]
```

對角項 = 5/6 ≈ 0.833，非對角項 = -1/6 ≈ -0.167。

### 4.2 物理意義

激勵 P1 時的磁通分佈：

```
K_I × [1,0,0,0,0,0]^T = [5/6, -1/6, -1/6, -1/6, -1/6, -1/6]^T
```

- P1 保留 5/6 的磁通（主通路，flux sink）
- 其餘 5 極各回流 -1/6（回流路，部分 flux source）

**磁通守恆**：每行加總 = 5/6 - 5×(1/6) = 0，即 ∑Phi_i = 0（div B = 0）。

### 4.3 K_I 的關鍵作用

K_I 是連結「哪個線圈通電」和「每個 pole 的磁荷強度」的矩陣：

```
激勵 Coil 1 (P1, I = [1,0,0,0,0,0]^T):
  Q_P1 = C × 5/6        → 主磁荷（dominant）
  Q_P2 = Q_P3 = ... = C × (-1/6)  → 回流磁荷（secondary）

  共 6 個 charges 同時存在，非只有 P1 一個
```

**這正是 Part II 單磁荷擬合缺失的部分**——它只擬合了 dominant charge，忽略了 5 個 secondary charges 的回流貢獻。

### 4.4 校準 K_I（Eq 2.19，論文實測值）

名義 K_I 假設完美對稱。論文通過 Hall sensor 校準得到的實測值顯示上下層不對稱：

```
Lower poles (P1, P3, P6): diagonal ≈ 0.60-0.63    (avg 0.61)
Upper poles (P2, P4, P5): diagonal ≈ 0.90-0.93    (avg 0.91)
```

原因：lower pole 被銑平（milled flat）以支撐培養皿，移除了部分材料，降低了磁通通過效率。

**我們的擬合使用名義 K_I**（假設對稱，diagonal = 5/6），尚未使用校準值。

---

## 5. 擬合策略

### 5.1 目標函數

給定 FEM 數據 {(p_j, B_j^FEM)}_{j=1}^{N}，最小化向量 B 的 SSE：

```
J(ell, R_a) = sum_{j=1}^{N} || B_model(p_j; ell, R_a) - B_FEM(p_j) ||^2
```

每個節點貢獻 3 個分量（Bx, By, Bz），殘差向量長度 = 3N。

### 5.2 關鍵觀察：B 對 R_a 是線性的

定義 C = N_c / (mu_0 × R_a)，則：

```
Q_i = -C × (K_I × I)_i          （含負號，見 §6）

B_model(p) = C × B_unit(p; ell)
```

其中 B_unit 是 C=1 時的模型場（只取決於 ell）。

**直觀理解**：改變 C（即改變 R_a）只改變所有磁荷的**振幅**，不改變**位置**。位置由 ell 決定。C 就像音量旋鈕——所有磁荷同時放大或縮小，場的空間形狀不變。

### 5.3 2D → 1D 降維

利用線性性，目標函數變成：

```
J(ell, C) = || C × b_unit(ell) - b_FEM ||^2
```

其中 b_unit 和 b_FEM 是 3N×1 向量（把所有節點的 Bx, By, Bz 串起來）。

對固定 ell，這是 scalar least-squares 問題。求導令其為零：

```
dJ/dC = 2 × (b_unit^T × b_unit) × C - 2 × (b_unit^T × b_FEM) = 0
```

解出：

```
                b_unit(ell)^T × b_FEM
C_opt(ell) = ─────────────────────────
                b_unit(ell)^T × b_unit(ell)
```

**這是解析解**——一個內積就算出，不需要迭代。

代入得到 ell 的一維 reduced cost：

```
                             [b_unit^T × b_FEM]^2
J*(ell) = ||b_FEM||^2  -  ──────────────────────────
                             b_unit^T × b_unit
```

### 5.4 一維搜尋 ell

**Phase 1：粗掃描**

在 ell ∈ [400, 2000] um 均勻取 300 點，對每個 ell：
1. 計算 6 個磁荷位置 c_i = ell × d_hat_i
2. 計算 B_unit（令 C=1，即 R_a = N_c/mu_0）
3. 用解析公式算 C_opt
4. 算 J* = ||C_opt × b_unit - b_FEM||^2
5. 找最小 → ell_0

**Phase 2：精煉**

用 MATLAB `fminbnd` 在 [ell_0 - 200um, ell_0 + 200um] 區間精煉。

**Phase 3：恢復 R_a**

```
R_a = N_c / (mu_0 × C_opt)
```

### 5.5 ell 和 R_a 的可分離性

直覺（燈泡類比）：

想像暗室裡一盞燈。你拿亮度計在不同位置量亮度，想反推：
- **瓦數** ↔ R_a（磁阻越小 → 磁荷越強 → 場越大）
- **距離** ↔ ell（磁荷到中心的距離）

只在一個位置量 → 無法區分「近處暗燈」和「遠處亮燈」。
在多個位置量 → 衰減模式不同（近燈衰減快，遠燈衰減慢）→ 可以分離。

數學上：場在兩點的比值消掉 R_a，只取決於 ell。所以 **ell 由空間衰減模式決定，R_a 由幅度決定**。

### 5.6 Joint 多 coil 擬合

將 6 個 coil 的數據堆疊起來，共享同一組 (ell, R_a)：

```
b_unit_all = [ b_unit_coil1; b_unit_coil2; ...; b_unit_coil6 ]
b_fem_all  = [ b_fem_coil1;  b_fem_coil2;  ...; b_fem_coil6  ]

C_opt = (b_unit_all^T × b_fem_all) / (b_unit_all^T × b_unit_all)
```

每個 coil 使用不同的 I_vec（只有該 coil 對應的 pole 通電），但 K_I、ell、R_a 共享。

---

## 6. 符號約定

### 6.1 我們的 Q vs 論文的 Q

**論文 Eq 2.4**：
```
Q = +(N_c / (mu_0 × R_a)) × K_I × I    （無負號）
```

**我們的程式碼**：
```
Q = -(N_c / (mu_0 × R_a)) × K_I × I    （有負號）
```

### 6.2 為什麼需要負號

在 APDL 模型中，激勵 Coil 1（P1）時，FEM 數據顯示 WP 中心的 B 場**指向 P1**——P1 是 **flux sink**。

正磁荷在 P1 位置 c_1 產生的場在原點：
```
B(0) = k_m × q × (0 - c_1) / |c_1|^3 = -k_m × q / ell^2 × d_hat_1
```
→ 指向遠離 P1 的方向。但 FEM 顯示 B 指向 P1，所以 Q_P1 必須是負的。

加上負號：Q_P1 = -(5/6)C < 0 → B 指向 P1 → 正確。

### 6.3 對擬合結果的影響

**數學上完全不影響**。不加負號時 C_opt 會變成負值來補償，最終 B_model 和 cost J* 完全相同。負號只確保 C_opt > 0 和 R_a > 0，使參數有物理可解釋性。

### 6.4 上下層的 coil_sign

上下層 pole 的 APDL 線圈繞向不同：
- Lower poles（P1, P3, P6）：正電流 → flux **INTO** pole → sink
- Upper poles（P2, P4, P5）：正電流 → flux **OUT OF** pole → source

為了統一模型（始終讓激勵 pole 為 sink），引入 coil_sign：

```
coil_sign = +1 (lower), -1 (upper)
I_eff = coil_sign × I_nominal
```

在 `fit_charge_model_6coil.m` 中：
```matlab
coil_sign = 1 - 2 * (1 - c.pole_is_lower(paper_idx));
I_vec(paper_idx) = coil_sign;
```

---

## 7. 程式碼實現

### 7.1 檔案結構

```
mt_constants.m          幾何常數、極頭位置、物理常數
filter_iron_nodes.m     幾何錐體模型排除鐵芯節點
point_charge_model.m    6-pole 模型 B 場計算（含 K_I）
fit_charge_model.m      單 coil (P1) 擬合
fit_charge_model_6coil.m  6 coils per-coil + joint 擬合
```

### 7.2 `point_charge_model.m` — 核心模型

```matlab
function [Bx, By, Bz] = point_charge_model(p_wp, ell, R_a, I_vec, K_I, c)
```

Step 1：計算 6 個磁荷位置
```
c_i = ell × [cos(theta_i)×sin(alpha), sin(theta_i)×sin(alpha), z_sign_i×cos(alpha)]
```

Step 2：計算磁荷向量
```
Q = -(N_c / (mu_0 × R_a)) × K_I × I_vec    （6×1 向量）
```

Step 3：6 個磁荷的 Coulomb 場疊加
```
B = k_m × sum_{i=1}^{6} Q_i × (p - c_i) / |p - c_i|^3
```

### 7.3 `fit_ell_cost` — 1D cost function

```matlab
function [cost, C_opt] = fit_ell_cost(ell, p_wp, b_fem, I_vec, K_I, c)
    R_a_unit = N_c / mu_0;               % 令 C = 1
    [bx, by, bz] = point_charge_model(p_wp, ell, R_a_unit, I_vec, K_I, c);
    b_unit = [bx; by; bz];               % 3N × 1

    C_opt = (b_unit' * b_fem) / (b_unit' * b_unit);   % 解析解
    residual = C_opt * b_unit - b_fem;
    cost = sum(residual.^2);
end
```

為什麼 R_a_unit = N_c/mu_0 使 C = 1？
```
C = N_c / (mu_0 × R_a) = N_c / (mu_0 × N_c/mu_0) = 1
```

### 7.4 鐵芯節點排除（`filter_iron_nodes.m`）

磁荷模型只適用於空氣節點。排除演算法：

1. 對每根極頭，計算從 tip 到每個節點的向量
2. 投影到極軸方向，算出沿軸距離 s 和垂直距離 r_perp
3. 錐體半徑 r_cone(s) = R_tip + s × (R_base - R_tip) / L_cone
4. 若 s > 0 且 r_perp < r_cone 且 s < L_cone → 鐵芯節點
5. 距極頭 < 100 um 的節點也排除

結果：排除 ~49,000 / 390,000 節點（~12.6%）。

---

## 8. 6-charge 擬合結果

### 8.1 單 Coil 擬合（fit_charge_model.m，僅 Coil 1 / P1）

```
擬合區域              ell [um]    R_a [A/Wb]    Mean err
Fit A (100 um cube):    835        9.21x10^8      4.87%
Fit B (R < 500 um):     799        1.06x10^9      7.16%
```

驗證（80 um cube）：Fit A mean error = 4.87%, Fit B = 7.16%。

### 8.2 6 Coil Per-coil 擬合（fit_charge_model_6coil.m）

每個 coil 獨立擬合 (ell, R_a)，使用完整 6-charge 模型 + K_I：

```
Pole   Layer   ell [um]   R_a [A/Wb]     Mean err
P1     Lower     836.1    9.273x10^8      4.88%
P3     Lower     836.2    9.286x10^8      4.86%
P6     Lower     835.5    9.234x10^8      4.91%
P5     Upper     748.8    1.186x10^9      4.57%
P2     Upper     747.5    1.188x10^9      4.61%
P4     Upper     746.5    1.184x10^9      4.62%
```

```
Lower avg (P1,P3,P6):  ell = 836 um,  R_a = 9.26x10^8
Upper avg (P2,P4,P5):  ell = 748 um,  R_a = 1.19x10^9
差異:                   delta_ell = -88 um (-10.5%)
```

### 8.3 Joint 擬合（共享 ell, R_a）

```
Joint (all 6):  ell = 835 um,  R_a = 9.22x10^8,  mean err = 4.87%
```

Joint 的 ell 偏向 lower poles（因為 lower 的 FEM 場強更大，在 least-squares 中權重更高）。

### 8.4 Cost Landscape 分析

Fit A（100 um cube）的 cost landscape 非常平坦：ell 從 705 到 1026 um，cost 變化不到 10%。這是因為 100 um cube 內的場近乎均勻，無法有效約束 ell（衰減模式不明顯）。

---

## 9. 與論文的對比

### 9.1 參數對比

```
                    ell [um]    R_a [A/Wb]    Center err
我們 (Joint):         835       9.22x10^8       ~5%
論文 (Long 2016):     900       6.30x10^8       <1%
```

### 9.2 差異的主要原因：material model

- **我們的 APDL**：murx = 280（線性，低場初始磁導率）
- **論文**：未明確指定，但操作點 (B~1T) 的 mu_r 可能 ~1400 或使用非線性 B-H

低 mu_r → 鐵芯導磁能力差 → 更多磁通側向洩漏 → 場的高階多極成分更大 → 6 磁單極近似精度下降。

### 9.3 擬合方法比較

```
                我們                          論文
策略:          1D scan + 解析 C             可能 2D lsqnonlin
K_I:           名義 eye(6) - ones(6)/6      名義 or 校準
Q 符號:        含負號（確保 R_a > 0）       不含負號（符號在 R_a 中吸收）
結果:          數學上等價
```

---

## 10. 已知限制

### 10.1 模型本身的限制

- 6 個磁單極是近似——真實場有高階多極成分
- 在 |p| → ell 時精度下降（場點接近磁荷，1/r^3 發散）
- 名義 K_I 假設完美對稱，但實際上下層不對稱

### 10.2 FEM 的限制

- murx = 280 可能不是最佳——導致場分佈偏離理想磁單極
- 線性材料模型（無飽和效應）

### 10.3 擬合區域 trade-off

```
小（100 um cube）：中心精度最高，但 ell 約束弱（cost 平坦）
大（R < 500 um）： ell 約束強，但遠處節點接近磁荷增加 error
```

---

# Part II：簡化單磁荷 3D 向量擬合（不含 K_I）

## 11. 動機與定位

### 11.1 這是什麼

一種**數據驅動的局部近似**：對每個 coil 獨立擬合一個 charge（自由 3D 位置 + 自由振幅），不使用 K_I，不考慮其他 5 個 pole 的耦合。

### 11.2 為什麼需要這個

6-charge 模型的 FEM 數據天然包含 6 極耦合效應。但在 WP 中心 100 um cube 內，dominant charge（5/6 的磁荷）壓倒性佔優，5 個 secondary charges（各 -1/6）的貢獻可忽略：

```
1-charge 和 6-charge 在 100 um cube 的 |B| error 幾乎相同（已用 test_simplified_model.m 驗證）
```

因此，1-charge 擬合可以作為：
1. 快速的局部近似（不需要知道 K_I）
2. 直接觀察每個 pole 的等效 charge 位置
3. 驗證 6-charge 模型的 dominant charge 假設

### 11.3 與 6-charge 模型的根本差異

```
6-charge（Part I）：                    1-charge（Part II）：
  ├── 物理模型（基於磁路理論）             ├── 數據擬合（純幾何 Coulomb fit）
  ├── K_I 定義 6 極耦合結構               ├── 不使用 K_I
  ├── 所有 charge 在半徑 ell 球面上        ├── Charge 在自由 3D 位置
  ├── 擬合 2 個參數 (ell, R_a)            ├── 擬合 4 個參數 (Q, cx, cy, cz)
  └── 結果是全局模型參數                   └── 結果只反映 dominant pole 特性
```

### 11.4 為什麼必須用向量 B 擬合

測試了三種方法：

```
方法                            向量誤差    問題
1D scalar (沿極軸, fit |B|)       4.95%    charge 方向被鎖死，限制精度
3D scalar (自由位置, fit |B|)     23.07%   方向完全崩壞（偏軸 10.6°）
3D vector (自由位置, fit B)        2.22%   大小和方向同時準確
```

scalar |B| 只攜帶「大小」資訊，無法約束 charge 在 3D 的位置。3D 自由 + scalar fit 會讓 charge 偏離極軸去 overfit magnitude。**向量 B 是唯一能正確約束 3D charge 位置的目標函數。**

---

## 12. 數學模型定義

### 12.1 模型

單一磁荷 Q 在位置 c = (cx, cy, cz) 產生的 Coulomb 向量場：

```
B(p) = Q × (p - c) / |p - c|^3
```

注意：這裡 Q 已經吸收了 k_m（不再另外乘），即 Q_1ch = k_m × Q_phys。

展開三個分量：
```
Bx_i = Q × (x_i - cx) / r_i^3
By_i = Q × (y_i - cy) / r_i^3
Bz_i = Q × (z_i - cz) / r_i^3

r_i = sqrt((x_i - cx)^2 + (y_i - cy)^2 + (z_i - cz)^2)
```

待擬合參數 = 4 個：Q（scalar）, c = (cx, cy, cz)。

---

## 13. 3D 向量擬合完整推導

### 13.1 參數分離：B 對 Q 線性，對 c 非線性

固定 c 時，Bx_i, By_i, Bz_i 都正比於 Q。定義 **basis vector** m(c)：

```
            ⎡ (x_1-cx)/r_1^3 ⎤
            ⎢      ...       ⎥
            ⎢ (x_N-cx)/r_N^3 ⎥   ← Bx 部分，N 個
            ⎢ (y_1-cy)/r_1^3 ⎥
m(c) =     ⎢      ...       ⎥   ← By 部分，N 個       (3N × 1)
            ⎢ (y_N-cy)/r_N^3 ⎥
            ⎢ (z_1-cz)/r_1^3 ⎥
            ⎢      ...       ⎥   ← Bz 部分，N 個
            ⎣ (z_N-cz)/r_N^3 ⎦
```

FEM 數據堆疊：
```
b = [ Bx_1, ..., Bx_N, By_1, ..., By_N, Bz_1, ..., Bz_N ]^T    (3N × 1)
```

模型：
```
B_model = Q × m(c)
```

### 13.2 Q 的解析解（Normal Equation）

目標函數：
```
J(Q, c) = || Q × m(c) - b ||^2
        = Q^2 × (m^T m) - 2Q × (m^T b) + (b^T b)
```

J 對 Q 是二次（m^T m > 0）。求導：
```
dJ/dQ = 2Q × (m^T m) - 2 × (m^T b) = 0
```

解出：
```
              m(c)^T × b
Q_opt(c) = ──────────────
              m(c)^T × m(c)
```

### 13.3 Reduced Cost Function

把 Q_opt 代回 J：
```
                             [m^T × b]^2
J*(c) = ||b||^2  -  ────────────────────
                         m^T × m
```

推導：
```
J* = (B_val^2/A^2) × A - 2 × (B_val/A) × B_val + ||b||^2
   = B_val^2/A - 2 B_val^2/A + ||b||^2
   = ||b||^2 - B_val^2/A

其中 A = m^T m,  B_val = m^T b
```

### 13.4 fminsearch 優化 c

用 MATLAB fminsearch 優化 c = (cx, cy, cz)，最小化 J*(c)：

```matlab
cost_fn = @(pos) vector_charge_cost(pos, px, py, pz, bvec_fem);
opts = optimset('TolX', 1e-10, 'TolFun', 1e-22, 'MaxIter', 10000);

% 兩組初始猜測（charge sign ambiguity）
[c1, f1] = fminsearch(cost_fn, +c_1d_initial, opts);
[c2, f2] = fminsearch(cost_fn, -c_1d_initial, opts);
charge_opt = argmin(f1, f2);
```

初始猜測來自 1D scalar 掃描（charge 沿極軸方向）。

### 13.5 完整演算法流程

```
輸入：N 個場點 {(x_i, y_i, z_i), (Bx_i, By_i, Bz_i)}^FEM

Step 1：初始猜測
  ├── 1D 掃描 b ∈ [0, 800] um，charge 沿極軸
  ├── 對每個 b，c = (R_norm + b) × d_pole
  ├── 算 scalar a^2 和 cost
  └── 取最小 → c_0 = 初始猜測

Step 2：3D 向量擬合
  ├── cost(c) = vector_charge_cost(c)
  │     ├── 算 m(c) ← 3N×1 basis vector
  │     ├── Q = (m^T × b) / (m^T × m)  ← 解析解
  │     └── return || Q×m - b ||^2
  ├── fminsearch(cost, +c_0) → (c_1, f_1)
  ├── fminsearch(cost, -c_0) → (c_2, f_2)
  └── charge_opt = argmin(f_1, f_2)

Step 3：恢復參數
  ├── Q_opt = (m^T × b) / (m^T × m)
  ├── ell = |c*|
  └── R_a = -(k_m × N_c × 5 × coil_sign) / (6 × mu_0 × Q_opt)

輸出：c*, Q_opt, ell, R_a
```

### 13.6 cost function 實現

```matlab
function cost = vector_charge_cost(pos, px, py, pz, bvec_fem)
    dx = px - pos(1);  dy = py - pos(2);  dz = pz - pos(3);
    r3 = (dx.^2 + dy.^2 + dz.^2).^(3/2);
    mvec = [dx ./ r3; dy ./ r3; dz ./ r3];     % 3N x 1 basis
    Q = (mvec' * bvec_fem) / (mvec' * mvec);    % analytical Q
    residual = Q * mvec - bvec_fem;
    cost = sum(residual.^2);
end
```

---

## 14. Variable Projection 等價性

### 14.1 問題

先解析求 Q 再優化 c（3 參數），和同時優化 (Q, cx, cy, cz)（4 參數），結果是否相同？

### 14.2 回答：數學上嚴格等價

這是 Variable Projection 方法（Golub & Pereyra, 1973）。

聯合問題：min_{Q,c} J(Q,c) = ||Q × m(c) - b||^2

最優解 (Q*, c*) 必須滿足：
```
dJ/dQ = 0  →  Q* = (m(c*)^T × b) / (m(c*)^T × m(c*))    ... (i)
dJ/dc = 0                                                    ... (ii)
```

條件 (i) 恰好就是解析解。因為 J 對 Q 正定（二次），Q 的最優解唯一：
```
min_c J*(c) = min_c [min_Q J(Q,c)] = min_{Q,c} J(Q,c)
```

### 14.3 數值驗證（Coil 1 / P1）

```
Separated (3 param):  c = [673.29, -1.56, -521.46] um,  Q = -6.290855e-09
Joint (4 param):      c = [673.29, -1.56, -521.46] um,  Q = -6.290855e-09

Cost 差異: 6.44 x 10^-20  →  完全等價
```

### 14.4 實務優勢

```
                 分離 (VarPro)              聯合 4 參數
搜索維度:         3D (cx, cy, cz)           4D (Q, cx, cy, cz)
Q 精度:           解析（機器精度）           數值解（受收斂影響）
尺度問題:         無（c ~ 10^-4）           Q ~ 10^-9 vs c ~ 10^-4
結果:             完全等價                   完全等價
```

---

## 15. 從 Q 反推 ell 和 R_a

### 15.1 ell（等效磁荷距離）

```
ell = |c*| = sqrt(cx^2 + cy^2 + cz^2)
```

### 15.2 R_a（空氣磁阻）

1-charge 的 Q 已吸收 k_m：
```
1-charge:  B = Q_1ch × (p - c) / |p - c|^3
6-charge:  B = k_m × Q_phys × (p - c) / |p - c|^3
→ Q_1ch = k_m × Q_phys
```

從 6-charge 模型（§3.5），激勵 pole j 時 dominant charge 的 Q_phys：
```
Q_phys = -(N_c / (mu_0 × R_a)) × (K_I × I)_j × coil_sign
       = -(N_c / (mu_0 × R_a)) × (5/6) × coil_sign
```

注意 coil_sign：
- Lower (P1,P3,P6): coil_sign = +1，Q_phys < 0（sink）
- Upper (P2,P4,P5): coil_sign = -1，Q_phys > 0（source）

代入：
```
Q_1ch = k_m × Q_phys = -(k_m × N_c × 5 × coil_sign) / (6 × mu_0 × R_a)
```

解出 R_a：
```
              k_m × N_c × 5 × coil_sign
R_a = - ─────────────────────────────────
                6 × mu_0 × Q_1ch
```

### 15.3 數值驗算（Coil 1 / P1, lower, coil_sign = +1）

```
Q_1ch = -6.290855 x 10^-9
k_m   = 1 x 10^-7
N_c   = 70
mu_0  = 4*pi*10^-7 = 1.2566 x 10^-6

分子 = -(10^-7 × 70 × 5 × 1) = -3.5 x 10^-5
分母 = 6 × 1.2566x10^-6 × (-6.2909x10^-9) = -4.744 x 10^-14

R_a = -3.5x10^-5 / (-4.744x10^-14) = 7.38 x 10^8 A/Wb
```

---

## 16. 全 6 coil 結果

### 16.1 Per-coil 結果（fit_3d_vector_all_coils.m）

```
Coil  Pole  Layer   ell [um]   R_a [A/Wb]   Vec err   Dev [°]
 1    P1    Lower    851.6     7.38x10^8     2.22%     2.50
 2    P3    Lower    854.2     7.34x10^8     2.16%     2.18
 3    P6    Lower    853.2     7.33x10^8     2.25%     2.59
 4    P5    Upper    746.3     9.31x10^8     0.67%     1.26
 5    P2    Upper    745.8     9.34x10^8     0.69%     1.40
 6    P4    Upper    746.5     9.28x10^8     0.68%     1.20
```

### 16.2 層平均

```
Lower (P1,P3,P6):  ell = 853.0 um,  R_a = 7.35x10^8,  vec_err = 2.21%
Upper (P2,P4,P5):  ell = 746.2 um,  R_a = 9.31x10^8,  vec_err = 0.68%
All 6:             ell = 799.6 um,  R_a = 8.33x10^8,  vec_err = 1.44%
```

### 16.3 Q 符號驗證

```
P1  Lower(sink)  Q = -6.2909e-09  (negative) ✓
P3  Lower(sink)  Q = -6.3196e-09  (negative) ✓
P6  Lower(sink)  Q = -6.3260e-09  (negative) ✓
P5  Upper(src)   Q = +4.9849e-09  (positive) ✓
P2  Upper(src)   Q = +4.9664e-09  (positive) ✓
P4  Upper(src)   Q = +4.9976e-09  (positive) ✓
```

Lower Q_avg = -6.31x10^-9, Upper Q_avg = +4.98x10^-9。
|Q_lower| > |Q_upper| 因為 lower 的磁通更集中（半錐體 tip area 更小）。

---

## 17. Charge 幾何位置分析

### 17.1 Charge 與 pole 錐體的關係

問題：擬合出來的 charge 位置在物理 pole 錐體內嗎？

**Upper poles (P2, P4, P5)：YES — charge 在錐體內部**

```
P5: 沿軸 s=247 um,  r_perp=19 um,  r_cone=89 um  →  inside (r_perp << r_cone)
P2: 沿軸 s=246 um,  r_perp=18 um,  r_cone=89 um  →  inside
P4: 沿軸 s=245 um,  r_perp=20 um,  r_cone=88 um  →  inside
```

Upper poles 的 tip-axis 方向（水平，指向 WP 中心的反方向）與 d_pole 方向（WP→tip）的夾角只有 1.3°，幾乎一致。charge 沿 d_pole 延伸自然落在錐體內。

**Lower poles (P1, P3, P6)：NO — charge 在錐體外部**

```
P1: 沿軸 s=265 um,  r_perp=233 um,  r_cone=92 um  →  outside (偏出 141 um)
```

Lower poles 的 pole 物理軸方向（水平）與 d_pole 方向（從 WP 中心向右下方）夾角 = 35.3°。charge 沿 d_pole 延伸偏離錐體。

### 17.2 為什麼 upper 在內、lower 在外

```
Upper pole:
  d_pole 方向 ≈ pole 物理軸方向    (夾角 1.3°)
  → charge 沿極軸延伸 → 在錐體內

Lower pole:
  d_pole 方向 ≠ pole 物理軸方向    (夾角 35.3°)
  → charge 向斜下方延伸 → 偏出錐體
```

Lower poles 的 tip 在 WP 中心下方（z_wp = -289 um），但極軸是水平的（沿 xy 平面向外），所以 d_pole（從中心到 tip）有很大的向下分量，但極軸沒有。

### 17.3 物理解釋

charge 位置是**數學抽象**，不是物理實體：

```
WP center ── 500 um ──→ Tip ── ~350 um ──→ Charge (ell ≈ 850 um)
(origin)                 (pole surface)       (equivalent point source)
```

等效位置在 tip 後方（更遠離 WP 中心），因為真實的磁通源不是 tip 上的一個點，而是整個尖端表面甚至更深的鐵芯。ell > R_norm（850 > 500）是合理的。

---

# Part III：兩種方法對比

## 18. 結構差異總表

```
                        6-charge (Part I)           1-charge (Part II)
─────────────────────────────────────────────────────────────────────
基礎              物理模型（磁路理論）         數據擬合（Coulomb fit）
K_I 耦合          有（6 極互相耦合）           無（忽略其他 5 極）
Charge 數量        6 個                         1 個（只取 dominant）
位置約束           全在半徑 ell 球面，           完全自由 3D
                  方向 = 極頭方向
振幅機制           K_I 固定 6 個 charge          Q 自由擬合
                  的相對比例
自由參數           2 (ell, R_a)                 4 (Q, cx, cy, cz)
                                               等效 3（Q 解析）
擬合函數           向量 B（3N 分量）            向量 B（3N 分量）
```

## 19. 為什麼結果不同

### 19.1 ell 的差異

```
                    6-charge      1-charge
Lower (P1,P3,P6):    836           853     (+17 um)
Upper (P2,P4,P5):    748           746     (-2 um)
Joint:                835            —
```

6-charge 的 per-coil ell 是在 **K_I 結構約束下**的最佳值：charge 必須在球面上，方向鎖定。
1-charge 的 ell 是 **無約束 3D 位置**的結果（ell = |c*|），可以偏軸。

Lower poles 差異更大（17 um）因為 1-charge 的 charge 有 2.5° 偏軸，|c*| 被偏軸分量推大。

### 19.2 R_a 的差異

```
                    6-charge      1-charge
Lower (P1,P3,P6):  9.26x10^8     7.35x10^8   (-21%)
Upper (P2,P4,P5):  1.19x10^9     9.31x10^8   (-22%)
```

1-charge 的 R_a 系統性更小（場更強），因為：

1. **缺少 secondary charges 的回流**：6-charge 模型中 5 個 secondary charges 部分抵消 dominant charge 的場（同向回流減弱淨場）。1-charge 沒有這些回流，所以需要更小的 R_a（更強的 Q）來產生同樣的 |B|。

2. **ell 的差異**：1-charge 的 ell 稍大 → 離 WP 更遠 → 需要更強的 Q 來補償距離。但 R_a ∝ 1/Q，所以 R_a 更小。

定量估算：
```
B ∝ Q / ell^2 ∝ 1 / (R_a × ell^2)

6-charge 的等效單 charge 貢獻 ≈ (5/6) × total → 等效 R_a_eff = R_a × 6/5
但 secondary charges 的抵消效應使 R_a 需要更大以匹配相同 |B|
```

### 19.3 向量誤差的差異

```
                    6-charge      1-charge
Lower (P1,P3,P6):   ~4.9%         2.21%
Upper (P2,P4,P5):   ~4.6%         0.68%
```

1-charge 的向量誤差**更低**，因為：
- 自由 3D 位置有 3 個 DOF vs 6-charge 的 1 個（ell 在球面上，方向鎖死）
- 1-charge 只需要擬合 dominant pattern，不需要兼顧 secondary charges 的回流方向

但這不代表 1-charge 是更好的模型——它只是在 100 um cube 的局部擬合更靈活。在遠場或多 coil 疊加時，1-charge 模型就不適用了。

---

## 20. 各自的意義

### 20.1 6-charge 模型（Part I）的用途

- **控制器設計**：K_I 結構允許任意電流組合 → 預測任意激勵下的 B 場
- **力模型**：F = quadratic form in I（Eq 2.6-2.7），需要 K_I
- **全局模型**：在整個 WP 球（R < 500 um）內有效
- **物理可解釋**：ell 和 R_a 有明確的磁路意義

### 20.2 1-charge 模型（Part II）的用途

- **局部診斷**：觀察每個 pole 的等效 charge 位置和強度
- **上下層差異分析**：直接看到 lower vs upper 的 ell 和 R_a 差異
- **方法驗證**：確認 dominant charge 假設在 WP 中心成立
- **快速近似**：不需要 K_I 就能估計單 coil 的中心場

### 20.3 互補關係

```
1-charge 揭示：
  → 每個 pole 的 dominant charge 方向偏軸 1-3°（6-charge 假設零偏軸）
  → Lower 和 Upper 的 ell 差異 107 um（6-charge joint 用單一 ell 折衷）
  → 中心場被 dominant charge 壓倒性主導（secondary 可忽略）

6-charge 提供：
  → K_I 結構使模型可預測任意電流組合
  → 全局適用性（不限於中心小區域）
  → 物理參數（R_a）可連結到磁路分析
```

---

# Part IV：參數空間探索與可辨識性實驗

> 2026-03-09 新增。使用 Coil 1 (P1) 數據和全 6 coil 數據，系統性測試不同約束層級對擬合精度和解唯一性的影響。

## 21. 動機

Model A（Part I）的 4.94% 向量誤差有多少來自球面約束（6 極共享 ell）？放鬆約束能否在不失唯一性的前提下提升精度？

## 22. 測試的方法

共 7 種，按自由度從少到多排列：

```
ID    方法                                    自由度    資料
─────────────────────────────────────────────────────────────
[A]   共用 ell, 球面, K_I                     2         1 coil
[A']  6-coil joint, 共用 ell+R_a              2         6 coils
[a]   per-pole ell_i, 固定方向, K_I           7         1 coil
[V]   3D vector 單電荷, 無 K_I                4         1 coil
[c]   共用 L + 角度偏移, K_I                  14        1 coil
[B]   自由 3D 位置, K_I                       19        1 coil
[J]   6-coil joint, 自由 3D, per-C_k          24        6 coils
```

### 22.1 Method [a]: Per-pole ell_i

每個電荷沿極頭方向自由滑動：`c_i = ell_i × d_hat_i`。方向固定，距離自由。

### 22.2 Method [c]: Shared L + angular offsets

6 個電荷固定在半徑 L 的球面上，但各自可偏離原方向：
`c_i = L × direction(theta_i + dtheta_i, alpha + dalpha_i)`

### 22.3 Method [J]: 6-coil joint free 3D

共用 18 個位置座標，6 組 coil 數據聯合擬合。每個 coil 有自己的 C_k（解析求解），KI_w_k 通過 apdl_to_paper_idx 映射正確建構。

## 23. 主結果表

```
ID   Method                             Npar  Mean err  Max err  Unique?
────────────────────────────────────────────────────────────────────────
[A]  Shared ell, sphere (1 coil)           2    4.94%    7.94%    YES
[A'] 6-coil joint, shared ell+R_a          2    5.23%    8.26%    YES
[a]  Per-pole ell_i, fixed dir (1 coil)    7    1.62%    3.29%    YES
[V]  3D vector, 1 charge, no K_I           4    2.22%    4.35%    YES
[c]  Shared L + angular offsets (1 coil)  14    0.36%    1.05%    NO
[B]  Free 3D, K_I (1 coil)               19    0.30%    1.00%    NO
[J]  6-coil joint, free 3D, per-C_k      24    1.11%    2.83%    YES
```

## 24. ell 對比（每極）

```
ID     P1(L)   P3(L)   P6(L)   P2(U)   P4(U)   P5(U)   L avg   U avg
─────────────────────────────────────────────────────────────────────
[A]    835.4   835.4   835.4   835.4   835.4   835.4   835.4   835.4
[A']   788.3   788.3   788.3   788.3   788.3   788.3   788.3   788.3
[a]    999.1   759.9   761.5   785.5   717.6   715.7   840.2   739.6
[V]    851.6    —       —       —       —       —       851.6    —
[J]    814.8   817.5   815.2   765.8   768.2   766.8   815.9   766.9
Diss.   —       —       —       —       —       —      ~900    ~900
                                                        [um]    [um]
```

## 25. R_a 對比

```
ID     Method                                |R_a| [A/Wb]
──────────────────────────────────────────────────────────
[A]    Shared ell, sphere (1 coil)           9.21 x 10^8
[A']   6-coil joint, shared ell+R_a          1.02 x 10^9
[a]    Per-pole ell_i (1 coil)               7.10 x 10^8
[V]    3D vector, no K_I (1 coil)            7.38 x 10^8
[J]    6-coil joint, free 3D (avg)           1.01 x 10^9
Diss.  Long 2016                            ~6.3  x 10^8
```

[J] per-coil |R_a| 極為一致（spread ~2%）：
- Lower (P1, P3, P6): 1.010~1.019 x 10^9
- Upper (P2, P4, P5): 1.001~1.008 x 10^9

## 26. 可辨識性分析

以 5 組不同初始條件測試，位置 spread (std) [um]：

```
ID   Method                      P1     P2     P3     P4     P5     P6
──────────────────────────────────────────────────────────────────────
[A]  Shared ell (1D scan)         0      0      0      0      0      0
[a]  Per-pole ell_i               0      0      0      0      0      0
[B]  Free 3D, 1 coil            522   1180   1257    732    586   1425
[J]  6-coil joint, free 3D      0.4    0.4    0.3    0.3    0.3    0.3
```

[B] 單 coil 不可辨識：P1 權重 5/6 = 0.833，其餘各 1/6 = 0.167。被動極位置不敏感。
[J] 6-coil joint 恢復唯一性：每個極在某組數據中是主導的（5/6 權重），spread 改善 1000~4000 倍。

### 26.1 為什麼徑向可辨識、角度不可辨識

- 改距離（方案 [a]）：直接改變 1/r^3 衰減率，一階效應，可辨識
- 改角度（方案 [c]）：遠處偏轉幾度對近場幾乎無影響，二階效應，不可辨識
- 加上 L 和 C 的 trade-off（推遠 → C 變大補償），角度模型徹底崩潰

## 27. [J] 電荷位置分析

### 27.1 方向偏差

```
Pole  Layer   Dev [deg]
P1    Lower     1.55
P3    Lower     1.28
P6    Lower     1.69
P2    Upper     5.49
P4    Upper     5.59
P5    Upper     5.59
```

Lower 幾乎不偏（~1.5°），Upper 偏 ~5.5°。

### 27.2 Charge 與 pole cone 的關係

```
Pole  Layer  ell[um]  s[um]  r_perp[um]  r_cone[um]  Inside cone?
P1    Lower   814.8   269.5    163.7       93.2          NO
P3    Lower   817.5   269.6    168.4       93.2          NO
P6    Lower   815.2   271.0    162.2       93.5          NO
P2    Upper   765.8   263.9     67.2       92.1         YES
P4    Upper   768.2   266.2     68.7       92.5         YES
P5    Upper   766.8   264.8     68.6       92.2         YES
```

- s = charge 在 tip 後方多遠（沿 pole axis）
- r_perp = 離 pole axis 的垂直距離
- r_cone = 該處 cone 半徑（r_perp < r_cone 才在 cone 內）

Upper poles 的 charge 在 cone 內部。Lower poles 的 charge 跑出 cone 外，垂直偏移 ~163 um，方向是 **-z（向下）**。

原因：Lower pole 被磨平（VSBV 切掉上半部），等效磁荷中心從原本的 cone 軸線向下偏移，遠離被移除的部分。Upper pole 是完整錐體，charge 留在 cone 內。

## 28. 誤差分解

```
Model A 的 4.94% 誤差來源：

  Point-charge model floor       ~0.30%  (模型本身極限)
  + Fixed direction cost         ~0.81%  (方向鎖定在極頭軸)
  + Shared ell cost              ~3.33%  (強迫 6 極等距)
  + Single-coil bias             ~0.50%  (只有 1 coil 約束 6 poles)
  ───────────────────────────────────────
  = Model A total                ~4.94%
```

各方法移除的誤差源：
- [a] per-pole ell_i → 1.62%：移除 shared ell cost (-3.33 pp, 67%)
- [V] free 3D + no K_I → 2.22%：移除 shared ell + direction cost
- [J] 6-coil joint 3D → 1.11%：移除 shared ell + single-coil bias
- [B] free 3D, 1 coil → 0.30%：移除全部（但解不唯一）

## 29. 關鍵結論

1. **球面約束是最大的誤差源**（3.33 pp / 67%），不是 K_I 耦合或模型本身
2. **per-pole ell_i [a] 是最佳甜蜜點**：+5 params, -67% error, 解唯一
3. **6-coil joint [J] 使 19-param 完全可辨識**：spread 從 ~1000 um → 0.3 um
4. **Lower/Upper 不對稱在所有方法中一致**：Lower ell ~816 um, Upper ~767 um
5. **Lower charge 偏離 cone（向下 ~163 um）**：磨平半錐體使等效中心下移
6. **R_a 全部高於論文值 (~6.3e8)**：與 murx=280 線性材料模型有關

---

# 附錄

## A. 程式碼與圖表索引

### A.1 程式碼

```
檔案                               功能                              相關章節
─────────────────────────────────────────────────────────────────────────
mt_constants.m                    幾何常數、極頭位置、物理常數       §2-3
filter_iron_nodes.m               鐵芯節點排除                       §7.4
point_charge_model.m              6-pole 模型 B 場計算（含 K_I）     §7.2
import_ansys_data.m               載入 ANSYS FEM 數據                §7
fit_charge_model.m                單 coil 6-charge 擬合              §5, §8.1
fit_charge_model_6coil.m          6 coils per-coil + joint 擬合      §5.6, §8.2-8.3
test_simplified_model.m           1-charge scalar vs 6-charge 比較   §11.2
test_3d_charge_fit.m              3D 向量擬合 Coil 1 + 驗證          §12-15, §17
fit_3d_vector_all_coils.m         3D 向量擬合 all 6 coils            §16
test_19param_fit.m                19-param 自由 3D vs 2-param 比較   §22, §26
test_joint_6coil_fit.m            6-coil joint 19-param 擬合 [J]     §22.3, §23-27
```

### A.2 數據

```
data/charge_model_fit.mat          Coil 1 的 Fit A / Fit B 結果
data/charge_model_6coil_fit.mat    6-coil per-coil + joint 結果
data/vecfit_3d_all_coils.mat       3D 向量擬合 all 6 coil 結果
data/joint_6coil_19param_fit.mat   6-coil joint 19-param [J] 結果
```

### A.3 圖表

```
fig2_6a.png                  FEM vs 6-charge quiver overlay
fig2_6b.png                  6-charge error scatter (R<400um)
fig2_6a_vecfit.png           FEM vs 1-charge vector quiver overlay
fig2_6b_vecfit.png           1-charge vector error scatter (100um cube)
fig_xz_geometry_coil1.png    xz-plane: pole shape, tip, charge positions
```

---

## B. 符號總表

```
符號          名稱                            值 / 定義                 來源
────────────────────────────────────────────────────────────────────────────
B             磁通密度                         [T]                      Eq 2.2
q_i           第 i 極等效磁荷                  = Phi_i / mu_0 [A*m]    Eq 2.1
Q             磁荷向量 [q_1..q_6]              [A*m]                    Eq 2.3
I             電流向量 [I_1..I_6]              [A]                      Eq 2.4
K_I           磁通分佈矩陣                     = I_6 - (1/6)×1_6       Eq 2.8
ell           等效磁荷距離                     ~900 um (論文)           fitted
R_a           集總空氣磁阻                     ~6.3x10^8 (論文)        fitted
N_c           每極線圈匝數                     70                       p.14
k_m           磁 Coulomb 常數                  = mu_0/(4*pi) = 10^-7   p.18
mu_0          真空磁導率                       = 4*pi*10^-7            —
R_norm        WP 到極頭的物理距離              500 um                   p.14
R_norm_xy     xy 平面投影                      408 um                   computed
R_norm_z      z 投影                           289 um                   computed
alpha         極角                             54.74°                   computed
d_hat_i       極頭方向單位向量                 [cos(th)sin(a), ...]    §3.3
coil_sign     上下層符號修正                   +1 lower, -1 upper      §6.4
C             等效振幅因子                     = N_c/(mu_0×R_a)        §5.2
Q_1ch         1-charge Q (含 k_m)              = k_m × Q_phys          §12.1
```
