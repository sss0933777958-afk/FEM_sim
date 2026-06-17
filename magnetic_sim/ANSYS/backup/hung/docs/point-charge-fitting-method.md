 # Point-Charge Model Fitting 方法說明

> 來源：Fei Long (龍飛) 2016 Dissertation, Section 2.2.3
> 整理目的：給人看的流程文件，方便套用到 Hung hexapole

---

## 1. 我們要做什麼

FEM 模擬會算出 workspace (WP) 附近數萬個節點的 B 場。但這些離散數據不能直接用於：
- 即時力計算（需要解析梯度）
- 控制器設計（需要電流 → 力的閉式關係）
- 快速掃描不同電流組合

**目標**：用一個只有 2 個參數的解析模型來近似 FEM 結果。

---

## 2. 物理直覺

### 為什麼可以用 point charge？

Pole tip 半徑 ≈ 40 um，到 WP 中心 ≈ 500 um。距離/尺寸比 = 12.5 >> 1。

從 WP 看過去，整根 pole tip「像一個點」——就像從遠處看一個帶電球體，場等效於所有電荷集中在球心。

### 模型概念圖

```
        WP 中心 (原點)
            ●
           ╱|╲
          ╱ | ╲
   q₁ ●  q₂●  ● q₃        ← 6 個磁單極 (point charges)
         ╲ | ╱              排列在半徑 ell 的球面上
          ╲|╱
   q₆ ●  q₅●  ● q₄
```

每個 charge 產生 Coulomb 型磁場（距離平方反比），6 個疊加 = 總場。

---

## 3. 數學模型

### 3.1 單個磁荷的場

磁荷 i 在位置 c_i，觀察點 p 的 B 場：

```
B_i(p) = k_m × q_i × (p − c_i) / |p − c_i|³
```

- k_m = mu_0 / (4 pi) = 10⁻⁷ [N/A²]（磁 Coulomb 常數）
- q_i = 磁荷大小 [A·m]
- c_i = 磁荷位置 [m]
- 分母是三次方（= 距離平方反比 × 方向向量）

### 3.2 6 極疊加

```
B(p) = Σᵢ k_m × q_i × (p − c_i) / |p − c_i|³
```

### 3.3 磁荷位置

6 個磁荷放在以 WP 中心為原點、半徑 ell 的球面上：

```
c_i = ell × d_hat_i
```

d_hat_i = 從 WP 指向第 i 根 pole 的單位方向向量：

```
d_hat_i = [ cos(theta_i) × sin(alpha),
            sin(theta_i) × sin(alpha),
            z_sign_i × cos(alpha) ]
```

- theta_i：方位角（P1=0°, P2=180°, P3=120°, P4=300°, P5=60°, P6=240°）
- alpha = 54.74°（magic angle）
- z_sign_i：Lower = −1, Upper = +1

**ell 不等於物理距離 (500 um)**。ell 是讓模型最佳逼近 FEM 的「等效距離」，擬合值 ~835-900 um。

### 3.4 電流 → 磁荷（Hopkinson 定律）

```
Q = (N_c / (mu_0 × R_a)) × K_I × I
```

- N_c = 70（線圈匝數）
- R_a = 空氣磁阻（待擬合）
- K_I = 磁通分佈矩陣（6×6，見下）
- I = 電流向量 [A]

### 3.5 完整模型

合併後，只有 **2 個待定參數**：

| 參數 | 意義 | 預期值 |
|------|------|--------|
| **ell** | 等效磁荷到 WP 的距離 | ~835-900 um |
| **R_a** | 空氣磁阻 | ~6.3×10⁸ ~ 1.0×10⁹ A/Wb |

其他都是已知量（N_c, k_m, alpha, theta_i, K_I, I）。

---

## 4. K_I 矩陣（磁通分佈矩陣）

### 怎麼來的

6 根 pole 共用一個 yoke，通電 1 根時磁通要從其他 5 根回來（磁通守恆 div B = 0）。

假設完美對稱：

```
通電 P1：
  P1 的磁通 = +5/6（主通路，flux sink）
  其他 5 根各 = −1/6（回流）
  合計 = 5/6 − 5×(1/6) = 0 ✓
```

### 矩陣形式

```
K_I = I₆ − (1/6) × 1₆

    = ┌  5/6  -1/6  -1/6  -1/6  -1/6  -1/6 ┐
      │ -1/6   5/6  -1/6  -1/6  -1/6  -1/6 │
      │ -1/6  -1/6   5/6  -1/6  -1/6  -1/6 │
      │ -1/6  -1/6  -1/6   5/6  -1/6  -1/6 │
      │ -1/6  -1/6  -1/6  -1/6   5/6  -1/6 │
      └ -1/6  -1/6  -1/6  -1/6  -1/6   5/6 ┘
```

### 物理意義

激勵 Coil k 時，6 個磁荷的權重 = K_I 的第 k 列。

例：激勵 P1 → w = [5/6, −1/6, −1/6, −1/6, −1/6, −1/6]

→ **P1 有大磁荷（主導），其他 5 根有小反向磁荷（回流）**。6 個 charges 同時存在，不是只有 1 個。

---

## 5. 擬合流程（Step by Step）

### Step 1：準備 FEM 數據

```
1. 跑 ANSYS 模擬（single coil, 1A, 70 turns）
2. 匯出 WP 附近的節點座標 (x, y, z) 和 B 場 (Bx, By, Bz)
3. 座標轉換到 WP 中心
4. 排除鐵心內的節點（不是空氣中的節點不能用）
5. 取小區域：|x|, |y|, |z| < 50 um → ~647 個節點
6. 堆疊成向量：b_fem = [Bx₁; By₁; Bz₁; Bx₂; By₂; Bz₂; ...] (1941 × 1)
```

### Step 2：建構 unit field

對給定的 ell，計算「C = 1 時的模型 B 場」：

```
b_unit(ell) = k_m × Σᵢ w_i × (p_j − ell × d_hat_i) / |p_j − ell × d_hat_i|³

堆疊成向量：b_unit = [Bx₁_model; By₁_model; Bz₁_model; ...] (1941 × 1)
```

### Step 3：Variable Projection（VarPro）解析求解 C

因為 B_model = C × b_unit，C 的最佳值有**解析解**：

```
              b_unit(ell)ᵀ × b_fem
C_opt(ell) = ─────────────────────────
              b_unit(ell)ᵀ × b_unit(ell)
```

→ **一個內積就算完，不用迭代**。

### Step 4：一維搜尋 ell

C 已解析搞定，只剩 ell 一個參數：

```
Step 4a: 粗掃
  ell = 400, 410, 420, ..., 2000 um
  對每個 ell → 算 C_opt → 算殘差 cost
  找到 cost 最低的區域

Step 4b: 精搜
  用 fminbnd 在粗掃最低附近精確搜尋
  → 得到最佳 ell*
```

### Step 5：回推 R_a

```
R_a = N_c / (mu_0 × |C_opt|)
```

### 流程圖

```
FEM 數據 (647 nodes × 3 components)
    │
    ▼
掃描 ell (400 ~ 2000 um)
    │
    ├── 每個 ell：
    │     ├── 算 6 個 c_i = ell × d_hat_i
    │     ├── 算 b_unit (unit field)
    │     ├── VarPro: C_opt = (b_unitᵀ b_fem)/(b_unitᵀ b_unit)  [解析]
    │     └── cost = ||C_opt × b_unit − b_fem||²
    │
    ▼
最佳 ell* (cost 最低)
    │
    ├── C_opt* → R_a = Nc/(mu0 × |C*|)
    └── 擬合完成：ell* 和 R_a* 就是你要的 2 個參數
```

---

## 6. 三種 fitting 方法

龍飛試了 3 種方法，逐步改進：

### [A] 基本版：共用 ell 球面（2 參數）

- 6 個 charge 都在半徑 ell 的球面上，方向固定
- 只用 1 組 coil 數據
- **自由度：2**（ell + R_a）
- 結果：ell = 835 um, 誤差 4.94%

### [J] 自由 3D：每極獨立位置（24 參數）

- 6 個 charge 的 (x,y,z) 完全自由
- 用全部 6 組 coil 數據
- 每組有自己的 C_k（6 個）
- **自由度：24**（18 座標 + 6 個 C）
- 結果：各極 ell 766~818 um, 誤差 1.11%

### [B-6x] 定案版：ell + delta（19 參數）★

- c_i = ell × d_hat_i + delta_i（球面基準 + 微調）
- ell 固定為 [A] 的 835 um，只擬合 18 個 delta
- 用 all6 疊加數據（所有 coil 同時通電）
- **自由度：19**（18 個 delta + 1 個 C）
- 結果：各極 |c| 734~783 um, 誤差 0.07%

### 對比表

```
方法      參數數   ell [um]        R_a          誤差    數據來源
───────────────────────────────────────────────────────────────
[A]          2    835 (共用)      9.2×10⁸      4.94%   1 coil
[J]         24    766~818 (各極)  1.0×10⁹      1.11%   6 coils
[B-6x] ★   19    734~783 (各極)  1.0×10⁹      0.07%   all6 疊加
───────────────────────────────────────────────────────────────
論文值      ~2    ~900            ~6.3×10⁸       ?       ?
```

---

## 7. 物理發現

### Lower vs Upper 不對稱

所有方法一致顯示：Lower pole 的等效磁荷距離 > Upper pole（差 ~50 um）。

原因：Lower pole 被銑平（支撐培養皿），移除部分材料 → 等效源中心外移。

### ell > R_norm

ell (~835 um) > 幾何距離 (500 um)。原因：真實 pole tip 不是點，磁通從整個尖端表面擴散。等效源在「更深處」。

### R_a 偏差

FEM 的 R_a (~10⁹) > 論文的 R_a (~6.3×10⁸)。主因：FEM 用線性材料 (mu_r=280)，低場近似使 WP 場強偏低 ~25%。

---

## 8. 套用到 Hung 模型

### 相同的部分
- 基本公式（Coulomb field, K_I, VarPro）
- alpha = 54.74°
- R_norm = 500 um
- 方位角 theta_i
- N_c = 70

### 需要改的部分

| 項目 | Long2016 | Hung | 怎麼改 |
|------|----------|------|--------|
| d_hat_i 方向 | 由 alpha 決定（對稱） | 由 TILT_UP/TILT_DN 決定（不對稱） | 用 Hung pole tip 座標算 |
| TILT_UP | N/A | 35° | 改 d_hat_i 公式 |
| TILT_DN | N/A | 5.71° | 改 d_hat_i 公式 |
| SPH_OFST | -12.711 mm | 需重算 | 從 APDL 提取 |
| FEM 數據路徑 | magnetic_sim/hexapole-long2016/results/ | magnetic_sim/hung/results/ | 改路徑 |
| 鐵心濾除 | cone filter | 需調整 | 依 Hung 幾何改 |

### 最小改動步驟
1. 複製 `magnetic_sim/hexapole-long2016/analysis/fit_charge_model.m` 到 `magnetic_sim/hung/analysis/`
2. 改 `results_dir` 路徑
3. 改 `d_hat_i`（用 Hung 的 tip 座標 / R_norm）
4. 改 `SPH_OFST`（Hung 的 WP 中心在 APDL 座標的 z 值）
5. 改鐵心濾除（cone filter 的幾何參數）
6. 跑 fit → 得到 Hung 的 ell 和 R_a

---

## 9. 關鍵公式速查

```
磁荷位置：      c_i = ell × d_hat_i
單荷 B 場：     B_i(p) = k_m × q_i × (p − c_i) / |p − c_i|³
6 極疊加：      B(p) = Σ B_i(p)
電流→磁荷：     Q = (N_c / (mu_0 × R_a)) × K_I × I
VarPro：        C_opt = (b_unitᵀ b_fem) / (b_unitᵀ b_unit)
回推 R_a：      R_a = N_c / (mu_0 × |C_opt|)
K_I（名義）：   K_I = I₆ − (1/6) × 1₆
```

---

## 10. 參考

- Long, F. (2016). Dissertation, Section 2.2.3, Eq. 2.1-2.4, 2.8, 2.19
- MATLAB 腳本：`magnetic_sim/hexapole-long2016/analysis/fit_charge_model.m`（[A] 方法）
- 完整技術文件：`magnetic_sim/hexapole-long2016/docs/charge-model-fitting.md`
- Fitting 方法對比：`magnetic_sim/hexapole-long2016/docs/fitting-methods.md`
