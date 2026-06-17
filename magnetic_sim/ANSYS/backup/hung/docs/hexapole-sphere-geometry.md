# Hung Hexapole — Sphere Geometry Framework

## Pole 編號與配對

**上層 poles:** P2, P4, P5
**下層 poles:** P1, P3, P6

| Pair | Lower Pole | Upper Pole |
|------|-----------|------------|
| I    | P1        | P2         |
| II   | P3        | P4         |
| III  | P6        | P5         |

## 座標系

| 座標系 | 符號 | 顏色 | 說明 |
|--------|------|------|------|
| Assembly | (x_a, y_a, z_a) | 藍 | CAD/ANSYS 座標系 |
| Magnetic | (x_m, y_m, z_m) | 紅 | 軸沿 3 對 pair 方向，正交 |

## 球面幾何條件

球心 = working point (WP)，6 個 tip 在球面上。

| # | 條件 | 效果 |
|---|------|------|
| 1 | 每組 pair 連線過球心 | L_i = -U_i（對徑點） |
| 2 | 上層 3 tip 共平面，下層 3 tip 共平面 | 上下兩平行平面，對稱於球心 |
| 3 | 球半徑 R = 0.5 mm | 每個 tip 到 WP 距離 = 0.5 mm |
| 4 | 方位角等距 120° | 3 對 pair 在水平投影上均勻分布 |

## 數學框架

### 座標系設定

z 軸垂直於上下平面（從下層指向上層），球心在原點。

### 參數化

上層 tip 位於球面與上層平面的交圓上：

```
圓半徑：r = R sin(θ)
平面高度：h = R cos(θ)
```

θ = 球面極角（z 軸到 tip 方向的夾角），α = 90° - θ = 仰角（水平面到 tip 的仰角）。

### 6 個 tip 座標

方位角 φ_1 = 0°, φ_2 = 120°, φ_3 = 240°（可整體旋轉）。

**上層：**

```
U_i = ( r cos(φ_i),  r sin(φ_i),  h )    i = 1, 2, 3
```

**下層（對徑）：**

```
L_i = ( -r cos(φ_i), -r sin(φ_i), -h )   i = 1, 2, 3
```

其中 r² + h² = R²。

### 自由度

條件 1-4 後剩餘 **1 個自由度**：θ（或等價地 h）。

| θ | α (仰角) | h (mm) | r (mm) | 特性 |
|---|----------|--------|--------|------|
| 小 | 大 | 大 | 小 | 上下平面遠離赤道，水平圓小 |
| 54.74° | **35.26°** | 0.2887 | 0.4082 | **Magic angle** |
| 大 | 小 | 小 | 大 | 上下平面靠近赤道，水平圓大 |

### Pair 方向向量

第 i 對 pair 的方向向量（U_i → L_i 或反向）：

```
d_i = U_i - L_i = ( 2r cos(φ_i), 2r sin(φ_i), 2h )
```

歸一化後：

```
d_hat_i = ( sin(θ) cos(φ_i), sin(θ) sin(φ_i), cos(θ) )
```

## 正交條件 → 固定 θ

d_hat_i · d_hat_j = 0 要求：

```
sin^2(θ) cos(φ_i - φ_j) + cos^2(θ) = 0
```

φ 間距 120° → cos(120°) = -1/2，代入得：

```
1 - 3cos^2(θ)/2 = 0  →  cos^2(θ) = 2/3  →  θ = 54.74°
```

α = 90° - 54.74° = **35.26° (magic angle)**

### MATLAB 驗證

```
Pair 1 vs 2: 90.00° ✓
Pair 1 vs 3: 90.00° ✓
Pair 2 vs 3: 90.00° ✓
Tip-to-tip: 1.000 mm ✓
```

---

## Pole 擺放

### 傾斜角

| 層 | Poles | 傾斜角 | 方向 |
|----|-------|--------|------|
| 上層 | P2, P4, P5 | 35° | 從 tip 沿軸往上往外延伸至 block |
| 下層 | P1, P3, P6 | 5.71° | 從 tip 沿軸往下往外延伸至 block |

### Pole 尺寸 (from Mag_Pole_Bottom.STEP)

| 參數 | 值 | 說明 |
|------|-----|------|
| 圓柱半徑 | 3.175 mm | Diameter 6.35mm |
| 錐段長度 | 15.875 mm | semi-angle 11.31° |
| 圓柱段長度 | 27.125 mm | 全圓 15mm + D型 12.125mm |
| 總長 | 43 mm | |
| 削平 | tip 起 28mm | 通過軸心的平面，去下半，flat face 朝上 |
| Tip | 尖錐頂點 | STEP 無 fillet |

### Tip 座標 (mm, magic angle)

| Pole | Type | Azimuth | Tip (x, y, z) mm |
|------|------|---------|-------------------|
| P1 | lower | 0° | (+0.4082, 0, -0.2887) |
| P2 | upper | 180° | (-0.4082, 0, +0.2887) |
| P3 | lower | 120° | (-0.2041, +0.3536, -0.2887) |
| P4 | upper | 300° | (+0.2041, -0.3536, +0.2887) |
| P5 | upper | 60° | (+0.2041, +0.3536, +0.2887) |
| P6 | lower | 240° | (-0.2041, -0.3536, -0.2887) |

> 座標與 APDL 程式碼一致：P1 在 +X 方向（azim=0°），upper 在 +Z，lower 在 -Z。

### APDL 模型

腳本：`magnetic_sim/ANSYS/backup/hung/apdl/variants/MT_Hung_SphereModel.txt`
圖檔：`magnetic_sim/ANSYS/backup/hung/figures/hexapole_sphere_6poles.png`

---

## 座標系轉換

**Assembly 座標系 = Measured 座標系**，不需要額外轉換。

| ANSYS 軸 | 方向 | 對應 pole |
|----------|------|----------|
| +X | P1 (0°) | Pair I 正方向 |
| +Y | P5 投影 (60°) | |
| +Z | 垂直向上 | upper tip 側 |

SolidWorks 開 IGES 後方向直接對應。

## 機構約束

| 參數 | 值 | 來源 |
|------|-----|------|
| TILT_UP | 35° | STEP 實測（upper pole 相對水平向上傾斜） |
| TILT_DN | 5.71° | STEP 實測（lower pole 相對水平向下傾斜） |
| TILT_UP + TILT_DN | 40.71° | 不滿足 70.53° 正交條件，需 pole 長度微調補償 |

正交條件要求 TILT_UP + TILT_DN = 70.53°（= 2 × magic angle）。
Hung 實際 40.71° < 70.53°，透過 pole 長度微調（見 `pole-geometry.md`）達到正交。
