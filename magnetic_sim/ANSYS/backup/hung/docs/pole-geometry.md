# Hung Hexapole — Pole Geometry

## Pole 編號與配對

**下層 poles:** 1, 3, 6
**上層 poles:** 2, 4, 5

| Pair | Lower Pole | Upper Pole |
|------|-----------|------------|
| I | 1 | 2 |
| II | 3 | 4 |
| III | 6 | 5 |

## Pole 傾斜角度（STEP 實測）

### Hung

| Group | Poles | Tilt from Horizontal | Direction |
|-------|-------|---------------------|-----------|
| Upper | 2, 4, 5 | **35°** | 向下傾斜（tip 在 block 下方） |
| Lower | 1, 3, 6 | **5.71°** | 向上傾斜（tip 在 block 上方） |

```
        Block (上, 外)
            \
             \  35° from horizontal
              \_________ 水平面
               \
                Tip (下, 內)          ← upper pole (2, 4, 5)

                Tip (上, 內)          ← lower pole (1, 3, 6)
               /  5.71°
              /_________ 水平面
             /
        Block (下, 外)
```

上下層 tip 從兩側匯聚到 working point。

**Source:** STEP 實測 `magnetic_sim/ANSYS/backup/hung/step_for_fem/組合件1-7.STEP`

### Long2016

| Group | Poles | Tilt from Horizontal | Direction |
|-------|-------|---------------------|-----------|
| Lower 3 | P1(0°), P3(120°), P6(240°) | **0°** (水平) | 完全水平 |
| Upper 3 | P5(60°), P2(180°), P4(300°) | **35°** | 向下傾斜（tip 在 block 下方） |

> 上層 35° 為 STEP 實測值（APDL 參數化模型計算為 36.6°，有 ~1.6° 差異）

```
        Block (上, 外)
            \
             \  35° from horizontal
              \_________ 水平面
               \
                Tip (下, 內)

    ─────────────────────── lower pole level
    Block (外) ────────── Tip (內)   ← 完全水平 (0°)
```

**Source:** STEP 實測 `01Long Fei Tweezer/01CAD Model/hall sensor assembly and Supporter Structure (2014-10-9)/ASSEM_The_Whole_Setup.STEP`

### Comparison

| | Long2016 | Hung |
|--|----------|------|
| Upper tilt | 35° 向下 | 35° 向下 |
| Lower tilt | 0° (水平) | 5.71° **向上** |

上層角度相同（35°），下層 Hung 比 Long 多了 5.71° 向上傾斜。

---

## Tip Apex 位置（from Hexapole_Assembly_FEM.STEP）

| Pole | Type | Apex (mm) |
|------|------|-----------|
| 1 | lower | (-41.613, -82.923, -18.607) |
| 2 | lower | (-41.304, -83.321, -18.841) |
| 3 | lower | (-41.865, -83.276, -18.731) |
| 4 | upper | (-41.864, -82.940, -18.105) |
| 5 | upper | (-41.274, -82.883, -18.027) |
| 6 | upper | (-41.619, -83.654, -17.925) |

Hex center: (-42.21, -83.10) mm
全部 6 個 tip 匯聚於 hex center 附近（r < 1mm）。

---

## Pair Tip-to-Tip 連線夾角

### 連線方向

| Pair | 連線 | Direction (normalized) |
|------|------|----------------------|
| I | Pole 1 → Pole 2 | (0.309, -0.398, -0.234) |
| II | Pole 3 → Pole 4 | (0.001, 0.336, 0.626) |
| III | Pole 6 → Pole 5 | (0.345, 0.771, -0.102) |

### 3D 夾角

| Pair 組合 | 銳角 | 鈍角 |
|-----------|------|------|
| I(1↔2) vs II(3↔4) | 77.0° | 103.0° |
| I(1↔2) vs III(6↔5) | 79.0° | **101.0°** |
| II(3↔4) vs III(6↔5) | 70.3° | 109.7° |

**三對 pole pair 的 tip-to-tip 連線在 3D 空間中不正交（≠90°）。**

鈍角範圍：101° ~ 110°，偏離正交 11° ~ 20°。

### 實測驗證

Pair I ↔ III 鈍角：
- 計算值：101.0°
- STEP 實測值：**100.98°** ✓

---

## 設計約束條件

| # | 約束 | 說明 |
|---|------|------|
| 1 | **3 對 pole pair 的 tip-to-tip 連線互相正交** | 3 條連線在 3D 空間中兩兩夾 90° |
| 2 | **每對 tip-to-tip 直線距離 = 1 mm** | opposing pair 的兩個 tip 端點間距 1 mm |

---

## Magic Angle 正交條件推導

### 座標系

- 原點 O：working point（6 tip 匯聚處）
- Z 軸：垂直向上
- XY 平面：水平面
- alpha：pair 軸與水平面的仰角
- phi：XY 平面上的方位角

### 設定

3 對 pair 軸以 120° 方位角排列，每對與水平面夾角都是 alpha。方向向量：

```
v1 = (cos(alpha)cos(phi1),  cos(alpha)sin(phi1),  sin(alpha))
v2 = (cos(alpha)cos(phi2),  cos(alpha)sin(phi2),  sin(alpha))
v3 = (cos(alpha)cos(phi3),  cos(alpha)sin(phi3),  sin(alpha))
```

其中 phi2 = phi1 + 120°，phi3 = phi1 + 240°。

### 正交條件

v1 · v2 = 0：

```
v1 · v2 = cos^2(alpha) × [cos(phi1)cos(phi2) + sin(phi1)sin(phi2)] + sin^2(alpha)
```

中括號項 = cos(phi2 - phi1) = cos(120°) = -1/2，所以：

```
v1 · v2 = cos^2(alpha) × (-1/2) + sin^2(alpha)
```

代入 sin^2(alpha) = 1 - cos^2(alpha)：

```
= -cos^2(alpha)/2 + 1 - cos^2(alpha)
= 1 - 3cos^2(alpha)/2
```

令 = 0：

```
1 - 3cos^2(alpha)/2 = 0
3cos^2(alpha)/2 = 1
cos^2(alpha) = 2/3
alpha = arccos(sqrt(2/3)) = 35.26°
```

### 結論

當 3 對 pair 以 **120° 等距方位排列**且**傾斜角相同**時，alpha = **35.26°（magic angle）**是唯一讓 3 軸正交的角度。

Long2016 設計選擇 ~35° 即為逼近此 magic angle。

### 正交時的 tip-to-tip 分量（d = 1mm 端到端距離）

```
水平分量 h  = d × cos(35.26°) = d × sqrt(2/3) = 0.8165 mm
垂直分量 dz = d × sin(35.26°) = d / sqrt(3)   = 0.5774 mm
```

---

## 維持正交的上下層角度組合

### 公式

正交條件要求 pair 的 tip-to-tip 向量滿足 dz/dh = 1/√2。

當上層傾斜 a_upper（向下）、下層傾斜 a_lower（向上）時：

```
(sin(a_upper) + sin(a_lower)) / (cos(a_upper) + cos(a_lower)) = 1/√2
```

簡化為：

```
a_lower = 70.53° - a_upper
```

（即上下層角度之和 = 70.53° = 2 × 35.26°）

### 可行組合

| 方案 | Upper (↓) | Lower (↑) | 合計 | 正交? | 備註 |
|------|----------|----------|------|-------|------|
| A 完美對稱 | 35.26° | 35.27° | 70.53° | ✓ | 最理想（場對稱） |
| B | 35° | 35.53° | 70.53° | ✓ | 接近對稱 |
| C | 40° | 30.53° | 70.53° | ✓ | |
| D | 30° | 40.53° | 70.53° | ✓ | |
| E | 45° | 25.53° | 70.53° | ✓ | |
| **F Long2016** | **35°** | **0°** | **35°** | **✗** | pair tilt=17.5°≠35.26° |
| **G Hung 現況** | **35°** | **5.71°** | **40.71°** | **✗** | pair tilt≠35.26° |

### 重要發現

**Long2016 的 upper=35° + lower=0° 並不滿足正交條件！**

正交需要 a_upper + a_lower = 70.53°，但 Long2016 只有 35° + 0° = 35°。

Long2016 之所以「近似正交」是因為 R_norm_xy 和 R_norm_z 的設定（tip 位置偏移）另外補償了角度不足。但嚴格來說，pair 軸傾斜只有 17.5°，不是 magic angle 35.26°。

### 對 Hung 的建議

若要從**角度**達到正交（不依賴 tip 位置微調）：
- 維持 upper = 35° → lower 需要 **35.53°** 向上
- 或 upper = 40° → lower = **30.53°** 向上
- 最簡單：兩層都設 **35.26°**（完美對稱）

---

## Pole 長度調整（正交 + 1mm 最佳化結果）

### 方法
保持 pole 傾斜角不變（upper 35°, lower 5.71°），僅調整每根 pole 長度（≤ 1.3% of 43mm = 0.54mm），使 3 對 pair 的 tip-to-tip 連線達到正交 + 1mm。

### 調整量

| STEP Pole | Type | dl (mm) | 動作 | 佔 pole 長度 |
|-----------|------|---------|------|-------------|
| 1 | lower | -0.074 | 縮短 | 0.17% |
| 2 | lower | -0.542 | 縮短 | 1.26% |
| 3 | lower | -0.437 | 縮短 | 1.02% |
| 4 | upper | +0.450 | 加長 | 1.05% |
| 5 | upper | +0.367 | 加長 | 0.85% |
| 6 | upper | +0.060 | 加長 | 0.14% |

### 調整後 tip 位置

| Pole | Type | New Apex (mm) |
|------|------|---------------|
| 1 | lower | (-41.6129, -82.8493, -18.6144) |
| 2 | lower | (-40.8377, -83.5914, -18.8949) |
| 3 | lower | (-42.2417, -83.4925, -18.7745) |
| 4 | upper | (-41.5447, -83.1250, -18.3634) |
| 5 | upper | (-41.5346, -83.0326, -18.2374) |
| 6 | upper | (-41.6189, -83.6049, -17.9594) |

### 驗證（MATLAB 計算）

```
Pair A (1↔6): d = 1.000000 mm ✓
Pair B (2↔4): d = 1.000000 mm ✓
Pair C (3↔5): d = 1.000000 mm ✓

Angle A-B: 90.0000° ✓
Angle A-C: 90.0000° ✓
Angle B-C: 90.0000° ✓
```

### ANSYS 驗證（已完成）

模型建構：`magnetic_sim/ANSYS/backup/hung/apdl/variants/MT_Hung_SphereModel.txt` → `Hung_SphereModel.iges`
SolidWorks 量測確認：

```
6 tips on R=0.5mm sphere: ✓
Tip-to-tip (Pair I,  P1↔P2): 1.000 mm ✓
Tip-to-tip (Pair II, P3↔P4): 1.000 mm ✓
Tip-to-tip (Pair III,P6↔P5): 1.000 mm ✓

Pair I  vs II:  90.00° ✓
Pair I  vs III: 90.00° ✓
Pair II vs III: 90.00° ✓
```

ANSYS 模擬比較（Coil1 激發, 70 A-turns）：

| 量測位置 | Hung | Long2016 |
|----------|------|----------|
| P1 tip (R=0.5mm) | 333 mT | 350 mT |
| WP 幾何中心 | 9.3 mT | 8.7 mT |

兩者在等效位置的 B-field 幾乎相同。

---

## Source
- Hung STEP: `magnetic_sim/ANSYS/backup/hung/step_for_fem/Hexapole_Assembly_FEM.STEP`, `magnetic_sim/ANSYS/backup/hung/step_for_fem/組合件1-7.STEP`
- Long2016 STEP: `01Long Fei Tweezer/01CAD Model/hall sensor assembly and Supporter Structure (2014-10-9)/ASSEM_The_Whole_Setup.STEP`
- Cone apex 提取方法：semi-angle matching (0.1974 rad = 11.31°)
