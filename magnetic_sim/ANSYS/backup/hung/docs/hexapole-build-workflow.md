# Hexapole Pole Tip 建模流程

## 步驟 1：決定 tip 位置

### 請使用者提供以下資訊

> **問題 1：球體半徑 R 是多少 mm？**
>
> 這是每個 pole tip 到 working point（球心）的距離。
> 同一 pair 的 tip-to-tip 距離 = 2×R。
>
> 例如：R = 0.5 mm → tip-to-tip = 1.0 mm

> **問題 2：相鄰 pole 的間距角度是多少度？**
>
> 從俯視圖看，相鄰兩根 pole 之間的方位角間距。
> 標準 hexapole = 60°（6 根 pole 均勻分布 360°）。
>
> 例如：spacing = 60°

### 輸入參數總表

| # | 參數 | 說明 | 範例 |
|---|------|------|------|
| 1 | 球體半徑 R (mm) | 每個 tip 到 working point 的距離 | 0.5 |
| 2 | Pole 間距角度 (°) | 相鄰 pole 在俯視圖的方位角間距 | 60 |

### Pole 編號與排列（俯視圖，逆時針）

X+ 方向固定指向 P1。

```
              P5 (upper, 60°)
             /
P3 (lower, 120°)       P1 (lower, 0°) → X+
             \         /
              ● center
             /         \
P2 (upper, 180°)       P4 (upper, 300°)
             \
              P6 (lower, 240°)
```

| 方位角 | Pole | 層 | Pair |
|--------|------|----|------|
| 0° | P1 | lower | I |
| 60° | P5 | upper | III |
| 120° | P3 | lower | II |
| 180° | P2 | upper | I |
| 240° | P6 | lower | III |
| 300° | P4 | upper | II |

對面 pole 相差 180°，同一 pair：P1↔P2、P3↔P4、P6↔P5。
上下層交替排列，每隔 60°。

### 座標系對齊

**ANSYS 座標系 = Measured 座標系**，不需要額外轉換：

| 軸 | ANSYS 方向 | Measured 方向 | 指向 |
|----|-----------|--------------|------|
| +X | +X | X_m | P1 (0°) |
| +Y | +Y | Y_m 投影 | P5 方向 (60°) |
| +Z | +Z | 垂直向上 | 上層 tip 側 |

SolidWorks 開 IGES 後看到的方向直接對應 measured 座標系：
- **+X 方向的 pole = P1**
- **右上 60° 方向 = P5**
- **左邊 180° = P2**

### 執行方式

使用者回答問題後，Claude 自動將值代入 MATLAB 程式碼中的 `[USER]` 變數並執行。

### 輸出

2 張 MATLAB 圖：

1. **3D 球面圖**（`step1_sphere_3D.png`）：6 個 tip 點在球面上的位置，標注 P1-P6 編號與 upper/lower
2. **俯視圖**（`step1_topview.png`）：X+ 指向 P1，標注編號、間距角度、upper/lower

### MATLAB 程式碼

```matlab
%% ===== Step 1: Pole Tip Positions =====
% ===== [USER] 使用者修改這兩行 =====
R = 0.5;          % [USER] 球體半徑 (mm)
spacing = 60;     % [USER] pole 間距角度 (°)

% Magic angle: cos^2(theta) = 1/3
theta_deg = acosd(1/sqrt(3));  % 54.74° (極角)
h = R*cosd(theta_deg);         % 上層平面高度
r = R*sind(theta_deg);         % 上層圓半徑

% Pole 定義 (X+ 指向 P1, 逆時針排列)
pole_names = {'P1','P5','P3','P2','P6','P4'};
pole_types = {'lower','upper','lower','upper','lower','upper'};
pole_az = [0, 60, 120, 180, 240, 300];

% Tip 座標
tips = zeros(6,3);
for i = 1:6
    phi = pole_az(i);
    if strcmp(pole_types{i},'upper')
        tips(i,:) = [r*cosd(phi), r*sind(phi), h];
    else
        tips(i,:) = [r*cosd(phi), r*sind(phi), -h];
    end
end
pair_idx = [1 4; 3 6; 5 2]; % P1↔P2, P3↔P4, P6↔P5

up_clr = [0 0.4 0.8]; lo_clr = [0.8 0.2 0.2];

%% 圖 1：3D 球面圖
fig1 = figure('Color','w','Position',[50 50 800 700]);
[xs,ys,zs] = sphere(40);
surf(xs*R, ys*R, zs*R, 'FaceAlpha',0.1, 'EdgeAlpha',0.06, 'FaceColor',[0.7 0.8 1]);
hold on; axis equal; grid on;
tt = linspace(0, 2*pi, 200);
plot3(r*cos(tt), r*sin(tt), h*ones(size(tt)), 'b--', 'LineWidth',1);
plot3(r*cos(tt), r*sin(tt), -h*ones(size(tt)), 'r--', 'LineWidth',1);
for i = 1:6
    if strcmp(pole_types{i},'upper')
        clr = up_clr; mk = 'v';
    else
        clr = lo_clr; mk = '^';
    end
    plot3(tips(i,1), tips(i,2), tips(i,3), mk, ...
        'MarkerSize',12, 'MarkerFaceColor',clr, 'Color',clr);
    text(tips(i,1)*1.5, tips(i,2)*1.5, tips(i,3)*1.2, ...
        sprintf('%s (%s)', pole_names{i}, pole_types{i}), ...
        'FontSize',11, 'FontWeight','bold', 'Color',clr, ...
        'HorizontalAlignment','center');
end
for k = 1:3
    i1 = pair_idx(k,1); i2 = pair_idx(k,2);
    plot3([tips(i1,1) tips(i2,1)], [tips(i1,2) tips(i2,2)], ...
        [tips(i1,3) tips(i2,3)], 'k--', 'LineWidth',1.2);
end
plot3(0,0,0, 'k+', 'MarkerSize',15, 'LineWidth',2.5);
xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
title(sprintf('Step 1: Pole Tip Positions (R = %.1f mm)', R));
view(35, 25);
print(fig1, 'step1_sphere_3D.png', '-dpng', '-r200');

%% 圖 2：俯視圖 (Top View) — Measured 座標系
fig2 = figure('Color','w','Position',[50 50 1000 1000]);
hold on; axis equal; axis off;
tt = linspace(0, 2*pi, 200);
line_len = R*2.0;
for i = 1:6
    phi = pole_az(i);
    if strcmp(pole_types{i},'upper')
        clr = up_clr;
    else
        clr = lo_clr;
    end
    plot([0, line_len*cosd(phi)], [0, line_len*sind(phi)], '-', ...
        'Color', clr, 'LineWidth', 4);
    plot(tips(i,1), tips(i,2), 'o', 'MarkerSize', 10, ...
        'MarkerFaceColor', clr, 'Color', 'k', 'LineWidth', 1.5);
    text(R*2.7*cosd(phi), R*2.7*sind(phi), ...
        sprintf('%s\n(%s)', pole_names{i}, pole_types{i}), ...
        'FontSize', 14, 'FontWeight', 'bold', 'Color', clr, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end
% Pair 虛線
for k = 1:3
    i1 = pair_idx(k,1); i2 = pair_idx(k,2);
    p1 = line_len*[cosd(pole_az(i1)), sind(pole_az(i1))];
    p2 = line_len*[cosd(pole_az(i2)), sind(pole_az(i2))];
    plot([p1(1) p2(1)], [p1(2) p2(2)], 'k--', 'LineWidth', 1);
end
% 間距弧線
arc_r = R*1.5;
for i = 1:6
    a1 = pole_az(i);
    if i < 6
        a2 = pole_az(i+1);
    else
        a2 = pole_az(1) + 360;
    end
    arc_a = linspace(a1, a2, 30);
    plot(arc_r*cosd(arc_a), arc_r*sind(arc_a), '-', 'Color', [0.5 0.5 0.5]);
    mid_a = (a1 + a2)/2;
    text(arc_r*1.12*cosd(mid_a), arc_r*1.12*sind(mid_a), ...
        sprintf('%d°', spacing), 'FontSize', 11, 'Color', [0.4 0.4 0.4], ...
        'HorizontalAlignment', 'center');
end
% Measured 座標系：X_m → P1 (0°), Y_m → P5 (60°)
ax_len = R*0.85;
red = [0.85 0 0];
plot([0 ax_len], [0 0], '-', 'Color', red, 'LineWidth', 2.5);
plot(ax_len, 0, '>', 'Color', red, 'MarkerFaceColor', red, 'MarkerSize', 8);
text(ax_len+0.05, -0.02, 'X_m', 'FontSize', 13, 'FontWeight', 'bold', 'Color', red);
plot([0 ax_len*cosd(60)], [0 ax_len*sind(60)], '-', 'Color', red, 'LineWidth', 2.5);
plot(ax_len*cosd(60), ax_len*sind(60), '^', 'Color', red, 'MarkerFaceColor', red, 'MarkerSize', 8);
text(ax_len*cosd(60)+0.04, ax_len*sind(60)+0.02, 'Y_m', ...
    'FontSize', 13, 'FontWeight', 'bold', 'Color', red);
plot(0, 0, 'k.', 'MarkerSize', 10);
% 圖例
lx0 = -R*2.8; ly0 = -R*2.5;
plot(lx0, ly0, 'o', 'MarkerSize', 8, 'MarkerFaceColor', up_clr, 'Color', 'k');
text(lx0+0.05, ly0, 'upper pole', 'FontSize', 12, 'Color', up_clr, ...
    'FontWeight', 'bold', 'VerticalAlignment', 'middle');
plot(lx0, ly0-0.12, 'o', 'MarkerSize', 8, 'MarkerFaceColor', lo_clr, 'Color', 'k');
text(lx0+0.05, ly0-0.12, 'lower pole', 'FontSize', 12, 'Color', lo_clr, ...
    'FontWeight', 'bold', 'VerticalAlignment', 'middle');
plot([lx0 lx0+0.06], [ly0-0.24 ly0-0.24], 'k--', 'LineWidth', 1);
text(lx0+0.08, ly0-0.24, 'pair connecting line', 'FontSize', 10, ...
    'Color', [0.3 0.3 0.3], 'VerticalAlignment', 'middle');
title(sprintf('Top View (Measured coord.): R = %.1f mm, spacing = %d°', R, spacing), 'FontSize', 14);
margin = R*3.2;
xlim([-margin margin]); ylim([-margin margin]);
print(fig2, 'step1_topview.png', '-dpng', '-r200');
```

---

## 步驟 2：擺放 pole 並建立 APDL 模型

### 請使用者提供以下資訊

> **問題：你的上層和下層 pole 要以什麼角度擺放？**
>
> 請參考下圖，回答兩個角度：

```
        Block (上, 外)
            \
             \  α_upper（上層傾斜角）
              \_________ 水平面
               \
                Tip (下, 內)          ← upper pole (P2, P4, P5)
                                       從 tip 沿軸往上、往外延伸到 block

                Tip (上, 內)          ← lower pole (P1, P3, P6)
               /  α_lower（下層傾斜角）
              /_________ 水平面
             /
        Block (下, 外)                 從 tip 沿軸往下、往外延伸到 block
```

| # | 參數 | 說明 | 範例 |
|---|------|------|------|
| 3 | α_upper (°) | 上層 pole 相對水平往上的傾斜角 | 35° |
| 4 | α_lower (°) | 下層 pole 相對水平往下的傾斜角 | 5.71° |

> **注意**：兩個角度都填正值。上層代表「往上幾度」，下層代表「往下幾度」。

### 執行方式

使用者回答問題後，Claude 自動將值代入 APDL 程式碼中的 `[USER]` 變數、執行 ANSYS、匯出 IGES。

### Pole 尺寸（from Mag_Pole_Bottom.STEP）

| 參數 | 值 | 說明 |
|------|-----|------|
| 圓柱半徑 | 3.175 mm | |
| 錐段長度 | 15.875 mm | semi-angle 11.31° |
| 總長 | 43 mm | |
| Tip | 尖錐頂點 | STEP 無 fillet |

### 輸出

- IGES 檔案 → SolidWorks 開啟量測驗證
- 使用者在 SolidWorks 量 tip-to-tip 距離，應為 2×R（例如 1.0 mm）

### APDL 注意事項

| # | 陷阱 | 解法 |
|---|------|------|
| 1 | 乘號 `*` 前後有空格 → APDL 當成註解 | `cos(theta)*R_sphere`（不留空格） |
| 2 | IGES 單位：ANSYS 寫 meters，SolidWorks 讀成 inches | 內部用 `MM = 1.0/25.4` 轉英寸 |
| 3 | Magic angle 搞混極角/仰角 | 極角 θ = acosd(1/√3) = 54.74° |
| 4 | VROTAT 接縫與削平面重合 | 截面偏移用 profile45 旋轉 45° |

### APDL 程式碼（已測試通過）

使用者需修改的參數標記為 `[USER]`。
含座標對齊（ANSYS +X = P1）、兩段式削平（前段 D 型 + 後段全圓）。

- **上層 pole**：削平面朝上
- **下層 pole**：削平面朝下（轉 180°）

完整腳本見 `magnetic_sim/ANSYS/backup/hung/apdl/variants/MT_Hung_SphereModel.txt`，以該檔案為準。

### 削平方法

不使用 VSBV（已知會失敗），改用**兩段式 VROTAT**：

| 段 | 範圍 | VROTAT | 結果 |
|----|------|--------|------|
| 前段 | tip → 28mm | 180°, 2 divisions | D 型半圓（削平） |
| 後段 | 28mm → 43mm | 360°, 4 divisions | 全圓柱 |

削平面方向由截面偏移決定：
- Upper pole: `-side_perp` → 保留下半 → 削平面朝上
- Lower pole: `+side_perp` → 保留上半 → 削平面朝下

### SolidWorks 驗證清單（步驟 2）

- [ ] P1 在 +X 方向（SolidWorks 右邊）
- [ ] tip-to-tip 距離 = 2×R mm
- [ ] 俯視 60° 等距
- [ ] 前段 D 型、後段全圓柱
- [ ] 上層削平面朝上、下層朝下

---

## 步驟 3：安裝 Block

### 請使用者提供以下資訊

> **問題 5：上層 pole 伸入 block 多少 mm？（沿 pole 軸量）**
> 例如：7mm
>
> **問題 6：下層 pole 伸入 block 多少 mm？（沿 pole 軸量）**
> 例如：4.5mm

### 輸入參數

| # | 參數 | 說明 | 範例 |
|---|------|------|------|
| 5 | 上層伸入量 (mm) | pole 尾端沿軸進入 block 的深度 | 7 |
| 6 | 下層伸入量 (mm) | 同上 | 4.5 |

### 執行方式

使用者回答後，Claude 自動代入 APDL 並執行。

### Block 尺寸（from STEP）

| | Upper Block | Lower Block |
|--|-------------|-------------|
| 外形 | 25 × 22 × 10 mm | 22 × 22 × 10 mm |
| L 型截面 Part 1 | side -0.25~12.75, up ±11 | side 2~9, up ±11 |
| L 型截面 Part 2 | side -12.25~-0.25, up ±5 | side -13~2, up ±5 |
| STEP 來源 | `upper_block.STEP` | `lower_block.STEP` |
| Setting 參考 | `upper_block_setting.STEP` | `lower_block_setting.STEP` |

### 上下層安裝邏輯不同

```
Upper pole (35°傾斜):          Lower pole (5.71°傾斜):
pole 從 block 底面進入           pole 從 block 側面穿過

    ┌──────────┐ outer          ┌──────────┐
    │          │                │  ●═══════╪═ pole (幾乎水平)
    │  ● end   │                │  pole 在  │
    │  ╱       │                │  block    │
    │ ╱ 7mm    │                │  中央     │
    ├╱─────────┤ entry          └──────────┘
    ╱ pole

Block Z:                        Block Z:
  bottom = entry_z              center = pole_tail_z
  top = entry_z + 10mm          ±5mm (居中)

Block XY 中心:                  Block XY 中心:
  = pole 軸在 outer face         = pole 尾端 XY
```

### Block 建模方法

- **不使用 LOCAL CS**（已知會造成位置錯誤）
- 用 **global 座標算 8 角點 + V 指令**建每個 sub-block
- L 型截面 = 2 個矩形 VADD 合併
- 寬邊朝 workspace（XY 方向用 azimuth+180° 旋轉）

### APDL 注意事項（Block 專用）

| # | 陷阱 | 解法 |
|---|------|------|
| 1 | CSKP + LOCAL CS 建 BLOCK 位置會錯 | 用 global 座標 V 指令 |
| 2 | Upper/Lower 安裝邏輯不同 | Upper: 從底面進入; Lower: 從側面穿過，居中 |
| 3 | STEP 截面 slot 不在原點 | Upper 平移 +0.25mm, Lower 平移 -2mm |
| 4 | VADD 可能失敗（不影響幾何） | 45 volumes 正常（36 pole + 9 block） |

### SolidWorks 驗證清單（步驟 3）

- [ ] 6 個 block 都在 pole 尾端
- [ ] 上層 block：pole 從底面進入，pole 尾端在 block 內
- [ ] 下層 block：pole 從側面穿過，pole 在 block 垂直中央，不凸出
- [ ] Block 寬邊朝 workspace（朝 hexapole 中心方向）

---

## 步驟 4：安裝 Yoke（環形鐵板）

### 參數

| 參數 | 值 | 說明 |
|------|-----|------|
| YOKE_RI | 38.0 mm | 內徑（from 導磁鐵環.STEP） |
| YOKE_RO | 62.5 mm | 外徑（from 導磁鐵環.STEP） |
| YOKE_T | 2.0 mm | 厚度 |
| COIL_H | 14.0 mm | Block 頂到 Yoke 底的間距（coil 空間） |

### Z 位置鏈

```
fc_up_endz = fc_tz(P2) + POLE_TOTAL_LEN*sin(TILT_UP)     ! 上層 pole 尾端 Z
fc_up_blk_top = fc_up_endz + BLK_T/2                       ! 上層 block 頂
fc_yoke_zbot = fc_up_blk_top + COIL_H                      ! Yoke 底面
fc_yoke_ztop = fc_yoke_zbot + YOKE_T                       ! Yoke 頂面
```

### 建模方法

截面矩形 KP（R_inner, zbot）→（R_outer, zbot）→（R_outer, ztop）→（R_inner, ztop），加軸 KP 在 Z 軸上，VROTAT 360°, 4 divisions。

---

## 步驟 5：安裝 Guide Post（3 根）+ Upper Core（3 根）

### Guide Post（下層 P1, P3, P6）

| 參數 | 值 | 說明 |
|------|-----|------|
| GP_R | 4.0 mm | 半徑 |
| GP 頂端 Z | fc_yoke_zbot | 接 yoke 底面 |
| GP 底端 Z | fc_lo_endz - BLK_T/2 + 3.56 mm | lower block 底 + 3.56mm |
| GP XY | pole 尾端 XY | = fc_tx(i) + POLE_TOTAL_LEN*cos(TILT_DN)*cos(azim) |

建法：截面矩形 KP + VROTAT 360°, 4 divisions。3 根迴圈建立（P1, P3, P6）。

### Upper Core（上層 P2, P4, P5）

| 參數 | 值 | 說明 |
|------|-----|------|
| UC_R | 4.0 mm | 與 GP 相同 |
| UC 底端 Z | fc_up_blk_top2 = fc_up_endz + BLK_T/2 | upper block 頂面 |
| UC 頂端 Z | fc_yoke_zbot | 接 yoke 底面 |
| UC XY | pole 尾端 XY | = fc_tx(i) + POLE_TOTAL_LEN*cos(TILT_UP)*cos(azim) |

建法同 GP。3 根迴圈建立（P2, P4, P5）。

---

## 步驟 6：VADD + 空氣域 + VOVLAP + VGLUE

### 合併鋼件

```apdl
VADD, P51X    ! 所有 steel volumes → 1 個合併體
```

### 空氣域

| 域 | 幾何 | 參數 |
|----|------|------|
| 球形精細區 | SPHERE | R = 7.0 mm，原點 |
| 圓柱外空氣 | CYL4 | R = 80 mm, H = 80 mm，Z 中心在模型中點 |

### Boolean 操作

```apdl
VOVLAP, ALL         ! 重疊所有 volumes
NUMCMP, VOLU         ! 壓縮 volume 編號
VGLUE, 8, 1,2,3,4,5,6  ! Glue 鋼件(V8)與 tips(V1-6)
```

### 結果：9 Volumes

| Volume | 內容 | 材料 |
|--------|------|------|
| V1-V6 | 6 根 pole tip（球內鐵） | MAT 2 (μr=280) |
| V7 | 球內空氣 | MAT 1 (μr=1) |
| V8 | 鋼主體（yoke + GP + cores + blocks + pole 外段） | MAT 2 (μr=280) |
| V9 | 圓柱外空氣 | MAT 3 (μr=1) |

---

## 步驟 7：材料 + Mesh + SOURC36 + 求解

### 材料與元素

```apdl
ET, 1, SOLID96       ! 磁純量勢元素
ET, 2, SOURC36       ! 電流源元素
MP, MURX, 1, 1       ! V7: 球內空氣
MP, MURX, 2, 280     ! V1-6, V8: 1018 鋼
MP, MURX, 3, 1       ! V9: 外空氣
```

### 分區 Mesh

| 區域 | Volumes | ESIZE | 說明 |
|------|---------|-------|------|
| Pole tips | V1-V6 | 0.3 mm | 最細（尖端需高解析） |
| 球內空氣 | V7 | 0.3 mm | 細（WP 區域） |
| 鋼主體 | V8 | 1.5 mm | 中等 |
| 外空氣 | V9 | 4.0 mm | 粗 |

目標：~2-3M elements。

### SOURC36 Coil 定義

| 參數 | 值 | 說明 |
|------|-----|------|
| TURNS | 70 | 匝數 |
| COIL_MEAN_R | 11.0 mm | (R_in=10 + R_out=12) / 2 |
| COIL_DY | 2.0 mm | R_out - R_in = 12 - 10 |
| COIL_DZ | 15.0 mm | 高度 |
| Coil 位置 Z | block_top + COIL_DZ/2 | 貼在 block 正上方 |
| Coil XY | pole 尾端 XY | 跟 GP/UC 相同 |
| 繞向 | 順時針（N1/N2 swap） | flux 向 -Z（toward block） |

Real constants：`R, i, 1, TURNS*CURR_ARRAY(i), COIL_DY, COIL_DZ`

CURR_ARRAY：每個檔案只激發 1 個 coil（coilN 檔 → CURR_ARRAY(N)=1，其餘=0）。

### 邊界條件

```apdl
! 圓柱頂面 MAG=0
NSEL, S, LOC, Z, z_top ± 1e-5
D, ALL, MAG, 0
! 圓柱底面 MAG=0
NSEL, S, LOC, Z, z_bot ± 1e-5
D, ALL, MAG, 0
! 圓柱側面 MAG=0 (用 LOCAL cylindrical CS)
LOCAL, 99, 1, 0, 0, 0
NSEL, S, LOC, X, AIR_CYL_R*0.99, AIR_CYL_R*1.01
D, ALL, MAG, 0
```

### 求解

```apdl
magsolv, 3, , , , , 1    ! DSP 方法
```

### SolidWorks 驗證清單（步驟 7 前）

- [ ] 9 volumes 正確（V1-6 tips, V7 sphere air, V8 steel, V9 cyl air）
- [ ] 材料指定正確（VATT）
- [ ] Mesh 完成無 error
- [ ] SOURC36 nodes 在正確位置（block 上方）
- [ ] MAG=0 BC 在外表面
