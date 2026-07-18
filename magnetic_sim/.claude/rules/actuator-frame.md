# 統一用 actuator 座標系（強制讀取）

**使用者拍板（2026-07-07）**：`matlab/long2016_hexapole_halfcut/Calibration_using_FEM_modeling/` 這整包校正專案，**一律在 actuator（磁極軸）座標系**做「載入 → 擬合 → 呈現」。**新寫 / 改寫任何 pipeline 都必須承接 actuator frame，不可停留在 measure / WP（ANSYS 原始）frame。**

當工作涉及：
- 寫 / 改任何 `load_*.m`、`main*.m`、`fit_*.m`、`solve_d.m`、`build_S.m`、`Pc*`/`ell` 相關腳本
- 在這包新增子流程 / 新變體 / 新 sensor 分析
- 討論「座標系 / measure vs actuator / R_act / 磁極軸」

**動手前先讀完此規則。**

對應 memory：`feedback_actuator_frame_convention.md`
相關規則：`charge-model-source-convention.md`（source/sign 慣例，**跟座標系是兩回事**）、`simulation-constraints.md`（alpha=54.74° / magic angle）
相關 memory：[[nonorthogonal-basis-projection]]、[[long2016-nofixl-bias]]、[[long2016-singlepole-plusx-geom]]

---

## 🔒 核心慣例

**measure / WP frame = ANSYS 原始座標（標準基底，資料 `.dat` 的 (x,y,z)）。**
**actuator frame = 三根磁極軸為基底的座標系。**

### 六極：用 `R_act` 旋轉

```matlab
tip   = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];   % 3x6 極尖端(measure/WP)
dhat  = tip ./ vecnorm(tip);                                       % 3x6 magic-angle 單位向量
R_act = [dhat(:,1), dhat(:,3), dhat(:,5)].';                       % 3x3：P1,P3,P5 單位向量當「列」
P_act = ( R_act * P_meas.' ).';                                    % 位置轉 actuator
B_act = ( R_act * B_meas.' ).';                                    % 磁場轉 actuator（norm 保持）
```

- `R_act` 是**正交旋轉**（magic angle 54.74° 讓 P1/P3/P5 兩兩正交 → `det(R_act)=1`、`R_act.'*R_act=I`）。
- `R_act * v` 每個分量 = 該極單位向量 · v（＝把 measure 向量對三根極軸做**內積投影**）。
- P1 的 actuator 軸（第 1 列）= `[cos35.26°, 0, −sin35.26°] ≈ [0.816, 0, −0.577]`（往下傾 35.26° 的 **3D** 向量，**不是** `[cos(-35.26) 0 0]`）；P5 同角度、z 轉正。
- 轉進 actuator 後，理想電荷位置落在 `Pc_base = R_act*dhat = [+u −u +v −v +w −w]`（正交格）。
- ⚠ magic-angle **理想**方向（離水平 35.26°）≠ 物理上極實際傾角 36.59°（`c.upper_incline` / `c.pole_axis`，別處用）；frame 轉換用的是前者。

### 單極：用 +x 建模（不旋轉）

`main_singlepole.m` / `load_singlepole.m`（tipcut sweep 那條）**不做 `R_act` 旋轉**——因為 mesh 幾何本身就把極軸建在 **+x**、WP 放原點：

```matlab
cnst.SPH_OFST = 0;      % WP 在原點
d_act = [1; 0; 0];      % 極軸 = +x
Pa = Pw;                % 資料一進來就在 actuator frame，無需旋轉
```

→ 六極靠 `R_act` 轉、單極靠 +x 建模，**殊途同歸**，最後都在 actuator frame。

---

## 現況（兩專案全 actuator，已驗；2026-07-15 四子夾合成 current_base/voltage_base）

| 專案 | 載入函式 | 進 actuator 的方式 |
|---|---|---|
| `current_base`（六極，USE_BIAS 統一 fix/no_fix） | `load_coils_actuator.m` | `R_act` 旋轉 P、B |
| `current_base`（單極） | `load_singlepole.m` | +x 建模、不旋轉 |
| `voltage_base`（USE_BIAS 統一） | 沿用 current_base `load_coils_actuator.m` | 同上 |

> 註：舊 `fix_dir` 的 `load_coils.m`（measure-frame per-coil）已隨 Varpar 統一淘汰；current_base 六極一律 `load_coils_actuator`。

---

## 強制規則

1. 這包**所有** pipeline 的載入 / 擬合 / 呈現一律在 **actuator frame**；不可停在 measure/WP frame。
2. 六極用 `R_act = [dhat(:,1),dhat(:,3),dhat(:,5)].'` 旋轉（P + B 都轉）；單極用 +x 建模（`SPH_OFST=0`, `d_act=[1,0,0]`）。
3. `R_act` 必為正交旋轉（`det≈1`）——若 `assert(det(R_act)-1<1e-9)` 失敗＝極軸選錯 / magic angle 沒守，停下查。
4. **座標系（actuator）跟 source/sign 慣例（flip-sink vs 全域 −B）是兩回事**，不可混為一談（sign 見 `charge-model-source-convention.md`）。
5. 新增子流程若讀別人的 `.mat` / 場，確認對方也是 actuator frame（本包預設都是）。

---

## 觸發片語

- 「measure 轉 actuator」/「R_act」/「磁極軸座標系」/「這包在哪個座標系」
- 在 `Calibration_using_FEM_modeling/` 下寫新 load/fit/主程式時

## 何時不適用

- 純 APDL 幾何 / mesh / 抽 `.dat`（那是 measure/WP raw，本就在 ANSYS 原始座標——**轉 actuator 是載入函式的事**）。
- 其他 topic（kuo / zhang quadrupole）有自己的 frame 慣例（見各自 mt_constants）。
