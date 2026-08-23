# step-to-apdl

把 STEP 拆成 APDL 重建所需的**參數表**(尺寸 + 軸 + 圓柱 + 每個值來自 STEP 哪行)。
**不**自己寫 APDL,那是 `apdl-geom-build.md` 的事。

## 何時用

- 新拿到 STEP 要重建 APDL 幾何
- STEP 改了想 diff 哪些尺寸動了

## 輸入

| 參數 | 範例 |
|---|---|
| `{topic}` | `zhang_quadrupole` |
| `{step}` | `magnetic_sim/ANSYS/main/CAD/{topic}/STEP/<basename>.STEP` |
| `{out}` | `magnetic_sim/ANSYS/main/apdl/{topic}/geom/<basename>_params.md`(預設) |

## 前置

- `cad-export.md` 跑完(STEP 已 export + model-check §1 通過)

## 方法論(一句話)

**3 維交叉檢查**:同一個 3D 圓柱必在 STEP 三類實體裡同時出現
—— `CARTESIAN_POINT`(軸端點)+ `CIRCLE`(底圓)+ `CYLINDRICAL_SURFACE`(側面)。
缺一就是 2D 弧 / 邊界曲線,**不可**當圓柱 import APDL。
(`zhang_quadrupole` 4 次 iteration 學到的教訓,見 memory `feedback_step_geom_extraction`)

## 步驟

1. **讀 header**:確認 `SI_UNIT(.METRE.)` 或 `.MILLI./.METRE.`、AP version(214 / 242)
2. **抓 CARTESIAN_POINT**:
   ```
   grep "CARTESIAN_POINT" {step}
   ```
   列 `(ID, x, y, z)`。
3. **抓 CIRCLE**:
   ```
   grep "CIRCLE" {step}
   ```
   列 `(ID, 中心ID, 半徑, 軸ID)`。
4. **抓 CYLINDRICAL_SURFACE**:
   ```
   grep "CYLINDRICAL_SURFACE" {step}
   ```
   列 `(ID, 軸ID, 半徑)`。
5. **三方交叉**:每個候選圓柱用 ID 串接,三類都有 → `trust=3D ✓`;少一類 → `trust=2D ⚠`
6. **寫參數表 `{out}`**(markdown 表格):
   ```markdown
   | name        | value  | unit | source           | trust |
   |-------------|--------|------|------------------|-------|
   | POLE_R      | 3.0e-3 | m    | line 412 CIRCLE  | 3D ✓  |
   | POLE_L      | 15e-3  | m    | line 388 CYL_SRF | 3D ✓  |
   | inner_arc   | 0.5e-3 | m    | line 244 CIRCLE  | 2D ⚠  |
   ```
7. ⏸ **檢查點 → [`model-check.md` §3](model-check.md#3-參數表-sanity-check)**

## 產物

- [ ] `{out}` 參數表(每行有 source line + trust 標記)
- [ ] `⚠` 行人工裁決後加註 `keep / drop / refine`

## 常見坑

- 2D 弧誤當 3D 圓柱 → 後面 APDL CYL 出來位置 / 半徑都對不上
- header `.MILLI./.METRE.`(=mm)沒看到 → 數值全錯 1000×
- 多 body STEP → 每 body 各有獨立 CARTESIAN_POINT 集,要分群處理
- AP242 + NURBS → 沒 CIRCLE / CYL_SRF 可抓 → 回 SolidWorks 改 export 為 AP214

## 適用 topic

`zhang_quadrupole` ✓(4 次)/ `kuo_quadrupole V2` ✓ / `long_fei` 部分手算可補。
與 pole 配置無關。
