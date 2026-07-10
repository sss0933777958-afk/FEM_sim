# ANSYS-CAD 尺寸對齊規則(強制讀取)

**原則**:**ANSYS 建的模型尺寸必須跟 CAD 檔(SolidWorks STEP / IGES)對齊**。CAD 是 source of truth,ANSYS 不可獨立設定不一致的尺寸。

當任何工作涉及:
- 寫/改 APDL 幾何腳本(`MT_Geom*.txt`、`MT_Sim_*.txt` 內 geometry block)
- 寫/改 `mt_constants.m`(或任何幾何常數檔)
- 跑新 sim 前的幾何 sanity check
- 看到 ANSYS sim 結果跟 CAD 視覺檢查不一致

**動手前必須先讀完此規則全文**。

對應 memory:`feedback_ansys_cad_alignment.md`
相關規則:`.claude/rules/main-workflows.md`(model-check SOP)、`.claude/rules/sim-cleanup.md`
相關 memory:`feedback_step_geom_extraction.md`(STEP 抽尺寸方法論)

---

## 🔒 核心強制規則

### 1. CAD STEP/IGES = source of truth

**順序**:SolidWorks 設計 → STEP/IGES export → APDL 幾何**對齊** CAD 數值

不可:
- 先寫 APDL 用「想當然耳」的數字(例如 1/8" 換算 3.175 mm)
- 假設 mt_constants 對而不驗證
- 「方便」用 round number(3.0 mm)而 CAD 是 3.047 mm

### 2. 改 ANSYS 幾何前必先量 CAD

**Workflow**:

```
   想改某個 dim
        ↓
   開 SolidWorks 量該 dim → 拿到 CAD 值
        ↓
   檢查 mt_constants.m 對應 var → 對比是否一致
        ↓
   若不一致 → ⚠ STOP, 報告使用者「mt_constants X = Y mm vs CAD = Z mm」
        ↓
   等使用者拍板用哪個值 → 才可改 APDL
```

### 3. 新 sim 第一步 = 從 CAD 抽參數表

不要直接抄舊 sim 腳本的數字。新 topic / 新 variant 第一步:
1. 量 CAD(SolidWorks 或 parse STEP)
2. 寫進 `magnetic_sim/ANSYS/main/apdl/<topic>/geom/<variant>_params.md` 參數表
3. APDL 腳本所有數字**只從參數表引用**,不寫死

對應 SOP:`magnetic_sim/ANSYS/main/doc/workflows/step-to-apdl.md`

### 4. 發現不一致時必通報

任何時候(寫 sim、debug、看圖、量 CAD)發現 ANSYS 跟 CAD 數字不對齊:

✅ **必做**:
- 立刻**停下手上工作**,把不一致報告給使用者
- 列出兩邊數值 + 差異百分比 + 物理影響估計
- 問使用者:「哪個是對的?要 (a) 改 ANSYS 對齊 CAD,(b) 改 CAD 對齊 ANSYS,(c) 查 dissertation/spec 確認?」

❌ **不可**:
- 自己選一個值繼續做
- 假設「差 4% 應該還好」直接 ignore
- 改任何一邊而沒先問

---

## STEP 抽尺寸方法論

對應 memory:[[step-geom-extraction]]

### 方法 1:量在 SolidWorks(直觀但要手動)

- 開 STEP / 原始 SLDPRT
- 用 Measure tool 量關鍵 dim
- 記入參數表

### 方法 2:Parse STEP 文字檔(可重複、可自動)

STEP 是純文字。常用 entity:
- `CARTESIAN_POINT` — 所有 vertex 座標(每個 mm)
- `CIRCLE` — analytic 圓(罕用,SolidWorks 多用 B-spline)
- `CYLINDRICAL_SURFACE` — analytic 圓柱(罕用)
- `CONICAL_SURFACE` — analytic 圓錐(罕用)
- `B_SPLINE_SURFACE_WITH_KNOTS` — SolidWorks 預設,沒明顯尺寸

**bash 抽 vertices**:
```bash
F=path/to/file.STEP
grep "CARTESIAN_POINT" "$F" | \
  sed -nE 's/.*\(\s*([+-]?[0-9.E+-]+)\s*,\s*([+-]?[0-9.E+-]+)\s*,\s*([+-]?[0-9.E+-]+)\s*\).*/\1 \2 \3/p' > pts.txt

# 抓特定區域(filter by x range, |y|<tolerance, z range)
awk '$1>X_MIN && $1<X_MAX && abs($2)<TOL && $3>Z_MIN' pts.txt | sort -u
```

### Trap(踩過的坑)

per [[step-geom-extraction]]:
- B-spline 控制點 ≠ surface 實際邊界(R_tip = 0.04 mm 設計值,STEP 抽到 0.033 是 control point 距離,不是 fillet 半徑)
- 邊界 corner vertices 可能是 bounding box 控制點,不在實際 surface 上(zhang_quadrupole 4 次 iteration 教訓)
- 多元件 assembly 不同 part 在同 x 重疊 → filter 必夾 |y| / z 範圍縮到單一 part

---

## 對齊不一致時的 3 種修法

對應記憶:[[ansys-cad-alignment]]

| Path | 動作 | 成本 | 何時用 |
|---|---|---|---|
| **A** | 改 mt_constants.m / APDL 對齊 CAD,re-run sim | 3-5 hr per sim batch | CAD 對的時候 |
| **B** | ANSYS IGESIN import CAD,放棄 primitives | 1-2 週重寫 pipeline | 想長期 CAD-driven |
| **C** | 改 CAD 對齊 ANSYS | 1-2 hr SolidWorks 編輯 | ANSYS 對的時候(罕見)|

**預設 Path A**,除非有明確理由走 B / C。

---

## 例外:何時可不對齊

唯一允許 ANSYS 數值跟 CAD 不同的情境:

1. **μ_r 等效法等「物質屬性」trick**(per [[long-fei-mueq-gap-approach]])— 不動幾何,改 μ_r 模擬氣隙
2. **明顯標記為「ideal spec」的 reference sim**(例如算 paper 對照,用 paper 的英制 round number)— 須在 sim 註解清楚標 `! IDEAL SPEC, not matching CAD`
3. **使用者明確要求**「不對齊也跑,我有別的考量」

以上任一情境,必須在 sim 註解 + memory 註明「intentional mismatch with CAD」+ 列原因。

---

## 應用到 kuo 各 topic

當前已知不一致(2026-06-01):

| Topic | ANSYS 用值 | CAD 實測 | 差異 | 動作 |
|---|---|---|---|---|
| `long2016_hexapole_halfcut` POLE_R | 3.175 mm | 3.047 mm | +4.2% | ⚠ 待對齊 |
| `long2016_hexapole_halfcut` POLE_CONE_LEN | 15.875 mm | 14.827 mm | +7.1% | ⚠ 待對齊 |
| (其他 topic 待查)| — | — | — | — |

新發現的不一致**必須加進此表**並在 [[ansys-cad-alignment]] memory 同步。

---

## 觸發片語

當使用者說以下任何句:

- 「對齊 CAD」/「對齊 SolidWorks」
- 「ANSYS 跟 CAD 不一致 / 差很多」
- 「為什麼 ANSYS 是 X mm / CAD 是 Y mm」
- 「改 mt_constants」/「改 POLE_R / POLE_CONE_LEN / 任何幾何 const」
- 「量一下這個 IGES / STEP」
- 「跑 model-check §1」(SolidWorks 視覺檢查)

→ 自動啟動此規則 + read 全文。

---

## 何時不適用

- 純後處理(MATLAB 從 .dat 抽 / 畫圖)— 不涉及幾何設定
- COMSOL 工作(COMSOL 有自己的 CAD import workflow)
- 純文字 / LaTeX 編輯
- references 文件閱讀
