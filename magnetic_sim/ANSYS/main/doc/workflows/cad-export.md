# cad-export

從 SolidWorks SLDPRT 一次出 3 個檔:STEP(紀錄/SEMu/COMSOL)、IGES(ANSYS import)、
IGES_converted(MKS 單位修正後給 APDL 用)。

## 何時用

- 新設計第一次定型
- 改了 SLDPRT 尺寸 / 拓樸 → 要重 export

## 輸入

| 參數 | 範例 |
|---|---|
| `{topic}` | `kuo_quadrupole` / `long_fei` / `zhang_quadrupole` |
| `{sldprt}` | `magnetic_sim/ANSYS/main/CAD/{topic}/SLDPRT/<name>.SLDPRT` |
| `{basename}` | `<name>`(不含副檔名,後續 3 個檔都用同名) |

## 步驟

1. **確認目錄存在**(不存在就建):
   ```
   magnetic_sim/ANSYS/main/CAD/{topic}/{SLDPRT,STEP}/
   magnetic_sim/ANSYS/main/IGES/{topic}/
   magnetic_sim/ANSYS/main/IGES_converted/{topic}/
   ```

2. **SolidWorks 匯出 STEP** → `magnetic_sim/ANSYS/main/CAD/{topic}/STEP/{basename}.STEP`
   - File → Save As → STEP AP214
   - Options → 單位 `Meter`

3. **SolidWorks 匯出 IGES** → `magnetic_sim/ANSYS/main/IGES/{topic}/{basename}.iges`
   - File → Save As → IGES
   - Options → Output coordinate system: 預設

4. ⏸ **檢查點 → [`model-check.md` §1](model-check.md#1-solidworks-視覺檢查)**
   不過就回 SLDPRT 修,不要 re-export。

5. **產生 `IGES_converted`**:
   **kuo 是 MKS-metre** — 不可抄 hung 的 sed 公式
   (`s/,1.0,6,,/,1.0,1,,/` 是 hung 1/25.4 inches 用)。
   kuo 的做法:**重新從 SolidWorks export IGES 時把單位改成 mm**,然後把該檔
   放進 `IGES_converted/{topic}/`。或在原 `.iges` header 改 IGES unit flag 為 2(=mm)。
   詳細根因見 memory `feedback_iges_unit_conversion`。

6. ⏸ **檢查點 → [`model-check.md` §2](model-check.md#2-iges-round-trip-檢查)**

7. **紀錄版本** — 在 `magnetic_sim/ANSYS/main/CAD/{topic}/STEP/_VERSIONS.md` 加一行:
   ```
   YYYY-MM-DD | {basename} | <一句改了什麼>
   ```

## 產物

- [ ] `magnetic_sim/ANSYS/main/CAD/{topic}/STEP/{basename}.STEP`
- [ ] `magnetic_sim/ANSYS/main/IGES/{topic}/{basename}.iges`
- [ ] `magnetic_sim/ANSYS/main/IGES_converted/{topic}/{basename}.iges`
- [ ] `_VERSIONS.md` 多一行

## 常見坑

- 抄 hung sed 公式 → kuo MKS 模型尺寸炸 1000×(memory `feedback_iges_unit_conversion`)
- SLDPRT 改了只 export 一個檔 → 三檔版本漂移
- STEP header 寫 m 但 instance 寫 mm → 後續 step-to-apdl 解析會錯 1000×

## 適用 topic

`kuo_quadrupole` ✓ / `long_fei` ✓ / `zhang_quadrupole` 待 / 任何新 topic。
與 pole 配置(quadrupole / hexapole / dipole)無關。
