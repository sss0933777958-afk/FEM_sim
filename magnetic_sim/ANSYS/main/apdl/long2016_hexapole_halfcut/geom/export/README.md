# apdl/long2016_hexapole_halfcut/geom/export/ — 幾何 only IGES 匯出腳本

**用途**：純建幾何（無 air domain、無 mesh、無 solve）→ `IGESOUT` 寫出 .iges + 3 張 PNG 視圖。供下游 mesh / CAD 對照用。

**內容**：`MT_Geom_Export*.txt`，代表檔：
- `MT_Geom_Export.txt` — 基準半切幾何（metre/MKS）。
- `MT_Geom_Export_mm.txt` — 同幾何的 mm 版（SolidWorks 相容單位）。
- `MT_Geom_Export_gap200um.txt` / `_mm_gap200um.txt` — protrusion 拆兩段、留 200 µm 氣隙的變體。
- `MT_Geom_BaseGap.txt`（**現行 6 極 base 磁隙**參數化產生器；**2026-07-13 由 `MT_Geom_gap_200um.txt` 改名**）— 全鋼 `VADD` 成 1 volume 後跑「6 極迴圈」：每極 `WPROTA` 到該極 + 建「**貼合該極 base 截面的 conformal cutter**」(cutter⊆鋼) + `VSBV,vsteel,ALL,,,KEEP`（鋼挖 conformal 空腔、KEEP cutter 當 conformal 實體 slab）。gap **都往 tip 端移**（下極 gxa 42.5、上極 40.95）、截面 conformal（下極＝T 形 3 子塊 `VADD`、上極＝y±11 z[9,19] 矩形）。⚠ VADD 釋放 volume 號、下極子塊選取**用 component 快照（`CM`/`CMSEL,U`）非 `gvb0+1..+3`**。
  - **[2026-07-13] 頂部 `GAP_CFG` preset（改一個數切模型）**：`0`=gap_300um（base 300µm、無導柱 gap）/ `1`=**gapdist**（base 100µm ＋ 導柱 protrusion 兩端 gap：下極 z=−7/0 厚 78.6µm、上極 z=2/9 厚 86.87µm、圓 conformal r=5mm＝截面 78.54mm²）/ `2`=**gap_400um**（base 400µm、無導柱 gap）。`BASE_T`/`ADD_PROT`/`PROT_T_LOW`/`PROT_T_UP` 由 preset 設。self-check：base 6 slab；`ADD_PROT=1` 另 +12 導柱圓盤＝18。
  - **輸出**：metre `IGESOUT` → **`ANSYS_data/.../db/geom_hexvariants/<name>.iges`（scratch，非 main/IGES/）** → `../scripts/make_gap_step.py <name>`(OCP) 轉 **mm STEP** → `model_check/long_fei/<name>.step`。已產 `long2016_hexapole_gap_300um` / `gapdist`（18 slab）/ `gap_400um`（6 slab）三 STEP。交付檢查 per deliver-step-for-check。
- ~~`MT_Geom_Export_P1BaseGap.txt`~~（**已移除、描述過時**）— 舊 VSBW 雙刀法；**現由上面 `MT_Geom_BaseGap.txt` 的 VSBV-KEEP conformal cutter 法取代**。（磁隙用實體 slab 不 void 的理由見 memory `project_long2016_p1_base_gap`。）
- `MT_Geom_Export_sphtip.txt` / `_HollowProt.txt` / `MT_Geom_mm_WithCoil.txt` / `MT_Geom_mm_hp_split.txt` — 球尖 / 中空柱 / 含 coil ring / 拆柱等變體。
- `MT_Geom_Export_mm_SinglePole.txt` — **單一下極：削平半錐填回完整圓錐 ＋ 4 塊支撐座 ＋ 1 根 protrusion**（略過 half-cut VSBV；無 yoke / 無上極 / 無 coil）。**只出 mm 版**（使用者要求）→ `model_check/long_fei/SinglePoleFilled.iges`。

**資料來源 / 流向**：幾何邏輯抄自 `long2016_hexapole_full/geom/MT_Geom_Export.txt` + `hexapole-long2016` 半切 logic；metre 版 IGES → `IGES/long2016_hexapole_halfcut/`，mm 版 → `model_check/long2016_hexapole_halfcut/`（units flag 6→2，兩邊必同步）；ANSYS scratch（.db/.out/.png）落在 `ANSYS_data/long2016_hexapole_halfcut/geom_export_metre*`。

**命名 / 慣例**：`MT_Geom_Export*`；`_mm` = mm 版、無後綴 = metre 版、變體 tag 接在後（`_gap200um`/`_sphtip`/`_HollowProt`）；改動標 `[ADDED]`/`[MODIFIED]`；R_norm=500e-6 等尺寸先對齊 CAD 才改。

**相關**：見 `../README.md`、`doc/workflows/{cad-export,iges-sync-quick}.md`、`.claude/rules/{apdl-editing,ansys-cad-alignment}.md`。
