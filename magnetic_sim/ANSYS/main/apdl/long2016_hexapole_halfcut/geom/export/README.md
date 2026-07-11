# apdl/long2016_hexapole_halfcut/geom/export/ — 幾何 only IGES 匯出腳本

**用途**：純建幾何（無 air domain、無 mesh、無 solve）→ `IGESOUT` 寫出 .iges + 3 張 PNG 視圖。供下游 mesh / CAD 對照用。

**內容**：`MT_Geom_Export*.txt`，代表檔：
- `MT_Geom_Export.txt` — 基準半切幾何（metre/MKS）。
- `MT_Geom_Export_mm.txt` — 同幾何的 mm 版（SolidWorks 相容單位）。
- `MT_Geom_Export_gap200um.txt` / `_mm_gap200um.txt` — protrusion 拆兩段、留 200 µm 氣隙的變體。
- `MT_Geom_Export_P1BaseGap.txt` / `_mm_P1BaseGap.txt`（**檔名 legacy，現含全 6 極**）— **6 極 base 各 200µm 磁隙 = 貼合 base 實形的 conformal 實體 slab（→ μ_r=1 於 sim）**。最終 VOVLAP 後跑「6 極迴圈」：每極用 `LOCAL,15,,,,ROTA`+`CSYS,15` 在該極 local frame `VSEL` 選 base（z-band 區分上/下極、排 prot/cone）＋ 旋轉 yz `WPLANE`（3 點法，local x=xc）＋ `VSBW,ALL` 雙刀切真實鋼。下極 P1/P3/P6(ROTA 0/120/240) local x[42.5,42.7]、每極 5 vol@30mm³（block1∪2∪3）；上極(ROTA 60/180/300) local x[41,41.2]、每極 1 vol@44mm³(blockB)。**鋼總量不變（VSBW 分割非移除）**；最終 91 volumes。→ metre `IGES/long2016_hexapole_halfcut/`、mm `model_check/long_fei/..._Geom_P1BaseGap.iges`。scratch `ANSYS_data/.../geom_P1basegap_{mm,metre}`。⚠ 磁隙用 VSBW 實體 slab、**不要 VSBV 挖洞**（void mesh 對不齊）；若真要 VSBV，tool(NV2)必須維持 selected。μ_r=1 於 sim deck 每極 `VSEL,LOC`+VATT 指派。
- `MT_Geom_Export_sphtip.txt` / `_HollowProt.txt` / `MT_Geom_mm_WithCoil.txt` / `MT_Geom_mm_hp_split.txt` — 球尖 / 中空柱 / 含 coil ring / 拆柱等變體。
- `MT_Geom_Export_mm_SinglePole.txt` — **單一下極：削平半錐填回完整圓錐 ＋ 4 塊支撐座 ＋ 1 根 protrusion**（略過 half-cut VSBV；無 yoke / 無上極 / 無 coil）。**只出 mm 版**（使用者要求）→ `model_check/long_fei/SinglePoleFilled.iges`。

**資料來源 / 流向**：幾何邏輯抄自 `long2016_hexapole_full/geom/MT_Geom_Export.txt` + `hexapole-long2016` 半切 logic；metre 版 IGES → `IGES/long2016_hexapole_halfcut/`，mm 版 → `model_check/long2016_hexapole_halfcut/`（units flag 6→2，兩邊必同步）；ANSYS scratch（.db/.out/.png）落在 `ANSYS_data/long2016_hexapole_halfcut/geom_export_metre*`。

**命名 / 慣例**：`MT_Geom_Export*`；`_mm` = mm 版、無後綴 = metre 版、變體 tag 接在後（`_gap200um`/`_sphtip`/`_HollowProt`）；改動標 `[ADDED]`/`[MODIFIED]`；R_norm=500e-6 等尺寸先對齊 CAD 才改。

**相關**：見 `../README.md`、`doc/workflows/{cad-export,iges-sync-quick}.md`、`.claude/rules/{apdl-editing,ansys-cad-alignment}.md`。
