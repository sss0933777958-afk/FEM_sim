# model_check/long_fei/ — Long Fei 六極半切 檢查用幾何（IGES + STEP）

**用途**：`../../IGES/long2016_hexapole_halfcut/` 經單位/旗標轉換後、可被 ANSYS `IGESIN` 正確讀入的 Long Fei 6 極下極半切 hexapole 幾何。**新交付一律出 STEP**（per deliver-step-for-check）。

**gap 變體 STEP（交付檢查用；皆由 `apdl/.../geom/export/MT_Geom_BaseGap.txt` 的 `GAP_CFG` preset 產、`../scripts/make_gap_step.py <name>` 轉 mm STEP，中繼 IGES 走 scratch `ANSYS_data/.../db/geom_hexvariants/`）**：
- `long2016_hexapole_gap_300um.step`（2026-07-12）：6 極 base 磁隙都往 tip 移、厚 **300µm**、截面 conformal（下極 T 形 45mm³/極、上極 y±11 矩形 66mm³/極）。7 solids（鋼 + 6 base slab）。
- **`long2016_hexapole_gapdist.step`（2026-07-13，分布式）**：把「均勻 400µm base 磁阻」拆成 **base 100µm ＋ 導柱 protrusion 兩端各一 gap**（下極 z=−7/0 厚 78.6µm、上極 z=2/9 厚 86.87µm、圓 r=5mm）。**19 solids**（鋼 32838 + 6 base slab[下15/上22mm³] + 12 導柱圓盤 conformal r=5＝**截面 78.54mm²**[下6.17/上6.82mm³]）。OCC 驗：導柱盤 @radial47.5 z=−7/0/2/9、面積 78.54mm²、base @x42.5/40.95。
- **`long2016_hexapole_gap_400um.step`（2026-07-13，均勻對照）**：只 base 磁隙、厚 **400µm**、無導柱 gap。7 solids（鋼 + 3 下極 T 60mm³ + 3 上極 88mm³）。
- gapdist 與 gap_400um 總磁阻等價（使用者串聯磁阻設計）；μ_r=1 屬 sim 階段（STEP 不帶材料）。

**內容**：`long2016_hexapole_no_gap.iges`（主）、`..._Geom_HollowProt.iges`、`..._Geom_HollowProt_Plain.iges`、`..._Geom_gap200um.iges`、`..._Geom_sphtip.iges`、`..._Geom_WPsphere.iges`（鐵件總成 ＋ WP 7mm 空氣球殼，raw 重疊、CAD 檢視用；mm/flag-2），另含 `Geom_WithCoil.iges`（6 coil rings）、`Geom_hp_split.iges`、`SinglePoleFilled.iges`（**單一下極：削平半錐填回完整圓錐** ＋ 4 塊支撐座 ＋ 1 根 protrusion 鐵柱；無 yoke / 無上極 / 無 coil 實體；mm/flag-2）。

**`..._Geom_WPsphere.iges` 出處**：由 `apdl/long2016_hexapole_halfcut/geom/export/MT_Geom_Export_mm_WPsphere.txt`（＝ `MT_Geom_Export_mm.txt` 鐵件 mm 建構 ＋ APDL `SPHERE` 加 WP 球，不做布林）`IGESOUT` → flag 6→2 patch。球心 z = −12.71 mm、R = 7 mm，與 FEM `V7` 同。

**`SinglePoleFilled.iges` 出處**：由 `apdl/long2016_hexapole_halfcut/geom/export/MT_Geom_Export_mm_SinglePole.txt`（＝ `MT_Geom_Export_mm.txt` 改：只建 1 下極 ＋ **略過 half-cut BLOCK+VSBV → 完整圓錐**、只建該極 1 根 protrusion、無 yoke、無上極、無 coil）`IGESOUT` → flag 6→2 patch。只出 mm 版（使用者要求）。

**`..._Geom_P1BaseGap.iges` 出處**（**檔名 legacy，現含全 6 極**）：由 `apdl/long2016_hexapole_halfcut/geom/export/MT_Geom_Export_mm_P1BaseGap.txt`（＝ `MT_Geom_Export_mm.txt` 於最終 VOVLAP 後加「6 極迴圈」：每極在旋轉 LOCAL frame `VSEL` 選該極 base ＋ 旋轉 yz `WPLANE` 兩面 ＋ `VSBW,ALL` 雙刀切）`IGESOUT` → flag 6→2 patch。**6 極 base 各 200µm 磁隙 = 貼合 base 實形的 conformal 實體 slab（之後 μ_r=1）**：下極 P1/P3/P6(ROTA 0/120/240) local x[42.5,42.7]、每極 5 vol@30mm³；上極(ROTA 60/180/300) local x[41,41.2]、每極 1 vol@44mm³(blockB)。**VSBW 分割真實鋼、非挖洞**（鋼總量不變）；slab 可日後每極 `VSEL,LOC` 指派 μ_r=1。⚠ protrusion 柱未切、仍跨接；只在 base 做隙。⚠ 磁隙用 VSBW 實體 slab、不要 VSBV void。圖示 `figures/long2016_hexapole_halfcut/geom_P1basegap/{fullmodel_top,fullmodel_iso,p1_base_gap_topview,p1_base_gap_iso}.png`。metre 版 → `IGES/long2016_hexapole_halfcut/..._Geom_P1BaseGap.iges`（flag 6）。μ_r=1 為 sim 階段（IGES 不帶材料）。

**資料來源 / 流向**：由 `../../IGES/long2016_hexapole_halfcut/`（+ `CAD_model/long_fei/STEP` 直轉的含 coil / split 版）轉換而來 → `apdl/.../geom` 用 `IGESIN` → ANSYS。

**命名 / 慣例**：⚠ 同一物理模型在 `model_check/` 用 topic 名 **`long_fei`**，但在 `IGES/`、`ANSYS_data/`、`apdl/` 仍用 **`long2016_hexapole_halfcut`**。與 `IGES/` 對應、須同步。MKS 轉換用 **flag 2 / 重 export**，不可抄 hung sed 6→1。`WithCoil`/`hp_split` 為由 STEP 直轉的額外變體（IGES/ 側無同名件）。

**相關**：`../README.md`、`../../IGES/long2016_hexapole_halfcut/`、`doc/workflows/iges-sync-quick.md`、`iges-model-id.md`。
