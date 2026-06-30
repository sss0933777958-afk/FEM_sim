# IGES_converted/long_fei/ — Long Fei 六極半切 單位轉換後 IGES

**用途**：`../../IGES/long2016_hexapole_halfcut/` 經單位/旗標轉換後、可被 ANSYS `IGESIN` 正確讀入的 Long Fei 6 極下極半切 hexapole 幾何。

**內容**：`Long2016_HexapoleHalfcut_Geom.iges`（主）、`..._Geom_HollowProt.iges`、`..._Geom_HollowProt_Plain.iges`、`..._Geom_gap200um.iges`、`..._Geom_sphtip.iges`、`..._Geom_WPsphere.iges`（鐵件總成 ＋ WP 7mm 空氣球殼，raw 重疊、CAD 檢視用；mm/flag-2），另含 `Geom_WithCoil.iges`（6 coil rings）、`Geom_hp_split.iges`、`SinglePoleFilled.iges`（**單一下極：削平半錐填回完整圓錐** ＋ 4 塊支撐座 ＋ 1 根 protrusion 鐵柱；無 yoke / 無上極 / 無 coil 實體；mm/flag-2）。

**`..._Geom_WPsphere.iges` 出處**：由 `apdl/long2016_hexapole_halfcut/geom/export/MT_Geom_Export_mm_WPsphere.txt`（＝ `MT_Geom_Export_mm.txt` 鐵件 mm 建構 ＋ APDL `SPHERE` 加 WP 球，不做布林）`IGESOUT` → flag 6→2 patch。球心 z = −12.71 mm、R = 7 mm，與 FEM `V7` 同。

**`SinglePoleFilled.iges` 出處**：由 `apdl/long2016_hexapole_halfcut/geom/export/MT_Geom_Export_mm_SinglePole.txt`（＝ `MT_Geom_Export_mm.txt` 改：只建 1 下極 ＋ **略過 half-cut BLOCK+VSBV → 完整圓錐**、只建該極 1 根 protrusion、無 yoke、無上極、無 coil）`IGESOUT` → flag 6→2 patch。只出 mm 版（使用者要求）。

**資料來源 / 流向**：由 `../../IGES/long2016_hexapole_halfcut/`（+ `CAD_model/long_fei/STEP` 直轉的含 coil / split 版）轉換而來 → `apdl/.../geom` 用 `IGESIN` → ANSYS。

**命名 / 慣例**：⚠ 同一物理模型在 `IGES_converted/` 用 topic 名 **`long_fei`**，但在 `IGES/`、`ANSYS_data/`、`apdl/` 仍用 **`long2016_hexapole_halfcut`**。與 `IGES/` 對應、須同步。MKS 轉換用 **flag 2 / 重 export**，不可抄 hung sed 6→1。`WithCoil`/`hp_split` 為由 STEP 直轉的額外變體（IGES/ 側無同名件）。

**相關**：`../README.md`、`../../IGES/long2016_hexapole_halfcut/`、`doc/workflows/iges-sync-quick.md`、`iges-model-id.md`。
