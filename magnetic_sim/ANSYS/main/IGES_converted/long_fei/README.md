# IGES_converted/long_fei/ — Long Fei 六極半切 單位轉換後 IGES

**用途**：`../../IGES/long2016_hexapole_halfcut/` 經單位/旗標轉換後、可被 ANSYS `IGESIN` 正確讀入的 Long Fei 6 極下極半切 hexapole 幾何。

**內容**：`Long2016_HexapoleHalfcut_Geom.iges`（主）、`..._Geom_HollowProt.iges`、`..._Geom_HollowProt_Plain.iges`、`..._Geom_gap200um.iges`、`..._Geom_sphtip.iges`，另含 `Geom_WithCoil.iges`（6 coil rings）、`Geom_hp_split.iges`。

**資料來源 / 流向**：由 `../../IGES/long2016_hexapole_halfcut/`（+ `CAD_model/long_fei/STEP` 直轉的含 coil / split 版）轉換而來 → `apdl/.../geom` 用 `IGESIN` → ANSYS。

**命名 / 慣例**：⚠ 同一物理模型在 `IGES_converted/` 用 topic 名 **`long_fei`**，但在 `IGES/`、`ANSYS_data/`、`apdl/` 仍用 **`long2016_hexapole_halfcut`**。與 `IGES/` 對應、須同步。MKS 轉換用 **flag 2 / 重 export**，不可抄 hung sed 6→1。`WithCoil`/`hp_split` 為由 STEP 直轉的額外變體（IGES/ 側無同名件）。

**相關**：`../README.md`、`../../IGES/long2016_hexapole_halfcut/`、`doc/workflows/iges-sync-quick.md`、`iges-model-id.md`。
