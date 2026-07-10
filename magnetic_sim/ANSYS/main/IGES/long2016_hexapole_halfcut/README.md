# IGES/long2016_hexapole_halfcut/ — Long Fei 六極半切 ANSYS 匯出 IGES

**用途**：Long Fei 6 極下極半切 hexapole（CAD topic 名 `long_fei`）的 ANSYS 匯出 IGES 幾何，建 mesh 用。**必須與 `../../IGES_converted/long2016_hexapole_halfcut/` 同步**。

**內容**：`long2016_hexapole_no_gap.iges`（主幾何）、`..._Geom_HollowProt.iges`（hollow protrusion）、`..._Geom_gap200um.iges`（protrusion 200µm 氣隙變體）、`..._Geom_P1BaseGap.iges`（**6 極 base 各 200µm 磁隙＝貼合 base 實形的 conformal 實體 slab**(下極 local x[42.5,42.7]、上極 local x[41,41.2]，VSBW 雙刀切真實鋼、非挖洞，→ μ_r=1 於 sim)；**檔名 legacy 現含全 6 極**；metre/flag-6）、`..._Geom_sphtip.iges`（球尖端變體）。

**資料來源 / 流向**：來自 `CAD_model/long_fei/STEP` → 此 IGES → `IGES_converted/long2016_hexapole_halfcut/` 供 `IGESIN` → `apdl/.../geom` → ANSYS。

**命名 / 慣例**：⚠ 此處 topic 用 `long2016_hexapole_halfcut`，CAD 端用 `long_fei`，**同一物理模型**。變體後綴 `HollowProt`/`gap200um`/`sphtip`。kuo/long 體系 MKS 轉換用 **flag 2 / 重 export**，不可用 hung sed。

**相關**：`../README.md`、`../../IGES_converted/long2016_hexapole_halfcut/`、`doc/workflows/iges-sync-quick.md`、`iges-model-id.md`。
