# CAD_model/long_fei/STEP/ — Long Fei hexapole STEP 中性格式

**用途**：Long Fei 六極半切 hexapole 的 STEP 匯出，供抽參數、COMSOL import、跨軟體交換。

**內容**：`long2016_hexapolehalfcut_geom.STEP`（主幾何）、`long2016_hexapolehalfcut_geom_hollowprot_plain.STEP`、`long2016_hexapolehalfcut_comsol.STEP`（COMSOL 用）、`geom_hp_split.STEP`、`geom_withcoil.STEP`（含 6 coil rings）、`Bottom pole.STEP`、`upper_pole_sensor_location.STEP`。

**資料來源 / 流向**：由 `../SLDPRT/` 匯出 → 抽參數（`doc/workflows/step-to-apdl.md`）/ 產 IGES / COMSOL import。

**命名 / 慣例**：檔名沿 `long2016_hexapolehalfcut_*`；`hollowprot` = hollow protrusion 變體、`withcoil` = 含線圈、`split` = 半切。topic `long_fei` = 下游 `long2016_hexapole_halfcut`（同模型）。

**相關**：`../README.md`、`ansys-cad-alignment.md`、`doc/workflows/step-to-apdl.md`、COMSOL 流程 `comsol-livelink.md`。
