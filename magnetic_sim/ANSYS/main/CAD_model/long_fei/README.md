# CAD_model/long_fei/ — Long Fei 六極半切 hexapole 原始幾何（source of truth）

**用途**：Long Fei（Long 2016 dissertation）6 極**下極半切** hexapole 的 SolidWorks 原檔與 STEP；此 model 幾何的 source of truth。

**內容**：兩個子夾 `SLDPRT/` 與 `STEP/`。
- `SLDPRT/` — `Bottom pole.SLDPRT`、`Magnetic pole Bottom_sensor_location.SLDPRT`、`upper_pole_sensor_location.SLDPRT`（上/下極 + Hall sensor 位置）。
- `STEP/` — `long2016_hexapolehalfcut_geom.STEP`（主幾何）、`..._geom_hollowprot_plain.STEP`、`..._comsol.STEP`、`geom_hp_split.STEP`、`geom_withcoil.STEP`（含 6 coil rings）、`Bottom pole.STEP`、`upper_pole_sensor_location.STEP`。

**資料來源 / 流向**：pipeline 起點。SolidWorks → STEP → `IGES/long2016_hexapole_halfcut/` → `IGES_converted/long2016_hexapole_halfcut/` → `apdl/.../geom` → ANSYS → `ANSYS_data/` → `matlab/`。

**命名 / 慣例**：⚠ **命名不一致**：CAD 這裡用 topic `long_fei`，但 `IGES/`、`IGES_converted/`、`ANSYS_data/`、`apdl/` 一律用 `long2016_hexapole_halfcut` —— **指同一個物理模型**（Long Fei 6 極下極半切 hexapole）。已知待對齊：POLE_R 3.175(ANSYS) vs 3.047(CAD) +4.2%、POLE_CONE_LEN 15.875 vs 14.827 +7.1%（per `ansys-cad-alignment.md`）。

**相關**：`../README.md`、`../../CLAUDE.md`、`ansys-cad-alignment.md`、`doc/workflows/cad-export.md`。
