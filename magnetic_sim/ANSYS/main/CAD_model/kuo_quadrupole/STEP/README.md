# CAD_model/kuo_quadrupole/STEP/ — Kuo Quadrupole STEP 中性格式

**用途**：Kuo Quadrupole 的 STEP 匯出，供解析尺寸（抽參數表）與跨軟體交換用。

**內容**：主線 baseline `0.46mm_baseline.STEP` 及 revision `0.46mm_baseline_V5.STEP`…`_V10.STEP`；變體 `0.46mm_20um.STEP`、`0.67mm.STEP`、`quadrupole.STEP`、`kuo_scale_down.STEP`、`kuo_scale_down_V2.STEP`、`kuo_scale_down_fabricate.STEP`。

**資料來源 / 流向**：由 `../SLDPRT/` 母檔匯出而來 → 供 `doc/workflows/step-to-apdl.md` 抽參數寫 `apdl/kuo_quadrupole/geom/<variant>_params.md`，並產 IGES。

**命名 / 慣例**：以極尖距（`0.46mm` / `0.67mm`）+ revision（`_V*`）命名；`_20um` 指 tip 圓角等特徵變體。STEP 抽尺寸注意 B-spline 控制點 ≠ 實際邊界（per `feedback_step_geom_extraction`）。

**相關**：`../README.md`、`ansys-cad-alignment.md`、`doc/workflows/step-to-apdl.md`。
