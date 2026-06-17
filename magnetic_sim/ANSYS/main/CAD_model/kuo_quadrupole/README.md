# CAD_model/kuo_quadrupole/ — Kuo 4-pole MEMS Quadrupole 原始幾何（source of truth）

**用途**：Kuo 自己的 4-pole MEMS Quadrupole 的 SolidWorks 原檔與 STEP 中性格式，是此 model 幾何的 **source of truth**；所有 ANSYS / mt_constants 幾何尺寸必須對齊這裡。

**內容**：兩個子夾 `SLDPRT/` 與 `STEP/`。
- `SLDPRT/` — 元件級 SolidWorks 原檔（`pole.SLDPRT`、`coil.SLDPRT`、`yoke.SLDPRT`、`guid_post.SLDPRT`）、組合件 `quadrupole.SLDASM`，以及縮小變體 `kuo_scale_down*.SLDPRT`（含 `_V2`、`_fabricate`）。
- `STEP/` — 對應 STEP 匯出，含主線 `0.46mm_baseline*.STEP`（多 revision V5–V10）、`0.46mm_20um.STEP`、`0.67mm.STEP`、`quadrupole.STEP`、`kuo_scale_down*.STEP`。

**資料來源 / 流向**：此處 = pipeline 起點。SolidWorks 設計 → STEP 匯出 → 下游 `IGES/kuo_quadrupole/` → `IGES_converted/kuo_quadrupole/` → `apdl/kuo_quadrupole/geom` → ANSYS solve → `ANSYS_data/` → `matlab/`。出檔 SOP 見 `doc/workflows/cad-export.md`；STEP 抽參數見 `doc/workflows/step-to-apdl.md`。

**命名 / 慣例**：variant 以幾何特徵命名（`Lp` 極尖距、`0.46mm` / `0.67mm` 變體、`ScaleDown`）。`.SLDPRT/.SLDASM` 為可編輯母檔、`.STEP` 為解析/交換用。改 ANSYS 幾何前**必先量對應 CAD**，不一致須停下通報使用者拍板（per `../../.claude/`→ `ansys-cad-alignment.md`）。

**相關**：上層 `../README.md`、`../../CLAUDE.md`；對齊規則 `ansys-cad-alignment.md`；出檔 / 抽參數 SOP 見 `doc/workflows/`。
