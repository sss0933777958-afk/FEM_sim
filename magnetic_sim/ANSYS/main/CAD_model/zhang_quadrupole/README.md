# CAD_model/zhang_quadrupole/ — Zhang 4-pole Quadrupole 原始幾何（source of truth）

**用途**：Zhang 4-pole（0/90/180/270° lozenge + yoke + pillars）Quadrupole 的 SolidWorks 原檔；此 model 幾何的 source of truth。

**內容**：子夾 `SLDPRT/` 與 `STEP/`。
- `SLDPRT/` — 主線 `0.46mm_baseline.SLDPRT` 及多 revision `_V2`…`_V10`、變體 `0.46mm_20um.SLDPRT`、`0.67mm_baseline.SLDPRT`、原始件 `zhipeng.SLDPRT`。
- `STEP/` — 目前為**空夾**（尚未匯出 STEP；需要時依 `doc/workflows/cad-export.md` 產出）。

**資料來源 / 流向**：pipeline 起點。SolidWorks → STEP → `IGES/zhang_quadrupole/` → `IGES_converted/zhang_quadrupole/` → `apdl/zhang_quadrupole/geom` → ANSYS。目前下游 IGES 只有 `Zhang_Quadrupole_Lp0p405.iges`。

**命名 / 慣例**：以極尖距（`0.46mm` / `0.67mm`）+ revision（`_V*`）命名；charge fit 已得 A≈493µm（err 4.4%，target 490）。`zhipeng` 為原作者原始件名。

**相關**：`../README.md`、`../../CLAUDE.md`、`ansys-cad-alignment.md`、`doc/workflows/cad-export.md`。
