# CAD_model/ — SolidWorks / STEP 原始幾何

幾何的 **source of truth**。所有 ANSYS 幾何尺寸必須對齊這裡（見 `.claude/rules/ansys-cad-alignment.md`）。

## 結構：`<topic>/{SLDPRT,STEP}/`

```
CAD_model/
├── kuo_quadrupole/   {SLDPRT, STEP}
├── long_fei/         {SLDPRT, STEP}   ← 注意：此處 topic 名用 long_fei
└── zhang_quadrupole/ {SLDPRT, STEP}
```

> ⚠ 命名不一致：CAD 這裡用 `long_fei`，但 `ANSYS_data/`、`IGES/`、`apdl/` 用 `long2016_hexapole_halfcut`。指同一個物理模型。

## 檔案類型
- `*.SLDPRT` / `*.SLDASM` — SolidWorks 原檔（可編輯的母檔）。
- `*.STEP` — 中性格式匯出（解析尺寸、跨軟體交換）。

## 規則
- 改 ANSYS 幾何 / `mt_constants` 前**必先量對應 CAD**，不一致須通報使用者拍板（per `ansys-cad-alignment.md`）。
- 出 STEP/IGES 流程：`doc/workflows/cad-export.md`。
- STEP 抽參數做參數表：`doc/workflows/step-to-apdl.md`。
