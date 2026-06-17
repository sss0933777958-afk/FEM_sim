# IGES/kuo_quadrupole/ — Kuo Quadrupole ANSYS 匯出 IGES

**用途**：Kuo 4-pole MEMS Quadrupole 由 ANSYS（或 SolidWorks）匯出的 IGES 幾何，建 mesh 的幾何來源。**必須與 `../../IGES_converted/kuo_quadrupole/` 同步**。

**內容**：
- `*.iges` 幾何 — `Quadrupole_Lp046_T55.iges`、`Quadrupole_Lp046_T100.iges`、`Quadrupole_Lp067*.iges`、`Quadrupole_R0500_*.iges`（含 `_with_coils`、`_STEEL_FROM_SIM`、`_F20`）、`Quadrupole_ScaleDown.iges`/`_V2`、`TEST_only_coils.iges`。
- log/暫存（可重生，**勿當交付檔**）— `build_*.out`、`*_geom*.out`/`.err`、`diagF20_*`、`kuo_steel_from_sim*`、`cleanup-ansys-*.bat`。

**資料來源 / 流向**：來自 `CAD_model/kuo_quadrupole/STEP` → ANSYS 匯出至此 → 轉換成 `IGES_converted/kuo_quadrupole/` 供 `IGESIN` → `apdl/kuo_quadrupole/geom`。

**命名 / 慣例**：variant 以幾何特徵命名（`Lp046`/`Lp067` 極尖距、`T55`/`T100` 厚度、`R0500` 等）。kuo（MKS）轉換用 **flag 2 或重 export**，**不可抄 hung 的 `sed 6→1`**（會把 mm 讀成 inch）。`.out/.err` 是 log，非幾何。

**相關**：`../README.md`、`../../IGES_converted/kuo_quadrupole/`、同步 SOP `doc/workflows/iges-sync-quick.md`、模型識別 `iges-model-id.md`。
