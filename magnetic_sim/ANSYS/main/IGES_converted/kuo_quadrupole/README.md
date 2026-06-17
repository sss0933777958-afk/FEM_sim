# IGES_converted/kuo_quadrupole/ — Kuo Quadrupole 單位轉換後 IGES（給 IGESIN）

**用途**：`../../IGES/kuo_quadrupole/` 經單位/旗標轉換後、可直接被 ANSYS `IGESIN` 正確讀入的 Kuo Quadrupole 幾何。

**內容**：`*.iges` 轉換後幾何，與 `IGES/` 一一對應 — `Quadrupole_Lp046_T55.iges`、`Quadrupole_Lp046_T100.iges`、`Quadrupole_Lp067*.iges`、`Quadrupole_R0500_*.iges`（含 `_with_coils`、`_STEEL_FROM_SIM`）、`Quadrupole_ScaleDown.iges`/`_V2`。（此處只放 `.iges`，無 log。）

**資料來源 / 流向**：由 `../../IGES/kuo_quadrupole/` 轉換而來 → `apdl/kuo_quadrupole/geom` 用 `IGESIN` 匯入建 mesh → ANSYS solve。

**命名 / 慣例**：檔名須與 `IGES/` 端**一一對應、同步**（刪/改名兩邊一起）。MKS（公尺）轉換用 **flag 2 或重 export**；**不可用 hung 的 `sed ,1.0,6,, → ,1.0,1,,`**（hung 英制專用，會把 mm 當 inch 讀成 200mm）。

**相關**：`../README.md`、`../../IGES/kuo_quadrupole/`、同步 SOP `doc/workflows/iges-sync-quick.md`、`iges-model-id.md`、`feedback_kuo_iges_export_workflow`。
