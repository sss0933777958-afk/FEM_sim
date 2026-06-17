# IGES/zhang_quadrupole/ — Zhang Quadrupole ANSYS 匯出 IGES

**用途**：Zhang 4-pole（0/90/180/270°）Quadrupole 由 ANSYS 匯出的 IGES 幾何，建 mesh 用。**必須與 `../../IGES_converted/zhang_quadrupole/` 同步**。

**內容**：`Zhang_Quadrupole_Lp0p405.iges`（目前唯一一份，L_P≈0.405mm 變體）。

**資料來源 / 流向**：來自 `CAD_model/zhang_quadrupole`（STEP 目前未匯出）→ 此 IGES → `IGES_converted/zhang_quadrupole/` 供 `IGESIN` → `apdl/zhang_quadrupole/geom`。

**命名 / 慣例**：`Lp0p405` = L_P 0.405mm（`p` 代小數點）。MKS 轉換 flag 2 / 重 export，不可抄 hung sed 6→1。

**相關**：`../README.md`、`../../IGES_converted/zhang_quadrupole/`、`doc/workflows/iges-sync-quick.md`、`iges-model-id.md`。
