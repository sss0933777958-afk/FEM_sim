# IGES/ — ANSYS 匯出的 IGES 幾何

由 ANSYS（或 SolidWorks）匯出的 IGES，作為建 mesh 的幾何來源。**必須與 `../IGES_converted/` 同步**。

## 結構：`<model>/`
```
IGES/
├── kuo_quadrupole/
├── long2016_hexapole_halfcut/
└── zhang_quadrupole/
```

## 檔案類型
- `*.iges` — 幾何（通常一份公尺 MKS、一份 mm；mm 供 SolidWorks 相容檢視）。
- `*.out` / `*.err` — ANSYS 匯出批次的 log（可重生，gitignore 可忽略）。
- `*.bat` — 匯出批次腳本。

## 規則
- **IGES/ 與 IGES_converted/ 必須同步**：任一 `.iges` 更新，對應 converted 版要重產（見 `doc/workflows/iges-sync-quick.md`）。
- kuo MKS 用 **flag 2 / 重 export**，**不可抄 hung 的 `sed 6→1`**。
- 看到使用者貼 `IGES/<...>.iges` 或 `IGES_converted/<...>.iges` 路徑 → 直接從路徑認出物理模型（per `.claude/rules/iges-model-id.md`），不要問「哪個模型」。
