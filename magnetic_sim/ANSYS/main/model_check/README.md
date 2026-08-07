# model_check/ — 模型檢查交付夾（STEP 主、legacy IGES）

**用途**：ANSYS 建完的幾何**交付使用者疊 CAD/SolidWorks 檢查**的檔案放這裡。
（2026-07-11 由 `IGES_converted/` 改名而來；原「單位轉換後給 ANSYS `IGESIN`」用途已 vestigial——`IGESIN` 2025R2 已棄用，CAD 匯入改走 STEP→SpaceClaim→`ac4para`→`.db`。）

## 🔒 交付一律出 STEP（規則 `deliver-step-for-check.md`）
- **檢查交付檔＝`.step`**（`SI_UNIT(.MILLI.,.METRE.)` 明確 mm，SolidWorks/OCC 讀對）。
- **不要交付 ANSYS `IGESOUT` 的 `.iges`**：被 SolidWorks/OCC **讀成英吋 ×25.4**（板厚 0.25→6.35mm），IGES units flag 6→2 + name `2HMM` **都救不了**。
- STEP 產法：OCC（OCP）讀原始 STEP solid + primitive 實體 compound + `write.step.unit=MM`；範本 `apdl/<model>/geom/scripts/make_*_step.py`。

## 結構：`<model>/`
```
model_check/
├── NTU_hexapole/     pole_assembly.step (交付檢查) + pole*.iges (legacy)
├── long_fei/         *.iges (legacy)
└── hung_hexapole/    *.iges (legacy)
```

## 檔案類型
- `*.step` — **交付檢查用**（mm、SolidWorks 讀對）。
- `*.iges` — legacy（ANSYS IGESOUT 出的、SW 讀成 ×25.4，**勿再交付**，僅留存）。

## 相關
- 規則 `.claude/rules/deliver-step-for-check.md`（交付一律 STEP）。貼 `model_check/<model>/*` 路徑時直接從第二層認模型（`long_fei` / `hung_hexapole` / `NTU_hexapole`），不必問。
- backup/hung 有自己的 `backup/hung/IGES_converted/`（**未改名**，另一設計）。
