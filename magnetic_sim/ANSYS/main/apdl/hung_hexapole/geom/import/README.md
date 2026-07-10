# hung_hexapole — 真實 CAD 幾何匯入 APDL deck (geom/import)

把 `CAD_model/hung_hexapole/STEP/Full_Assembly.STEP`(97 body、B-spline、mm)的真實幾何匯入 ANSYS。
**依 apdl 規則此夾只放 `.txt`**；轉檔用的 `.py`/`.bat`（OCP / SpaceClaim / ac4para wrapper）在 `../scripts/`。

## 內容
- `MT_Input_ANF.txt` — **canonical 匯入 deck**：MAPDL `/INPUT` ac4para 產出的 `.anf`（APDL 幾何檔）→ 建幾何 → 產 `.db`。
  輸出：`ANSYS_data/hung_hexapole/db/from_parasolid/hung_hexapole_full.db`（97 vol/540 area/680 kp、公尺 MKS）。

## Canonical 管線（全文見 `.claude/rules/cad-import-ansys.md`）
STEP →(`../scripts/step_to_iges.py`)→ IGES(mm 檢查)；
STEP →(`../scripts/sc_step_to_parasolid.py` SpaceClaim)→ `.x_t` →(`../scripts/ac4para252.bat`/`ac4para.exe`)→ `.anf` →(**本夾 `MT_Input_ANF.txt`** `/INPUT`)→ `.db`。

## 走過但不通（此環境，已不留檔）
- MAPDL `IGESIN`：2025 R2 只剩 SMOOTH，對複雜 trimmed NURBS 面「cannot project lines to surface」→ 匯 0。
- MAPDL `~PARAIN`：呼叫 translator bat 機制不穩、ANF 產不出 → 改「ac4para 直接產 ANF + `/INPUT`」繞過。

## 單位
- IGES：mm（檢查）。ANSYS `.db`：**公尺（MKS）**（0.0625m=62.5mm，符合本專案模擬慣例）。ANF 座標內嵌 kpt/lcurv/asurf，不可文字硬 scale。
