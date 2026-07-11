# apdl/NTU_hexapole/geom/export/ — 幾何匯出建模腳本

NTU_hexapole 由 CAD（`CAD_model/NTU_hexapole/STEP/`）**參數化重建**的 geometry-only APDL 腳本 + 尺寸表。

## 現有檔案
| 檔 | 是什麼 |
|---|---|
| `pole_params.md` | **尺寸提取表**：磁極(`pole.STEP`) 2D 輪廓 9 段邊界環 + **短導柱(`短導柱.STEP`)** 階梯圓柱/扁尖 + 組合(native 同框)說明。 |
| `MT_Geom_Pole.txt` | **磁極單件**建模腳本：K→L/LARC→AL→VEXT 0.25mm → 1 volume。`UNIT_MM` 開關。原點=STEP 原點。 |
| `MT_Geom_PoleAssembly.txt` | **磁極+短導柱 組合**建模腳本：磁極段(profile+VEXT) + 短導柱(**實心** 2×CYLIND r2.5/r2.49、軸線 x=20.5、z0.25→17.25) → **3 volumes、暫不 boolean**。短導柱位置＝OCC 烘焙 pole組合 transform 得偏移 (+20.5,0,+11.25)（**非** raw 同框）。`UNIT_MM` 開關。 |

## ⚠ IGES 單位坑（SolidWorks 讀成英吋 ×25.4）
ANSYS `IGESOUT` 產的 IGES，即使 units flag=2(mm) + 補上 units name `2HMM`，**SolidWorks（與 OpenCascade）仍讀成英吋、模型 ×25.4 太大**（板厚 0.25→6.35mm）。flag/name 修都無效。
→ **組合件改交 STEP**（`pole_assembly.step`，OCC 從原始 `pole.STEP` + 實心導柱組出、`SI_UNIT(.MILLI.,.METRE.)` 明確 mm）。STEP 單位穩、SolidWorks 讀對。**日後給使用者檢查的組合幾何一律出 STEP，不用 ANSYS IGES。**

## 產物
- `model_check/NTU_hexapole/pole.iges`（磁極單件；**mm flag=2**）。
- **`model_check/NTU_hexapole/pole_assembly.step`**（磁極+實心短導柱 組合；**STEP mm，SolidWorks 讀對**，供疊 `pole組合.STEP` 檢查）← 檢查用這個。
- `model_check/NTU_hexapole/pole_assembly.iges`（同幾何 IGES；**SolidWorks 讀成 ×25.4 英吋、勿用**，僅留存）。
- `ANSYS_data/NTU_hexapole/db/geom_export/{pole_geom,poleasm}.db`（模型 .db；scratch half-clean）。

## 慣例 / 坑
- **IGES 單位 flag**：mm 版 IGESOUT 出 flag 6，跑後外部 `sed 's/,1.0,6,,/,1.0,2,,/'` патch 成 **2（mm）**。**不可**抄 hung 的 `6→1`（英吋）。
- `pole.STEP` 遠端 r5 是**凸 major arc（213°）**，頂點 (25.5,0)，故 LARC 拆兩段（見 pole_params.md）。全長 ~25.15mm、板厚 0.25mm。
- 磁極 + 短導柱(實心) 已建。**組裝位置＝用 OpenCascade(OCP) 烘焙 pole組合.STEP 的 transform 取得**（各零件在組裝檔是 local frame + transform，**不可**用 raw 座標當同框；短導柱偏移 (+20.5,0,+11.25)）。
- 組合暫不 boolean（保留兩件供位置檢查）；位置確認後可 VADD/VGLUE 併件。
- **襯套(bushing，transform #449 含旋轉~7.2°) 與「複製成 6 極 hexapole」為後續，尚未建**（屆時同樣用 OCC 烘焙 transform 定位）。
