# apdl/NTU_hexapole/geom/export/ — 幾何匯出建模腳本

NTU_hexapole 由 CAD（`CAD_model/NTU_hexapole/STEP/`）**參數化重建**的 geometry-only APDL 腳本 + 尺寸表。

## 現有檔案
| 檔 | 是什麼 |
|---|---|
| `pole_params.md` | **尺寸提取表**：磁極(`pole.STEP`) 2D 輪廓 9 段邊界環 + **短導柱(`短導柱.STEP`)** 階梯圓柱/扁尖 + 組合(native 同框)說明。 |
| `MT_Geom_Pole.txt` | **磁極單件**建模腳本：K→L/LARC→AL→VEXT 0.25mm → 1 volume。`UNIT_MM` 開關。原點=STEP 原點。 |
| `MT_Geom_PoleAssembly.txt` | **磁極+短導柱 組合**建模腳本：磁極段(profile+VEXT) + 短導柱(**實心** 2×CYLIND r2.5/r2.49、軸線 x=20.5、z0.25→17.25) → **3 volumes、暫不 boolean**。短導柱位置＝OCC 烘焙 pole組合 transform 得偏移 (+20.5,0,+11.25)（**非** raw 同框）。`UNIT_MM` 開關。 |
| `MT_Geom_UpperAssembly.txt` | **upper_assembly = 3 極 + 6 導柱 + yoke**（來源 `建棋模擬_上層/{上層組合,YOKE}.step`）。3 極 WPROTA 0/120/240°；6 導柱 radius20.5 @0/60/…/300°（base z 交替 0.25/−0.5、step z11.25、r2.5/r2.49）；yoke=外盤(sharp 59.8×53.8)−中心 cutout(31×20.5)−6 孔 r2.5、厚5.5 @z[11.25,16.75]（一次 VSBV ALL）→ **16 volumes**。座標框 T(Rz180+z−9.55)。尺寸 OCC 抽自 CAD。跑通 0 error。**IGESOUT 改出 `IGES/`（不進 model_check）**。 |
| `scripts/make_upper_assembly_step.py` | **upper_assembly 交付 STEP 腳本**（OCC）：讀 上層組合.STEP **整個 shape 套 T** → `model_check/NTU_hexapole/upper_assembly.step`（**全 10 solids exact CAD**：3 極 Z[0,0.25]/yoke Z[11.25,16.75]/6 導柱；腳本讀回自驗）。 |
| `full_assembly_params.md` | **full_assembly 尺寸表**（從 `總組合.STEP` OCC 世界座標量測，**對齊 STEP 原點、不平移**）：6 磁極兩層 bbox（下 z[8.80,9.05]@0/120/240°、上 z[9.55,9.80]@60/180/300°）+ **6 導柱**（r2.5、徑向 20.5、頂 z26.80，長導柱 z[9.05,26.80]/短導柱 z[9.80,26.80]）。總組合 19 solids 對照。 |
| `MT_Geom_FullAssembly.txt` | **full_assembly = 整體 model 階段性建構**（覆蓋擴充「此」deck）。**階段① 6 磁極**（兩層 WPROTA loop、重用 UpperAssembly PART A 輪廓、對齊 總組合.STEP 原點不平移）；**階段② 6 導柱**（兩層 WPROTA loop、`WPOFFS,20.5,0,0`+`CYLIND,r2.5,,ZBP,26.80`、尖端簡化 r2.5 實心）；**階段③ yoke**（rounded-rect 外盤 59.8×53.8 − cutout 31×20.5 − 6 孔 r2.5、**忠實 r5 圓角**、厚5.5 @z[20.80,26.30]、一次 VSBV → 1 vol）→ **13 volumes、不 boolean**。`UNIT_MM=1`(mm，5µm 圓角公尺尺度建不起來)。跑通 0 error。**不出 IGES**（deliver-step-for-check）。coil 之後階段加。 |
| `scripts/make_full_assembly_step.py` | **full_assembly 交付 STEP 腳本**（OCC）：讀 `pole.STEP` 擺 6 磁極 + 6 導柱 OCC primitive 圓柱 r2.5 + **yoke 直接嵌 `總組合.STEP` solid06（真 CAD 含 r5 圓角）** → `model_check/NTU_hexapole/full_assembly.step`（**13 solids**，bbox 與 CAD `總組合.STEP` solid0–5+07–12+06 逐一吻合；腳本讀回自驗）。 |

## ⚠ IGES 單位坑（SolidWorks 讀成英吋 ×25.4）
ANSYS `IGESOUT` 產的 IGES，即使 units flag=2(mm) + 補上 units name `2HMM`，**SolidWorks（與 OpenCascade）仍讀成英吋、模型 ×25.4 太大**（板厚 0.25→6.35mm）。flag/name 修都無效。
→ **組合件改交 STEP**（`pole_assembly.step`，OCC 從原始 `pole.STEP` + 實心導柱組出、`SI_UNIT(.MILLI.,.METRE.)` 明確 mm）。STEP 單位穩、SolidWorks 讀對。**日後給使用者檢查的組合幾何一律出 STEP，不用 ANSYS IGES。**

## 產物（⚠ model_check 只放 STEP、不放 IGES；使用者拍板 2026-07-11）
- **`model_check/NTU_hexapole/pole_assembly.step`**（磁極+實心短導柱 單極組合；供疊 `pole組合.STEP` 檢查）。
- **`model_check/NTU_hexapole/upper_assembly.step`**（**upper_assembly = 3 極+6 導柱+yoke**，全 10 solids exact CAD 套 T；供疊 `上層組合.STEP` 檢查）← 上層檢查用這個。
- **`model_check/NTU_hexapole/full_assembly.step`**（**full_assembly = 6 磁極+6 導柱+yoke**，13 solids，對齊 `總組合.STEP` 原點；供疊 `總組合.STEP` 檢查）← 整體 model 階段性建構用這個；`db/geom/full_assembly/full_assembly.db`（13 vol）。
- `IGES/NTU_hexapole/upper_assembly_mm.iges`（deck IGESOUT，**不進 model_check**；SolidWorks 讀成 ×25.4，勿用、僅 legacy）。
- `ANSYS_data/NTU_hexapole/db/geom_export/uppergeom.db`（模型 .db）。
- （legacy `pole.iges` / `pole_assembly.iges` 已移出 model_check——model_check 只留 STEP。）

## 慣例 / 坑
- **IGES 單位 flag**：mm 版 IGESOUT 出 flag 6，跑後外部 `sed 's/,1.0,6,,/,1.0,2,,/'` патch 成 **2（mm）**。**不可**抄 hung 的 `6→1`（英吋）。
- `pole.STEP` 遠端 r5 是**凸 major arc（213°）**，頂點 (25.5,0)，故 LARC 拆兩段（見 pole_params.md）。全長 ~25.15mm、板厚 0.25mm。
- 磁極 + 短導柱(實心) 已建。**組裝位置＝用 OpenCascade(OCP) 烘焙 pole組合.STEP 的 transform 取得**（各零件在組裝檔是 local frame + transform，**不可**用 raw 座標當同框；短導柱偏移 (+20.5,0,+11.25)）。
- 組合暫不 boolean（保留兩件供位置檢查）；位置確認後可 VADD/VGLUE 併件。
- **襯套(bushing) [已建 2026-07-17]**：階段④ 6 個凸緣頂帽狀空心襯套（bore r2.5 / tube 外徑 r3.5 / flange 外徑 r7.5，用**總組合.STEP 尺寸**；⚠ 單件 襯套.STEP 是舊版 r3.0/7.0）；6 個在 r=20.5 @ 0/60/…/300°、全垂直、z[15.60,20.80]（flange z[15.60,16.60]、tube z[16.60,20.80]），套上 6 導柱（導柱填內孔）。APDL 用空心 CYLIND(外,內) + VADD 建（每襯套後 NUMCMP,VOLU 壓實避 boolean 號洞）；STEP 從總組合 solid13–18 splice。full_assembly 現 **19 volumes**。
