# 交付檢查一律出 STEP（不出 ANSYS IGES）——強制讀取

**使用者拍板（2026-07-11）**：**ANSYS 建完的幾何要交付給使用者疊 SolidWorks/CAD 檢查時，一律出 `.step`、不要出 ANSYS `IGESOUT` 的 `.iges`。**

當工作涉及：
- 把 ANSYS 建好的模型 / 組合幾何**交給使用者目視檢查**（疊 CAD、量尺寸、確認位置）
- 「出檔給你看 / 給你檢查 / 交付幾何 / 疊 CAD 對一下」

**動手前先讀完此規則。**

對應 memory：`feedback_deliver_step_for_check.md`
相關規則：`ansys-cad-alignment.md`（CAD=source of truth）、`cad-import-ansys.md`（OCC/OCP 已在用）。
相關 SOP：`doc/workflows/{cad-export,iges-sync-quick,model-check}.md`（那些是 **SolidWorks→export** 方向，跟本規則的 **ANSYS→交付** 方向不同）。

---

## 🔒 核心規則

1. **交付檢查＝STEP**。ANSYS `IGESOUT` 的 IGES **不可**直接給使用者在 SolidWorks 檢查。
2. STEP 用 **OCC（OCP python）** 產、**明確 mm 單位**（`SI_UNIT(.MILLI.,.METRE.)`）。
3. 交付前**用 OCC 讀回 STEP 驗 per-solid bbox = 真實 mm**（例：磁極板厚 0.25mm、非 6.35mm）才交付。

## ❌ 為什麼不出 ANSYS IGES（別再踩、別再重試修 IGES 單位）

ANSYS `IGESOUT` 產的 IGES，**SolidWorks 與 OpenCascade 都讀成英吋 → 模型 ×25.4 太大**（磁極板厚 0.25→**6.35mm**＝0.25×25.4）。本 session 已窮舉驗證**下列全部無效**：
- IGES global units **flag 6(m)→2(mm)** 補丁 ——無效。
- 再補 **units-name `2HMM`**（`,1.0,2,2HMM,`）——**仍無效**（SW/OCC 照樣 ×25.4）。
- 根因：ANSYS IGES 的單位宣告，SW/OCC 不吃 → fallback 英吋。**flag/name 都救不了，不要再花時間試。**

→ **直接改交 STEP**（STEP 的 `SI_UNIT` 宣告明確、SW/OCC 都讀對）。

## ✅ 怎麼產 mm STEP（OCC，範本已在 repo）

不從 ANSYS IGES 轉（OCC 讀它就 ×25.4）。改用 **原始 STEP 零件 solid ＋ OCC primitive 實體**，在正確組裝位置 compound，寫 mm STEP：

```python
from OCP.STEPControl import STEPControl_Reader, STEPControl_Writer, STEPControl_AsIs
from OCP.Interface import Interface_Static
from OCP.BRepPrimAPI import BRepPrimAPI_MakeCylinder
from OCP.BRepAlgoAPI import BRepAlgoAPI_Fuse
from OCP.gp import gp_Pnt, gp_Dir, gp_Ax2
from OCP.TopoDS import TopoDS_Compound
from OCP.BRep import BRep_Builder
# 原始零件 solid（OCC 讀 STEP 單位正確）
rd=STEPControl_Reader(); rd.ReadFile(src_step); rd.TransferRoots(); part=rd.OneShape()
# 需要的實體用 primitive 建（例：實心圓柱），擺到已知組裝位置（見 ansys-cad-alignment / OCC-baked transform）
cyl=BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(x,y,z), gp_Dir(0,0,1)), R, H).Shape()
# compound + 明確 mm 輸出
comp=TopoDS_Compound(); bb=BRep_Builder(); bb.MakeCompound(comp); bb.Add(comp,part); bb.Add(comp,cyl)
Interface_Static.SetCVal_s("write.step.unit","MM")
w=STEPControl_Writer(); w.Transfer(comp, STEPControl_AsIs); w.Write(out_step)
```

**範本腳本（可重現）**：`apdl/NTU_hexapole/geom/scripts/make_pole_assembly_step.py`
（本 session 產 `model_check/NTU_hexapole/pole_assembly.step` 用的，OCC 讀回板厚 0.25mm 已驗）。

## 落點
- 交付 STEP → `model_check/<model>/<name>.step`（沿用本次；日後要不要另開 STEP 交付夾另議）。
- ANSYS 產物 python 腳本 → `apdl/<model>/geom/scripts/`。

## 觸發片語（任一即啟動）
- 「出檔給你/使用者檢查」/「交付幾何」/「疊 CAD 檢查」/「模型建完給你看」/「輸出 IGES 給我檢查」（→ 改出 STEP）

## 何時不適用
- **ANSYS 自己 `IGESIN` 匯入**用的 `model_check/`（機器→機器，非給 SolidWorks 看）——仍走 `doc/workflows/{cad-export,iges-sync-quick}.md` 的 IGES flow。
- FEM 求解 / 後處理用 ANSYS `.db`（不涉交付檢查）。
- 使用者明確要 IGES（並自行處理單位）。
