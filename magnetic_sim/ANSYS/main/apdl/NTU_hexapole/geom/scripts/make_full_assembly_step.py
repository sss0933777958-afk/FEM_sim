# -*- coding: utf-8 -*-
"""
NTU_hexapole full_assembly — model_check STEP（階段性建構）。
階段①：6 磁極（讀 CAD 單極 solid pole.STEP，OCC 擺 6 份到 總組合.STEP 量到的世界座標）：
  下層 z=8.80 @ 0/120/240°；上層 z=9.55 @ 60/180/300°（局部 +x=凸端朝該方位，尖端朝中心 WP=原點）。
階段②[ADDED]：6 導柱（guide posts）——OCC primitive 實心圓柱 r2.5（簡化 CAD r2.49 尖段），
  沿各磁極方位軸 r=20.5，垂直 +z，底坐磁極板頂面、頂齊平 z=26.80：
  長導柱 ×3 @0/120/240° z[9.05,26.80]（長17.75）；短導柱 ×3 @60/180/300° z[9.80,26.80]（長17.00）。
階段③[ADDED]：yoke——直接嵌入 CAD 真 yoke（總組合.STEP solid06，含 r5 圓角），已在世界座標對齊原點、直接 add。
  rounded-rect 盤 59.8×53.8 − 中心 cutout 31×20.5 − 6 孔 r2.5，厚 5.5 @ z[20.80,26.30]。
產乾淨 mm STEP → model_check/NTU_hexapole/full_assembly.step（deliver-step-for-check：出 STEP 不出 ANSYS IGES）。
之後階段擴充 coil 時再往 compound 加 solid。
"""
import math
from OCP.STEPControl import STEPControl_Reader, STEPControl_Writer, STEPControl_AsIs
from OCP.Interface import Interface_Static
from OCP.TopoDS import TopoDS_Compound
from OCP.BRep import BRep_Builder
from OCP.gp import gp_Trsf, gp_Ax1, gp_Ax2, gp_Pnt, gp_Dir, gp_Vec
from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCP.BRepPrimAPI import BRepPrimAPI_MakeCylinder
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_SOLID
from OCP.Bnd import Bnd_Box
from OCP.BRepBndLib import BRepBndLib

ROOT = r"G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main"
POLE = ROOT + r"\CAD_model\NTU_hexapole\STEP\pole.STEP"
ASM  = ROOT + r"\CAD_model\NTU_hexapole\STEP\總組合.STEP"   # [ADDED 階段③] yoke = solid06
OUT  = ROOT + r"\model_check\NTU_hexapole\full_assembly.step"

# (方位角 deg, 層 z_bottom mm) — 下層 0/120/240 @ z=8.80、上層 60/180/300 @ z=9.55
PLACE = [(0,8.80),(120,8.80),(240,8.80),(60,9.55),(180,9.55),(300,9.55)]
# [ADDED 階段②] 6 導柱 (方位角 deg, 底 z mm)；徑向 20.5、半徑 2.5、頂 z=26.80
RPOST, RPOST_R, ZTOP = 20.5, 2.5, 26.80
POSTS = [(0,9.05),(120,9.05),(240,9.05),(60,9.80),(180,9.80),(300,9.80)]

rd = STEPControl_Reader(); rd.ReadFile(POLE); rd.TransferRoots(); pole = rd.OneShape()
comp = TopoDS_Compound(); bb = BRep_Builder(); bb.MakeCompound(comp)
for ang, z in PLACE:
    rot = gp_Trsf(); rot.SetRotation(gp_Ax1(gp_Pnt(0,0,0), gp_Dir(0,0,1)), math.radians(ang))
    tra = gp_Trsf(); tra.SetTranslation(gp_Vec(0,0,z))
    full = tra.Multiplied(rot)                      # 先繞 z 轉 ang，再平移 z
    shp = BRepBuilderAPI_Transform(pole, full, True).Shape()
    bb.Add(comp, shp)

# [ADDED 階段②] 6 導柱：垂直實心圓柱，中心 (20.5cosθ,20.5sinθ,zbot)、高 = 26.80-zbot
for ang, zb in POSTS:
    x = RPOST*math.cos(math.radians(ang)); y = RPOST*math.sin(math.radians(ang))
    cyl = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(x,y,zb), gp_Dir(0,0,1)), RPOST_R, ZTOP-zb).Shape()
    bb.Add(comp, cyl)

# [ADDED 階段③] yoke = 總組合.STEP solid06；[ADDED 階段④] 6 襯套 = solid13..18（真 CAD，含 r3.5/7.5 尺寸）
#   兩者已在世界座標對齊原點，直接 add。總組合 solid 序：00–05 極、06 yoke、07–12 導柱、13–18 襯套。
rda = STEPControl_Reader(); rda.ReadFile(ASM); rda.TransferRoots(); asm = rda.OneShape()
WANT = {6} | set(range(13, 19))      # yoke(solid06) + 6 襯套(solid13..18)
exa = TopExp_Explorer(asm, TopAbs_SOLID); k = 0
while exa.More():
    if k in WANT:
        bb.Add(comp, exa.Current())
    k += 1; exa.Next()

Interface_Static.SetCVal_s("write.step.unit","MM")
w = STEPControl_Writer(); w.Transfer(comp, STEPControl_AsIs); w.Write(OUT)
print("wrote", OUT)

# 讀回驗每 solid bbox（mm）
rd2 = STEPControl_Reader(); rd2.ReadFile(OUT); rd2.TransferRoots(); sh = rd2.OneShape()
ex = TopExp_Explorer(sh, TopAbs_SOLID); i = 0
while ex.More():
    b = Bnd_Box(); BRepBndLib.Add_s(ex.Current(), b); xm,ym,zm,xM,yM,zM = b.Get()
    print("  solid%d: X[%.2f,%.2f] Y[%.2f,%.2f] Z[%.2f,%.2f]"%(i,xm,xM,ym,yM,zm,zM)); i+=1; ex.Next()
print("EXPECT 19 solids: 6 磁極 + 6 導柱(r2.5) + 1 yoke(X±29.9/Y±26.9/z[20.80,26.30]) + 6 襯套(r3.5/7.5, z[15.60,20.80])")
