"""
Build NTU_hexapole UPPER-layer FULL ASSEMBLY (3 poles + 6 posts + yoke) as a mm STEP.

Source CAD = 建棋模擬_上層/上層組合.STEP (the whole upper assembly: 10 solids =
3 flat magnetic poles + 6 guide posts + 1 yoke). Per user: deliver the EXACT CAD.
Just apply the SAME transform used for the approved single pole to the WHOLE shape:
  T = Rz(180deg about world z) + translate z by -9.55
so pole0 stays canonical (tip +x, plate z[0,0.25]); the assembly ends up:
  3 poles Z[0,0.25] at 0/120/240deg ; yoke Z[11.25,16.75] (rounded-rect, 6 post holes,
  central cutout) ; 6 posts (radius 20.5) up to Z=17.25 (3 base z0.25, 3 base z-0.5).

Why STEP (not ANSYS IGES): ANSYS IGESOUT reads as INCH x25.4 in SolidWorks/OCC.
STEP declares SI_UNIT(.MILLI.,.METRE.) -> reads correct. Rule: deliver-step-for-check.md
model_check gets STEP only (no IGES) per user.

Run:  python make_upper_assembly_step.py   (requires OCP / cadquery's OpenCascade)
"""
import math
from OCP.STEPControl import STEPControl_Reader, STEPControl_Writer, STEPControl_AsIs
from OCP.Interface import Interface_Static
from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCP.gp import gp_Pnt, gp_Dir, gp_Ax1, gp_Trsf, gp_Vec
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_SOLID
from OCP.Bnd import Bnd_Box
from OCP.BRepBndLib import BRepBndLib

ROOT     = r"G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main"
SRC_STEP = ROOT + r"\CAD_model\NTU_hexapole\STEP\建棋模擬_上層\上層組合.STEP"
OUT_STEP = ROOT + r"\model_check\NTU_hexapole\upper_assembly.step"

# --- read whole upper assembly ---
rd = STEPControl_Reader(); rd.ReadFile(SRC_STEP); rd.TransferRoots(); sh = rd.OneShape()

# --- transform whole shape to canonical frame: Rz(180) about world z, then z -= 9.55 ---
tr = gp_Trsf(); tr.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), math.pi)
tz = gp_Trsf(); tz.SetTranslation(gp_Vec(0, 0, -9.55))
T  = gp_Trsf(); T.Multiply(tz); T.Multiply(tr)
shT = BRepBuilderAPI_Transform(sh, T, True).Shape()

# --- write mm STEP ---
Interface_Static.SetCVal_s("write.step.unit", "MM")
w = STEPControl_Writer(); w.Transfer(shT, STEPControl_AsIs); w.Write(OUT_STEP)
print("wrote", OUT_STEP)

# --- verify: read back, per-solid bbox (10 solids) ---
rd2 = STEPControl_Reader(); rd2.ReadFile(OUT_STEP); rd2.TransferRoots(); s2 = rd2.OneShape()
ex = TopExp_Explorer(s2, TopAbs_SOLID); i = 0
while ex.More():
    b = Bnd_Box(); BRepBndLib.Add_s(ex.Current(), b); xm, ym, zm, xM, yM, zM = b.Get()
    print("  solid%d: X[%.2f,%.2f] Y[%.2f,%.2f] Z[%.2f,%.2f]" % (i, xm, xM, ym, yM, zm, zM)); i += 1; ex.Next()
print("  EXPECT 10 solids: 3 poles Z[0,0.25] ; yoke Z[11.25,16.75] ; 6 posts up to Z=17.25")
