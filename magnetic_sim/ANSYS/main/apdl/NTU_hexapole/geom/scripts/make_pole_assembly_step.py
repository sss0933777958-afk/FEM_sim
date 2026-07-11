"""
Build NTU_hexapole pole+post ASSEMBLY as a mm STEP for user check (SolidWorks).

Why STEP (not ANSYS IGES): ANSYS IGESOUT .iges is read as INCH (x25.4) by
SolidWorks & OpenCascade even after units flag 6->2 and units-name '2HMM' patches
(both fail). STEP declares SI_UNIT(.MILLI.,.METRE.) unambiguously -> reads correct.
See rule: magnetic_sim/.claude/rules/deliver-step-for-check.md

Geometry (pole/magnet local frame = the frame pole.STEP / pole.iges use):
  - 磁極 (magnetic pole)  : original pole.STEP solid (flat planar pole, 0.25mm thick).
  - 短導柱 (guide post)    : SOLID stepped cylinder (recess omitted per user),
                             axis (x=20.5, y=0), r2.5 z0.25->11.25, r2.49 z11.25->17.25.
      post placement = 短導柱 native + assembly offset (+20.5,0,+11.25), where the
      offset came from OCC-baking pole組合.STEP transforms: T_post(#104)-T_pole(#589).

Run:  python make_pole_assembly_step.py
Requires: OCP (cadquery's OpenCascade). Read-only on inputs; writes the .step.
"""
from OCP.STEPControl import STEPControl_Reader, STEPControl_Writer, STEPControl_AsIs
from OCP.Interface import Interface_Static
from OCP.BRepPrimAPI import BRepPrimAPI_MakeCylinder
from OCP.BRepAlgoAPI import BRepAlgoAPI_Fuse
from OCP.gp import gp_Pnt, gp_Dir, gp_Ax2
from OCP.TopoDS import TopoDS_Compound
from OCP.BRep import BRep_Builder
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_SOLID
from OCP.Bnd import Bnd_Box
from OCP.BRepBndLib import BRepBndLib

ROOT = r"G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main"
POLE_STEP = ROOT + r"\CAD_model\NTU_hexapole\STEP\pole.STEP"
OUT_STEP  = ROOT + r"\model_check\NTU_hexapole\pole_assembly.step"

# --- magnetic pole = original pole.STEP solid (OCC reads STEP units correctly) ---
rd = STEPControl_Reader(); rd.ReadFile(POLE_STEP); rd.TransferRoots()
pole = rd.OneShape()

# --- SOLID guide post at pole frame: axis (20.5,0), stepped r2.5 / r2.49 ---
lo = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(20.5, 0, 0.25),  gp_Dir(0, 0, 1)), 2.50, 11.0).Shape()  # z 0.25->11.25
hi = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(20.5, 0, 11.25), gp_Dir(0, 0, 1)), 2.49, 6.0).Shape()   # z 11.25->17.25
post = BRepAlgoAPI_Fuse(lo, hi).Shape()

# --- compound (pole + post as separate solids) ---
comp = TopoDS_Compound(); bb = BRep_Builder(); bb.MakeCompound(comp)
bb.Add(comp, pole); bb.Add(comp, post)

# --- write mm STEP ---
Interface_Static.SetCVal_s("write.step.unit", "MM")
w = STEPControl_Writer(); w.Transfer(comp, STEPControl_AsIs); w.Write(OUT_STEP)
print("wrote", OUT_STEP)

# --- verify: read back, per-solid bbox must be real mm ---
rd2 = STEPControl_Reader(); rd2.ReadFile(OUT_STEP); rd2.TransferRoots(); sh = rd2.OneShape()
ex = TopExp_Explorer(sh, TopAbs_SOLID); i = 0
while ex.More():
    b = Bnd_Box(); BRepBndLib.Add_s(ex.Current(), b); xm, ym, zm, xM, yM, zM = b.Get()
    print("  solid%d: X[%.3f,%.3f] Y[%.3f,%.3f] Z[%.3f,%.3f]" % (i, xm, xM, ym, yM, zm, zM)); i += 1; ex.Next()
print("  EXPECT mm: pole Z-thickness 0.25 ; post axis x=20.5, Z[0.25,17.25]")
