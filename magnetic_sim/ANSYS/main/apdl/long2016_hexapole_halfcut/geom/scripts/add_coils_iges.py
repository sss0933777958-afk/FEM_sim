#!/usr/bin/env python
"""
add_coils_iges.py -- long2016 half-cut hexapole: add 6 solid COIL RINGS around the
6 magnetic-core posts, on top of the already-made GAP IGES, and write a mm IGES.

Task (user, 2026-07-28): wrap each of the 6 cores with a coil = annular cylinder
inner R=5 mm, outer R=8 mm, height 7 mm -- i.e. the exact APDL SOURC36 coil
(COIL_IN_R=5, COIL_OUT_R=8, COIL_H=PROT_H=7), so it can be assigned as a 70-turn
coil in Maxwell. The coils are kept as separate solid bodies (no boolean with steel),
to be meshed as stranded windings.

Placement (measured on the STEP; the STEP shares the APDL z-origin: yoke at z[0,2]):
  - post centres at radius 47.5 mm (= APDL YOKE_MID_R) along each pole azimuth
  - LOWER poles az 0/120/240 : coil z in [-7, 0]  (7 mm just below the yoke bottom)
  - UPPER poles az 60/180/300: coil z in [ 2, 9]  (7 mm just above the yoke top)
  This reproduces the APDL coil rings exactly (cf. MT_Geom_mm_WithCoil.txt, but with
  the exact 5/8/7 dims, no 0.5 mm COMSOL air gap).

Base geometry = the gap IGES (user picked). Output = a NEW IGES (keeps gap-only intact).
Reuses the OCC read/write pattern of make_gap_iges.py.
"""
import sys, os, math
from OCP.STEPControl import STEPControl_Reader
from OCP.IGESControl import IGESControl_Writer, IGESControl_Reader
from OCP.IFSelect import IFSelect_RetDone
from OCP.Interface import Interface_Static
from OCP.Bnd import Bnd_Box
from OCP.BRepBndLib import BRepBndLib
from OCP.BRepPrimAPI import BRepPrimAPI_MakeCylinder
from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut
from OCP.gp import gp_Pnt, gp_Ax2, gp_Dir
from OCP.GProp import GProp_GProps
from OCP.BRepGProp import BRepGProp
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_SOLID
from OCP.TopoDS import TopoDS_Compound
from OCP.BRep import BRep_Builder

# Base geometry + output. Override via argv: python add_coils_iges.py <base> <out>
# base may be STEP or IGES (auto-detected by extension). Defaults = gap version.
BASE = r"G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\model_check\long_fei\long2016_hexapole_gap.iges"
OUT  = r"G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\model_check\long_fei\long2016_hexapole_gap_coil.iges"
if len(sys.argv) >= 3:
    BASE, OUT = sys.argv[1], sys.argv[2]

def read_shape(path):
    """read STEP or IGES by extension -> OCC shape."""
    ext = os.path.splitext(path)[1].lower()
    rd = STEPControl_Reader() if ext in (".step", ".stp") else IGESControl_Reader()
    if rd.ReadFile(path) != IFSelect_RetDone:
        print("ERROR: cannot read base geometry %s" % path); sys.exit(1)
    rd.TransferRoots(); return rd.OneShape()

RLOC   = 47.5          # post centre radius (= APDL YOKE_MID_R)
COIL_IN, COIL_OUT, COIL_H = 5.0, 8.0, 7.0
LOWER_AZ = [0.0, 120.0, 240.0]   # coil z in [-7, 0]
UPPER_AZ = [60.0, 180.0, 300.0]  # coil z in [ 2, 9]
LOW_Z0, UP_Z0 = -7.0, 2.0

def total_vol(s):
    g = GProp_GProps(); BRepGProp.VolumeProperties_s(s, g); return g.Mass()
def bbox(s):
    b = Bnd_Box(); BRepBndLib.Add_s(s, b); return b.Get()
def nsolids(s):
    e = TopExp_Explorer(s, TopAbs_SOLID); n = 0
    while e.More(): n += 1; e.Next()
    return n

def coil_ring(cx, cy, z0):
    """annular solid: inner R=5, outer R=8, z in [z0, z0+7]; inner cylinder extended
    +-0.5 in z so the Cut has no coincident end faces with the outer."""
    outer = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(cx, cy, z0),      gp_Dir(0,0,1)), COIL_OUT, COIL_H).Shape()
    inner = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(cx, cy, z0-0.5),  gp_Dir(0,0,1)), COIL_IN,  COIL_H+1.0).Shape()
    return BRepAlgoAPI_Cut(outer, inner).Shape()

def main():
    base = read_shape(BASE)
    bb0 = bbox(base); n0 = nsolids(base)
    print("=== INPUT base geometry: %s ===" % os.path.basename(BASE))
    print("  bbox X[%.2f,%.2f] Y[%.2f,%.2f] Z[%.2f,%.2f]  solids=%d" %
          (bb0[0],bb0[3],bb0[1],bb0[4],bb0[2],bb0[5], n0))

    coils = []
    exp_area = math.pi*(COIL_OUT**2 - COIL_IN**2)   # 122.52 mm^2
    exp_vol  = exp_area*COIL_H                        # 857.7 mm^3
    print("=== BUILD 6 coil rings (Rin=%.0f Rout=%.0f H=%.0f, ring area=%.2f mm^2) ===" %
          (COIL_IN, COIL_OUT, COIL_H, exp_area))
    for az in LOWER_AZ:
        cx, cy = RLOC*math.cos(math.radians(az)), RLOC*math.sin(math.radians(az))
        c = coil_ring(cx, cy, LOW_Z0); coils.append(c)
        b = bbox(c)
        print("  LOWER az=%3d c=(%7.2f,%7.2f) z[%.2f,%.2f] vol=%.2f (expect %.1f)" %
              (az, cx, cy, b[2], b[5], total_vol(c), exp_vol))
    for az in UPPER_AZ:
        cx, cy = RLOC*math.cos(math.radians(az)), RLOC*math.sin(math.radians(az))
        c = coil_ring(cx, cy, UP_Z0); coils.append(c)
        b = bbox(c)
        print("  UPPER az=%3d c=(%7.2f,%7.2f) z[%.2f,%.2f] vol=%.2f (expect %.1f)" %
              (az, cx, cy, b[2], b[5], total_vol(c), exp_vol))

    builder = BRep_Builder(); comp = TopoDS_Compound(); builder.MakeCompound(comp)
    builder.Add(comp, base)
    for c in coils: builder.Add(comp, c)
    n1 = nsolids(comp)
    print("  solids: base=%d + 6 coils = %d" % (n0, n1))

    Interface_Static.SetCVal_s("write.iges.unit", "MM")
    Interface_Static.SetIVal_s("write.iges.brep.mode", 1)
    w = IGESControl_Writer("MM", 1); w.AddShape(comp)
    ok = w.Write(OUT)
    print("=== WROTE IGES ===\n  %s  ok=%s" % (OUT, bool(ok)))

    # readback fingerprint
    rb = IGESControl_Reader()
    if rb.ReadFile(OUT) != IFSelect_RetDone:
        print("WARN: cannot read back"); return
    rb.TransferRoots(); back = rb.OneShape()
    bb2 = bbox(back); n2 = nsolids(back)
    print("=== READBACK ===")
    print("  bbox X[%.2f,%.2f] Y[%.2f,%.2f] Z[%.2f,%.2f]  solids=%d (expect %d)" %
          (bb2[0],bb2[3],bb2[1],bb2[4],bb2[2],bb2[5], n2, n1))

if __name__ == "__main__":
    main()
