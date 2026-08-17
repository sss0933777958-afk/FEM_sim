"""
Convert a Long Fei hexapole base-gap ANSYS IGES → mm STEP for CAD check.

Source = metre IGES from MT_Geom_BaseGap.txt (steel + conformal gap slabs), written
to ANSYS_data/.../db/geom_hexvariants/<name>.iges. Presets: gap_300um / gapdist / gap_400um.
Per deliver-step-for-check.md: deliver a mm STEP (not ANSYS IGES). OCC reads this
ANSYS IGES already in mm — but we still measure bbox and auto-calibrate scale to be safe.

Run:  python make_gap_step.py [<name>]   (default name = long2016_hexapole_gap_300um)
      e.g. python make_gap_step.py long2016_hexapole_gapdist
Requires OCP / cadquery's OpenCascade.
"""
import sys
from OCP.IGESControl import IGESControl_Reader
from OCP.STEPControl import STEPControl_Reader, STEPControl_Writer, STEPControl_AsIs
from OCP.Interface import Interface_Static
from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCP.gp import gp_Trsf
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_SOLID
from OCP.Bnd import Bnd_Box
from OCP.BRepBndLib import BRepBndLib
from OCP.GProp import GProp_GProps
from OCP.BRepGProp import BRepGProp

ROOT = r"G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main"
NAME = sys.argv[1] if len(sys.argv) > 1 else "long2016_hexapole_gap_300um"
IGES = ROOT + r"\model_check\long_fei\%s.iges" % NAME   # [MODIFIED 2026-08-17] db/geom/ retired; IGESOUT now writes here
OUT  = ROOT + r"\model_check\long_fei\%s.step" % NAME

def bbox(shape):
    b = Bnd_Box(); BRepBndLib.Add_s(shape, b); return b.Get()  # xm,ym,zm,xM,yM,zM

def solids(shape):
    out = []; e = TopExp_Explorer(shape, TopAbs_SOLID)
    while e.More(): out.append(e.Current()); e.Next()
    return out

# --- read ANSYS IGES ---
r = IGESControl_Reader(); r.ReadFile(IGES); r.TransferRoots(); sh = r.OneShape()
xm, ym, zm, xM, yM, zM = bbox(sh)
dz = zM - zm
print("read: solids=%d  bbox Z[%.3f,%.3f] (dz=%.3f)  X[%.3f,%.3f]" % (len(solids(sh)), zm, zM, dz, xm, xM))

# --- auto-calibrate to mm (model true z-extent = 36mm: lower base -17 .. upper base +19) ---
if dz < 1.0:        # read in metre-magnitude
    sc = 1000.0
elif dz > 300.0:    # read ×25.4 inch-inflated (would be ~914)
    sc = 1.0 / 25.4
else:               # already mm
    sc = 1.0
if sc != 1.0:
    t = gp_Trsf(); t.SetScaleFactor(sc); sh = BRepBuilderAPI_Transform(sh, t, True).Shape()
print("scale applied = %.4f  (target: mm, z-extent -> ~36)" % sc)

# --- write mm STEP ---
Interface_Static.SetCVal_s("write.step.unit", "MM")
w = STEPControl_Writer(); w.Transfer(sh, STEPControl_AsIs); w.Write(OUT)
print("wrote", OUT)

# --- verify: read back, per-solid bbox; gap slabs = solids with a ~0.3mm thin dim ---
r2 = IGESControl_Reader  # noqa (placeholder to keep import list tidy)
rd = STEPControl_Reader(); rd.ReadFile(OUT); rd.TransferRoots(); s2 = rd.OneShape()
sv = solids(s2)
Xm, Ym, Zm, XM, YM, ZM = bbox(s2)
print("STEP: solids=%d  bbox X[%.2f,%.2f] Y[%.2f,%.2f] Z[%.2f,%.2f]" % (len(sv), Xm, XM, Ym, YM, Zm, ZM))

# gap slabs = solids with volume << steel (rotation-independent; expect 6 @ ~45mm3 lower-T / ~66mm3 upper)
def vol(s):
    g = GProp_GProps(); BRepGProp.VolumeProperties_s(s, g); return g.Mass()
vols = sorted(vol(s) for s in sv)
gaps = [v for v in vols if v < 200.0]   # mm^3; steel is ~1e5
print("solid volumes [mm^3] = %s" % ", ".join("%.1f" % v for v in vols))
print("gap slabs (vol<200mm3) = %d  (gap_300um/gap_400um: 6 base; gapdist: 18 = 6 base + 12 protrusion disks)" % len(gaps))
