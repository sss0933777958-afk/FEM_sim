"""
Convert the sensor-refine ANSYS IGES -> mm STEP for CAD check (deliver-step-for-check.md).

Source = ANSYS IGESOUT from MT_Geom_SensorRefine.txt (steel MAT_MT + 6 sensor R0.3 spheres),
written to the claude scratchpad. ANSYS IGES is read by OCC/SolidWorks as inches (x25.4);
we measure bbox and auto-calibrate scale to mm, then write an explicit mm STEP.

Deliverable: model_check/long_fei/sensor_refine.step  (poles + 6 R0.3 spheres for placement check).

Run:  python make_sensor_refine_step.py
Requires OCP / cadquery's OpenCascade.
"""
from OCP.IGESControl import IGESControl_Reader
from OCP.STEPControl import STEPControl_Reader, STEPControl_Writer, STEPControl_AsIs
from OCP.Interface import Interface_Static
from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCP.gp import gp_Trsf, gp_Pnt
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_SOLID
from OCP.Bnd import Bnd_Box
from OCP.BRepBndLib import BRepBndLib
from OCP.GProp import GProp_GProps
from OCP.BRepGProp import BRepGProp

ROOT    = r"G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main"
SCRATCH = r"C:\Users\Kuo\AppData\Local\Temp\claude\G--my-workspace\b77e29dd-bd9c-479e-93f7-d3abb6bbd5e1\scratchpad"
IGES = SCRATCH + r"\sensor_refine.iges"
OUT  = ROOT + r"\model_check\long_fei\sensor_refine.step"

# expected 6 sensor sphere centers (mm, ANSYS frame) + R
SC = [(4.8012, 0.0, -14.3477), (-3.1322, 0.0, -8.7220), (-2.4006, 4.1580, -14.3477),
      (1.5661, -2.7125, -8.7220), (1.5661, 2.7125, -8.7220), (-2.4006, -4.1580, -14.3477)]
R_SENS = 0.3

def bbox(shape):
    b = Bnd_Box(); BRepBndLib.Add_s(shape, b); return b.Get()  # xm,ym,zm,xM,yM,zM

def solids(shape):
    out = []; e = TopExp_Explorer(shape, TopAbs_SOLID)
    while e.More(): out.append(e.Current()); e.Next()
    return out

def vol(s):
    g = GProp_GProps(); BRepGProp.VolumeProperties_s(s, g); return g.Mass()

def centroid(s):
    g = GProp_GProps(); BRepGProp.VolumeProperties_s(s, g); p = g.CentreOfMass(); return (p.X(), p.Y(), p.Z())

# --- read ANSYS IGES ---
r = IGESControl_Reader(); r.ReadFile(IGES); r.TransferRoots(); sh = r.OneShape()
xm, ym, zm, xM, yM, zM = bbox(sh)
span = max(xM - xm, yM - ym, zM - zm)
print("read: solids=%d  bbox X[%.3f,%.3f] Y[%.3f,%.3f] Z[%.3f,%.3f]  span=%.3f"
      % (len(solids(sh)), xm, xM, ym, yM, zm, zM, span))

# --- auto-calibrate to mm (steel yoke Dout=106mm -> true span ~106mm) ---
if span < 1.0:        # metre-magnitude
    sc = 1000.0
elif span > 400.0:    # x25.4 inch-inflated (~2700)
    sc = 1.0 / 25.4
else:                 # already mm
    sc = 1.0
if sc != 1.0:
    t = gp_Trsf(); t.SetScaleFactor(sc); sh = BRepBuilderAPI_Transform(sh, t, True).Shape()
print("scale applied = %.5f  (target mm; steel span -> ~106)" % sc)

# --- write mm STEP ---
Interface_Static.SetCVal_s("write.step.unit", "MM")
w = STEPControl_Writer(); w.Transfer(sh, STEPControl_AsIs); w.Write(OUT)
print("wrote", OUT)

# --- verify: read back, list solids; 6 spheres = tiny solids (vol ~0.113 mm^3) near SC ---
rd = STEPControl_Reader(); rd.ReadFile(OUT); rd.TransferRoots(); s2 = rd.OneShape()
sv = solids(s2)
Xm, Ym, Zm, XM, YM, ZM = bbox(s2)
print("STEP: solids=%d  bbox X[%.2f,%.2f] Y[%.2f,%.2f] Z[%.2f,%.2f]" % (len(sv), Xm, XM, Ym, YM, Zm, ZM))
sph_vol = 4.0/3.0*3.141592653589793*R_SENS**3
small = [s for s in sv if vol(s) < 5*sph_vol]   # spheres ~0.113 mm^3
print("expected sphere vol = %.4f mm^3 ; found %d small solids" % (sph_vol, len(small)))
for s in small:
    cx, cy, cz = centroid(s)
    # nearest expected center
    d, best = min(((((cx-a)**2+(cy-b)**2+(cz-c)**2)**0.5), (a, b, c)) for a, b, c in SC)
    print("  sphere @ (%+.3f,%+.3f,%+.3f) mm  vol=%.4f  -> nearest SC (%+.3f,%+.3f,%+.3f) dist=%.3fmm"
          % (cx, cy, cz, vol(s), best[0], best[1], best[2], d))
