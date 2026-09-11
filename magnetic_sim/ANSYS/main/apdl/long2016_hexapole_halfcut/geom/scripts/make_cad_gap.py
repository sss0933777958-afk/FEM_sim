"""
Build the 100 um air-gap variant of the Long Fei half-cut hexapole from the CAD STEP.

Source : CAD_model/long_fei/STEP/long2016_hexapolehalfcut_geom.STEP  (SolidWorks, 73 solids, mm)
Output : CAD_model/long_fei/STEP/long2016_hexapolehalfcut_gap.STEP   (82 solids, mm)

METHOD = CUT, not insert.  [MODIFIED 2026-08-29 per user] Nothing is translated: the gap slab is
subtracted from the steel body on the + side of each interface, and the same slab is added back
as its own solid so Maxwell/ANSYS has a body to assign vacuum to.  Consequence: the assembly
keeps its original dimensions -- pole tips, yoke height and every outer face stay exactly where
the source CAD put them; only the steel on the + side of each interface loses 0.1 mm.
(The earlier "insert" version pushed the + side out by 0.1 mm, which lifted the yoke and with it
the upper pole tips.  That side effect is gone.)

Three gap families.  B and C are cut on the + side (branch-local +x / +z); A grows the other
way -- [MODIFIED 2026-08-29 per user] the lower-pole gap runs from the interface toward -x, so
it is cut out of the INNER block (the half that carries the pole), leaving the outer block and
the post untouched.  The inner pieces are plain boxes in x, so the T cross-section at x' = 42.4
is identical to the one at 42.5.

  A) lower block radial gap   x' = 42.4 .. 42.5, cut from the inner block pieces
     contact cross-section (measured, 130 mm^2) = T shape:
        y' in [-10,10] x z in [-17,-12]   (100 mm^2)
        y' in [ -5, 5] x z in [-12, -9]   ( 30 mm^2)
     branches at azimuth 0/120/240 deg (lower poles)

  B) lower post axial gap     z = -0.1 .. 0, cut from the LOWER POST (its top end)
     contact cross-section = circle R=5 at branch-local (47.5,0) -> 78.540 mm^2
     branches at azimuth 0/120/240 deg -- the post ends up 6.9 mm long, the yoke plate is
     left intact  [MODIFIED 2026-08-29 per user: cut the posts, not the plate]

  C) upper block radial gap   x' = 41.0 .. 41.1, cut from the outer block
     contact cross-section (measured, 220 mm^2) = rectangle y' in [-11,11] x z in [9,19]
     branches at azimuth 60/180/300 deg (upper poles)

Only the three LOWER posts get an axial gap [MODIFIED 2026-08-29 per user]: the upper posts
(z = 2 .. 9, azimuth 60/180/300) are left untouched.
The posts are exact Dia-10 x 7 cylinders (V = pi*5^2*7 = 549.7787, side = pi*10*7); the
"78.8209" end-face area the loose integrator reports is numerical, not a flare.

Run:  python make_cad_gap.py [gap_mm]      (default 0.1)
Requires OCP (cadquery's OpenCascade).
"""
import sys
import math

from OCP.STEPControl import STEPControl_Reader, STEPControl_Writer, STEPControl_AsIs
from OCP.Interface import Interface_Static
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_SOLID
from OCP.TopoDS import TopoDS_Compound
from OCP.BRep import BRep_Builder
from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder
from OCP.BRepAlgoAPI import BRepAlgoAPI_Fuse, BRepAlgoAPI_Cut, BRepAlgoAPI_Common
from OCP.gp import gp_Trsf, gp_Pnt, gp_Dir, gp_Ax1, gp_Ax2
from OCP.Bnd import Bnd_Box
from OCP.BRepBndLib import BRepBndLib
from OCP.GProp import GProp_GProps
from OCP.BRepGProp import BRepGProp

ROOT = r"G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main"
SRC = ROOT + r"\CAD_model\long_fei\STEP\long2016_hexapolehalfcut_geom.STEP"
OUT = ROOT + r"\CAD_model\long_fei\STEP\long2016_hexapolehalfcut_gap.STEP"

GAP = float(sys.argv[1]) if len(sys.argv) > 1 else 0.1   # [mm]
BRANCH_AZ = (0.0, 120.0, 240.0)                          # lower-pole branches
BRANCH_AZ_UP = (60.0, 180.0, 300.0)                      # upper-pole branches

# --- interface geometry on the 0 deg branch (mm, measured from the source STEP) ------------
X_IF = 42.5                          # lower block interface plane (gap A grows toward -x)
T_LO = (-10.0, 10.0, -17.0, -12.0)   # y0,y1,z0,z1  wide part of the T
T_HI = (-5.0, 5.0, -12.0, -9.0)      # y0,y1,z0,z1  narrow part of the T
POST_C = (47.5, 0.0)                 # post axis (same for the lower and the upper post)
POST_R = 5.0                         # post radius -> end face pi*R^2 = 78.540 mm^2
POST_Z_UP = 2.0                      # upper post / yoke plate interface plane
POSTS_LO = [0, 7, 8]                 # lower post solid indices (azimuth 0/120/240); the upper
                                     # posts [15,16,17] are deliberately left without a gap
X_IF_UP = 41.0                       # upper block interface plane
R_UP = (-11.0, 11.0, 9.0, 19.0)      # y0,y1,z0,z1  upper contact rectangle (220 mm^2)


def solids(shape):
    out = []
    e = TopExp_Explorer(shape, TopAbs_SOLID)
    while e.More():
        out.append(e.Current())
        e.Next()
    return out


def vol(s, eps=1e-9):
    """Volume with a tight integration tolerance.

    The default (loose) BRepGProp accuracy is not enough here: on the B-spline yoke plate it
    misreports the cut volume by ~0.4 mm^3, i.e. ~2% of a gap body, which makes the volume
    bookkeeping below look wrong even though the geometry is exact.
    """
    g = GProp_GProps()
    BRepGProp.VolumeProperties_s(s, g, eps, True)
    return g.Mass()


def bbox(s):
    b = Bnd_Box()
    BRepBndLib.AddOptimal_s(s, b, True, False)
    return b.Get()


def rotated(shape, deg):
    if deg == 0.0:
        return shape
    t = gp_Trsf()
    t.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), math.radians(deg))
    return BRepBuilderAPI_Transform(shape, t, True).Shape()


def hits(a, b, tol=1e-9):
    """Cheap bbox overlap test, used to skip boolean ops that cannot do anything."""
    xa, ya, za, XA, YA, ZA = bbox(a)
    xb, yb, zb, XB, YB, ZB = bbox(b)
    return not (XA < xb - tol or xa > XB + tol or YA < yb - tol or
                ya > YB + tol or ZA < zb - tol or za > ZB + tol)


# --- read source ---------------------------------------------------------------------------
r = STEPControl_Reader()
r.ReadFile(SRC)
r.TransferRoots()
src = r.OneShape()
sol = solids(src)
xm, ym, zm, xM, yM, zM = bbox(src)
print("source: %s" % SRC)
print("  solids = %d   bbox = X[%.3f,%.3f] Y[%.3f,%.3f] Z[%.3f,%.3f]"
      % (len(sol), xm, xM, ym, yM, zm, zM))
assert len(sol) == 73, "unexpected solid count -- the constants above are tied to this CAD file"

# --- build the gap bodies (one solid each), at their as-designed interfaces -----------------
y0, y1, z0, z1 = T_LO
b1 = BRepPrimAPI_MakeBox(gp_Pnt(X_IF - GAP, y0, z0), GAP, y1 - y0, z1 - z0).Shape()
y0, y1, z0, z1 = T_HI
b2 = BRepPrimAPI_MakeBox(gp_Pnt(X_IF - GAP, y0, z0), GAP, y1 - y0, z1 - z0).Shape()
f = BRepAlgoAPI_Fuse(b1, b2)
f.Build()
gapA0 = f.Shape()                                        # T-shaped slab, one body

gapB0 = BRepPrimAPI_MakeCylinder(
    gp_Ax2(gp_Pnt(POST_C[0], POST_C[1], -GAP), gp_Dir(0, 0, 1)), POST_R, GAP).Shape()

y0, y1, z0, z1 = R_UP
gapC0 = BRepPrimAPI_MakeBox(gp_Pnt(X_IF_UP, y0, z0), GAP, y1 - y0, z1 - z0).Shape()

# The post gap needs a separate, taller cutting tool.  A gap body that is exactly GAP thick ends
# flush with the post's top face, and on that degenerate face-on-face coincidence the OCC booleans
# return an empty result (a fuzzy value does not help -- measured).  Overshooting the tool past the
# joint removes the coincidence; the overshoot would then eat into the yoke plate, so the tool is
# restricted to the posts.  The gap body itself stays exactly GAP thick.
gapB0_tool = BRepPrimAPI_MakeCylinder(
    gp_Ax2(gp_Pnt(POST_C[0], POST_C[1], -GAP), gp_Dir(0, 0, 1)), POST_R, 4 * GAP).Shape()

# plan entries: (gap body to keep, tool to cut with, target solid indices or None = every solid)
plan = []
for az in BRANCH_AZ:
    plan.append((rotated(gapA0, az), None, None))
    plan.append((rotated(gapB0, az), rotated(gapB0_tool, az), POSTS_LO))
for az in BRANCH_AZ_UP:
    plan.append((rotated(gapC0, az), None, None))
gaps = [b for b, _, _ in plan]

# --- cut every gap slab out of the steel ----------------------------------------------------
steel = list(sol)
ncut = 0
removed = 0.0
for body, tool, targets in plan:
    g = body if tool is None else tool
    for i, s in enumerate(steel):
        if targets is not None and i not in targets:
            continue
        if not hits(s, g):
            continue
        # Only cut where there is real material to remove.  Testing the common volume first
        # (instead of cutting and looking at the volume change) keeps the boolean away from
        # bodies that merely touch the slab face -- those came back with sub-micron slivers
        # shaved off them.
        cm = BRepAlgoAPI_Common(s, g)
        cm.Build()
        if vol(cm.Shape()) < 1e-6:
            continue
        c = BRepAlgoAPI_Cut(s, g)
        c.Build()
        res = c.Shape()
        assert len(solids(res)) == 1, "cut split solid %d into pieces" % i
        removed += vol(s) - vol(res)
        steel[i] = res
        ncut += 1
print("  cuts: %d solids, %.5f mm^3 of steel removed" % (ncut, removed))

# --- write mm STEP --------------------------------------------------------------------------
comp = TopoDS_Compound()
bb = BRep_Builder()
bb.MakeCompound(comp)
for s in steel + gaps:
    bb.Add(comp, s)

Interface_Static.SetCVal_s("write.step.unit", "MM")
w = STEPControl_Writer()
w.Transfer(comp, STEPControl_AsIs)
w.Write(OUT)
print("wrote: %s" % OUT)

# --- verify: read back ----------------------------------------------------------------------
rd = STEPControl_Reader()
rd.ReadFile(OUT)
rd.TransferRoots()
chk = rd.OneShape()
cs = solids(chk)
xm, ym, zm, xM, yM, zM = bbox(chk)
print()
print("verify (read back):")
print("  solids = %d  (73 steel + 9 gap)" % len(cs))
print("  bbox   = X[%.3f,%.3f] Y[%.3f,%.3f] Z[%.3f,%.3f]   (must equal the source bbox)"
      % (xm, xM, ym, yM, zm, zM))
vA = (100.0 + 30.0) * GAP                                  # lower T contact 130 mm^2
vB = math.pi * POST_R ** 2 * GAP                           # post disc  78.540 mm^2
vC = (R_UP[1] - R_UP[0]) * (R_UP[3] - R_UP[2]) * GAP       # upper rect 220 mm^2
print("  gap bodies expected: 3 x %.4f (lower T) + 3 x %.4f (upper rect) + 3 x %.4f (post disc) mm^3"
      % (vA, vC, vB))
for v, i in sorted((vol(s), i) for i, s in enumerate(cs) if vol(s) < 30.0):
    xm, ym, zm, xM, yM, zM = bbox(cs[i])
    print("    V=%9.5f  X[%9.4f,%9.4f] Y[%9.4f,%9.4f] Z[%8.4f,%8.4f]" % (v, xm, xM, ym, yM, zm, zM))
tot_src = sum(vol(s) for s in sol)
tot_out = sum(vol(s) for s in cs)
gap_tot = 3 * (vA + vB + vC)
print("  total volume  source %.4f -> output %.4f  (delta %+.4f, must be ~0: steel -> gap bodies)"
      % (tot_src, tot_out, tot_out - tot_src))
print("  steel volume  source %.4f -> output %.4f  (removed %.4f, expect %.4f)"
      % (tot_src, tot_out - gap_tot, tot_src - (tot_out - gap_tot), gap_tot))
