#!/usr/bin/env python
"""
make_gap_iges.py -- long2016 half-cut hexapole: add a thin AIR GAP into the pole
support base by an OCC boolean cut on the clean SolidWorks STEP, then write a
millimetre IGES for delivery.

Task (user, 2026-07-28): at x = -41 mm insert a 10 um air-gap slab, growing +x, through
the 220 mm^2 pole support base, as a SEPARATE SOLID BLOCK (not a void) -- i.e. the
"VSBV KEEP" idea: carve the slab out of the steel and keep it as its own solid, to be
assigned mu_r=1 (air) in ANSYS. Verified geometry:
  - solid #63 (support base)  x -44 .. -41,  y +-11, z 9..19  (its +x face at x=-41)
  - solid #71 (pole)          x -41 .. -10.86, y +-11, z 9..19 base  (its -x face at x=-41)
  -> the 220 mm^2 cross-section (22 x 10) meets face-to-face at x=-41.
So: block = (steel INTERSECT slab-box) = the pole base 10 um (220 mm^2); rest = steel
MINUS the box; deliver Compound(rest, block) -> the 10 um slab is a distinct solid
sitting between support (x<=-41) and pole (x>=-40.99). Total volume unchanged.

OCC-on-STEP route (user picked): edits the provided STEP directly -> geometry exact,
mm IGES read correctly by SW/OCC/ANSYS (NOT the ANSYS-IGESOUT inch-inflated route).
Reuses: STEP read = make_gap_step.py ; mm IGES write = hung step_to_iges.py.
"""
import sys, math
from OCP.STEPControl import STEPControl_Reader
from OCP.IGESControl import IGESControl_Writer, IGESControl_Reader
from OCP.IFSelect import IFSelect_RetDone
from OCP.Interface import Interface_Static
from OCP.Bnd import Bnd_Box
from OCP.BRepBndLib import BRepBndLib
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder
from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut, BRepAlgoAPI_Common, BRepAlgoAPI_Fuse
from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCP.gp import gp_Pnt, gp_Trsf, gp_Ax1, gp_Ax2, gp_Dir
from OCP.GProp import GProp_GProps
from OCP.BRepGProp import BRepGProp
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_SOLID
from OCP.TopoDS import TopoDS, TopoDS_Compound
from OCP.BRep import BRep_Builder

STEP = r"G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\CAD_model\long_fei\STEP\long2016_hexapolehalfcut_geom.STEP"
OUT  = r"G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\model_check\long_fei\long2016_hexapole_gap.iges"

# --- gap cutter box (mm) ---------------------------------------------------
# All 3 UPPER poles (azimuth 60/180/300) have an identical 220 mm^2 support base at
# radius 41; build the box at the 180-deg pole (axis-aligned) and rotate it about z
# by (az-180) for the other two. (Lower poles 0/120/240 have no base there.)
POLES_AZ = [180.0, 60.0, 300.0]
GAP_X0 = -41.00          # gap plane (support-base / pole interface) at the 180-deg pole
GAP_T  = 0.010           # 10 um, grows +x (radially inward, toward centre)
# CUT box: slightly oversized in y,z (air margin) -> clean cut of the 220 mm^2 steel,
#   WITHOUT touching the base plate (z 0-2) or solid #15 (z 2-9).
BOX_Y  = (-12.0, 12.0)
BOX_Z  = (  8.0, 20.0)
# BLOCK: exact support-base cross-section (y+-11, z9-19 = 220 mm^2) -> the kept air-gap solid.
BLK_Y  = (-11.0, 11.0)
BLK_Z  = (  9.0, 19.0)

# --- extra cylindrical gap on ALL 3 upper poles' protrusion tops --------------
# z-axis disk, R=2 mm, 50 um thick, grows -z from z=19 (post top). Centre sits at
# radius 47.5 along each upper-pole azimuth (180/60/300); reuses POLES_AZ.
CYL_RLOC = 47.5          # radial distance of the cylinder centre from the z-axis
CYL_R    = 2.0
CYL_ZTOP = 19.0
CYL_T    = 0.050         # 50 um

# --- LOWER-pole T-shape support-base gap (130 mm^2), all 3 lower poles ---------
# lower poles 0/120/240; base at x=42.5 on the 0-deg pole, grows +x (radially out).
# T cross-section = 150 mm^2 (full pole base), total height 10 mm (z[-17,-7]), width 20 mm:
#   wide  y+-10 z[-17,-12] = 5mm tall = 100 mm^2 (the two bottom wings are 5mm)
#   narrow y+-5 z[-12, -7] = 5mm tall =  50 mm^2
LOWER_AZ  = [0.0, 120.0, 240.0]
LOW_X0    = 42.50
LOW_T     = 0.010        # 10 um, +x
LOW_WIDE  = (-10.0, 10.0, -17.0, -12.0)   # (ylo,yhi,zlo,zhi)  5mm tall
LOW_NARROW= ( -5.0,  5.0, -12.0,  -7.0)   # 5mm tall

# --- cylinder gap on the YOKE top, at all 3 LOWER-pole positions -------------
# z-axis disk R=2 mm, 50 um thick, at z=2 (yoke top) growing -z (z[1.95,2]);
# centre at radius 47.5 along each lower-pole azimuth (0/120/240); reuses LOWER_AZ.
YK_RLOC = 47.5
YK_R    = 2.0
YK_ZTOP = 2.0
YK_T    = 0.050         # 50 um, grows -z from z=2 into the yoke

def bbox(shape):
    b = Bnd_Box(); BRepBndLib.Add_s(shape, b)
    return b.Get()  # xmin,ymin,zmin,xmax,ymax,zmax

def total_vol(shape):
    g = GProp_GProps(); BRepGProp.VolumeProperties_s(shape, g); return g.Mass()

def nsolids(shape):
    e = TopExp_Explorer(shape, TopAbs_SOLID); n = 0
    while e.More(): n += 1; e.Next()
    return n

def steel_in_box(shape, x0, x1, y0, y1, z0, z1):
    """total steel volume inside an axis-aligned probe box (for gap fingerprinting)."""
    pb = BRepPrimAPI_MakeBox(gp_Pnt(x0, y0, z0), gp_Pnt(x1, y1, z1)).Shape()
    return total_vol(BRepAlgoAPI_Common(shape, pb).Shape())

def steel_in_box_shape(shape, box):
    """total steel volume inside a prebuilt (possibly rotated) probe box."""
    return total_vol(BRepAlgoAPI_Common(shape, box).Shape())

def main():
    r = STEPControl_Reader()
    if r.ReadFile(STEP) != IFSelect_RetDone:
        print("ERROR: cannot read STEP"); sys.exit(1)
    r.TransferRoots(); shape = r.OneShape()

    bb0 = bbox(shape); v0 = total_vol(shape); n0 = nsolids(shape)
    print("=== INPUT STEP ===")
    print("  bbox X[%.3f,%.3f] Y[%.3f,%.3f] Z[%.3f,%.3f]" % (bb0[0],bb0[3],bb0[1],bb0[4],bb0[2],bb0[5]))
    print("  solids=%d  total_vol=%.3f mm^3" % (n0, v0))
    slab0 = steel_in_box(shape, GAP_X0, GAP_X0+GAP_T, BOX_Y[0],BOX_Y[1], BOX_Z[0],BOX_Z[1])
    print("  steel in gap slab (before) = %.4f mm^3  (=> %.2f mm^2 x %g mm)" % (slab0, slab0/GAP_T, GAP_T))

    # --- for each upper pole: carve the 10 um slab (oversized cut box) and KEEP the
    #     exact 220mm^2 box as a solid block. Box built at 180-deg pole, rotated about z. ---
    zax = gp_Ax1(gp_Pnt(0,0,0), gp_Dir(0,0,1))
    base_cut = BRepPrimAPI_MakeBox(gp_Pnt(GAP_X0, BOX_Y[0], BOX_Z[0]),
                                   gp_Pnt(GAP_X0+GAP_T, BOX_Y[1], BOX_Z[1])).Shape()
    base_blk = BRepPrimAPI_MakeBox(gp_Pnt(GAP_X0, BLK_Y[0], BLK_Z[0]),
                                   gp_Pnt(GAP_X0+GAP_T, BLK_Y[1], BLK_Z[1])).Shape()
    rest = shape; blocks = []
    print("=== CARVE + KEEP BLOCK on %d upper poles ===" % len(POLES_AZ))
    for az in POLES_AZ:
        t = gp_Trsf(); t.SetRotation(zax, math.radians(az-180.0))
        cutbox = BRepBuilderAPI_Transform(base_cut, t, True).Shape()
        block  = BRepBuilderAPI_Transform(base_blk, t, True).Shape()
        carved = steel_in_box_shape(rest, cutbox)          # steel about to be carved (should be ~2.2)
        rest = BRepAlgoAPI_Cut(rest, cutbox).Shape()        # remove that pole's base slab
        blocks.append(block)
        bbk = bbox(block); vblk = total_vol(block)
        print("  pole az=%3d: carved=%.4f mm^3 ; block vol=%.4f (expect 2.2) bbox X[%.2f,%.2f] Y[%.2f,%.2f] Z[%.2f,%.2f]"
              % (az, carved, vblk, bbk[0],bbk[3],bbk[1],bbk[4],bbk[2],bbk[5]))
    # --- extra: cylindrical air-gap disk (kept as a solid) on ALL 3 upper poles' tops ---
    zb = CYL_ZTOP - CYL_T
    print("=== CARVE + KEEP CYL gap on %d upper-pole protrusion tops (R=%.1f, %g um, -z) ===" % (len(POLES_AZ), CYL_R, CYL_T*1e3))
    for az in POLES_AZ:
        cx = CYL_RLOC*math.cos(math.radians(az)); cy = CYL_RLOC*math.sin(math.radians(az))
        cutcyl = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(cx, cy, zb), gp_Dir(0,0,1)), CYL_R, CYL_T + 0.15).Shape()
        cylblk = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(cx, cy, zb), gp_Dir(0,0,1)), CYL_R, CYL_T).Shape()
        rest = BRepAlgoAPI_Cut(rest, cutcyl).Shape()
        blocks.append(cylblk)
        disk_probe  = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(cx,cy,zb),      gp_Dir(0,0,1)), CYL_R, CYL_T).Shape()
        below_probe = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(cx,cy,zb-CYL_T),gp_Dir(0,0,1)), CYL_R, CYL_T).Shape()
        print("  az=%3d c=(%.2f,%.2f): block vol=%.4f (expect %.4f) ; rest disk=%.6f (0=clean) ; below=%.4f (0.6283 intact)"
              % (az, cx, cy, total_vol(cylblk), math.pi*CYL_R**2*CYL_T,
                 steel_in_box_shape(rest, disk_probe), steel_in_box_shape(rest, below_probe)))

    # --- extra: LOWER-pole T-shape base gap (130 mm^2, kept solid) on all 3 lower poles ---
    bw = BRepPrimAPI_MakeBox(gp_Pnt(LOW_X0, LOW_WIDE[0], LOW_WIDE[2]),   gp_Pnt(LOW_X0+LOW_T, LOW_WIDE[1], LOW_WIDE[3])).Shape()
    bn = BRepPrimAPI_MakeBox(gp_Pnt(LOW_X0, LOW_NARROW[0], LOW_NARROW[2]),gp_Pnt(LOW_X0+LOW_T, LOW_NARROW[1], LOW_NARROW[3])).Shape()
    base_T = BRepAlgoAPI_Fuse(bw, bn).Shape()
    T_area = (LOW_WIDE[1]-LOW_WIDE[0])*(LOW_WIDE[3]-LOW_WIDE[2]) + (LOW_NARROW[1]-LOW_NARROW[0])*(LOW_NARROW[3]-LOW_NARROW[2])
    print("=== CARVE + KEEP lower-pole T-gap (%.0f mm^2) on %d lower poles ===" % (T_area, len(LOWER_AZ)))
    for az in LOWER_AZ:
        t = gp_Trsf(); t.SetRotation(zax, math.radians(az))
        Tblk = BRepBuilderAPI_Transform(base_T, t, True).Shape()
        carved = steel_in_box_shape(rest, Tblk)
        rest = BRepAlgoAPI_Cut(rest, Tblk).Shape()
        blocks.append(Tblk)
        left = steel_in_box_shape(rest, Tblk)   # steel still inside the T after cut (expect 0)
        print("  az=%3d: block vol=%.4f (expect %.4f) ; carved=%.4f ; rest-in-T=%.6f (0=clean)"
              % (az, total_vol(Tblk), T_area*LOW_T, carved, left))

    # --- extra: YOKE-top cylinder gap (kept solid) on all 3 lower-pole positions, z=2, -z 50 um ---
    yzb = YK_ZTOP - YK_T
    print("=== YOKE cyl gap @ z[%.3f,%.3f] R=%.1f on %d lower-pole positions ===" % (yzb, YK_ZTOP, YK_R, len(LOWER_AZ)))
    for az in LOWER_AZ:
        cx = YK_RLOC*math.cos(math.radians(az)); cy = YK_RLOC*math.sin(math.radians(az))
        yk_cut = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(cx,cy,yzb), gp_Dir(0,0,1)), YK_R, YK_T + 0.15).Shape()  # extend up past z=2 top
        yk_blk = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(cx,cy,yzb), gp_Dir(0,0,1)), YK_R, YK_T).Shape()
        rest = BRepAlgoAPI_Cut(rest, yk_cut).Shape()
        blocks.append(yk_blk)
        yk_disk = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(cx,cy,yzb),      gp_Dir(0,0,1)), YK_R, YK_T).Shape()
        yk_below= BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(cx,cy,yzb-YK_T), gp_Dir(0,0,1)), YK_R, YK_T).Shape()
        print("  az=%3d c=(%.2f,%.2f): block vol=%.4f (expect %.4f) ; rest disk=%.6f (0=clean) ; below=%.4f (intact)"
              % (az, cx, cy, total_vol(yk_blk), math.pi*YK_R**2*YK_T, steel_in_box_shape(rest, yk_disk), steel_in_box_shape(rest, yk_below)))

    builder = BRep_Builder(); comp = TopoDS_Compound(); builder.MakeCompound(comp)
    builder.Add(comp, rest)
    for b in blocks: builder.Add(comp, b)
    res = comp

    vrest = total_vol(rest); vblocks = sum(total_vol(b) for b in blocks)
    n1 = nsolids(res)
    print("  tile check: vol(rest)=%.4f + blocks=%.4f = %.4f  vs whole %.4f  (dV=%.5f expect ~0)"
          % (vrest, vblocks, vrest+vblocks, v0, vrest+vblocks-v0))
    print("  solids: rest=%d + %d blocks = %d (was %d, +%d)" % (nsolids(rest), len(blocks), n1, n0, n1-n0))

    # --- write mm IGES (BRep solids), same pattern as hung step_to_iges.py ---
    Interface_Static.SetCVal_s("write.iges.unit", "MM")
    Interface_Static.SetIVal_s("write.iges.brep.mode", 1)
    w = IGESControl_Writer("MM", 1); w.AddShape(res)
    ok = w.Write(OUT)
    print("=== WROTE IGES ===")
    print("  %s  ok=%s" % (OUT, bool(ok)))

    # --- read the IGES back and re-verify (delivery fingerprint) ---
    rr = IGESControl_Reader()
    if rr.ReadFile(OUT) != IFSelect_RetDone:
        print("WARN: cannot read back IGES"); return
    rr.TransferRoots(); rb = rr.OneShape()
    bb2 = bbox(rb); n2 = nsolids(rb)
    # probe the block INTERIOR (shrunk, away from coincident faces) to confirm the block survived
    blk = steel_in_box(rb, GAP_X0+0.0005, GAP_X0+GAP_T-0.0005, BLK_Y[0]+0.1,BLK_Y[1]-0.1, BLK_Z[0]+0.1,BLK_Z[1]-0.1)
    print("=== IGES READBACK ===")
    print("  bbox X[%.3f,%.3f] Y[%.3f,%.3f] Z[%.3f,%.3f]  (mm; expect ~+-74.79 / z -42.39..39.74)"
          % (bb2[0],bb2[3],bb2[1],bb2[4],bb2[2],bb2[5]))
    print("  solids=%d (expect %d incl. gap block) ; block interior has solid=%s (expect True -> block survived round-trip)"
          % (n2, n1, blk > 1e-4))

if __name__ == "__main__":
    main()
