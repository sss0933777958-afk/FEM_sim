#!/usr/bin/env python
# iges_to_step.py — 用 OCC(OCP) 把 ANSYS IGES 轉成明確 mm 單位的 STEP，供 SolidWorks 檢查。
# ---------------------------------------------------------------------------
# 依 .claude/rules/deliver-step-for-check.md：交付檢查一律出 STEP、不出 ANSYS IGES
#   (ANSYS IGESOUT 的 IGES 被 SW/OCC 讀成英吋 ×25.4)。此腳本讀「已修 unit flag=1(inch)」的
#   ANSYS IGES(hung 用 MM=1/25.4 建、flag 改 1 後 OCC 讀回正確 mm)，寫 write.step.unit=MM 的 STEP，
#   並讀回驗 bbox + 6 個最內側 tip 距 WP 原點 = 目標 R(mm)。
# 用法：  python iges_to_step.py <in.iges> <out.step> [R_expect_mm]
# ---------------------------------------------------------------------------
import sys, math
from OCP.IGESControl import IGESControl_Reader
from OCP.STEPControl import STEPControl_Reader, STEPControl_Writer, STEPControl_AsIs
from OCP.Interface import Interface_Static
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_VERTEX, TopAbs_SOLID
from OCP.TopoDS import TopoDS
from OCP.BRep import BRep_Tool
from OCP.Bnd import Bnd_Box
from OCP.BRepBndLib import BRepBndLib

def main(in_iges, out_step, r_expect=None):
    r = IGESControl_Reader(); r.ReadFile(in_iges); r.TransferRoots(); sh = r.OneShape()
    Interface_Static.SetCVal_s("write.step.unit", "MM")
    w = STEPControl_Writer(); w.Transfer(sh, STEPControl_AsIs); w.Write(out_step)
    # verify by reading STEP back
    rr = STEPControl_Reader(); rr.ReadFile(out_step); rr.TransferRoots(); s2 = rr.OneShape()
    bb = Bnd_Box(); BRepBndLib.Add_s(s2, bb); xmin, ymin, zmin, xmax, ymax, zmax = bb.Get()
    seen = set(); ds = []
    e = TopExp_Explorer(s2, TopAbs_VERTEX)
    while e.More():
        p = BRep_Tool.Pnt_s(TopoDS.Vertex_s(e.Current()))
        k = (round(p.X(), 5), round(p.Y(), 5), round(p.Z(), 5))
        if k not in seen:
            seen.add(k); ds.append(math.sqrt(p.X()**2 + p.Y()**2 + p.Z()**2))
        e.Next()
    ds.sort()
    ns = 0; es = TopExp_Explorer(s2, TopAbs_SOLID)
    while es.More(): ns += 1; es.Next()
    tip = [round(d, 4) for d in ds[:6]]
    msg = "wrote %s | solids=%d | bbox x[%.2f,%.2f] z[%.2f,%.2f] mm | 6 tip dist=%s mm" % (
        out_step, ns, xmin, xmax, zmin, zmax, tip)
    if r_expect is not None:
        ok = all(abs(d - float(r_expect)) < 1e-3 for d in ds[:6])
        msg += " | expect %.3f -> %s" % (float(r_expect), "OK" if ok else "MISMATCH")
    print(msg)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else None)
