#!/usr/bin/env python
"""
step_to_iges.py -- 把 hung_hexapole Full_Assembly.STEP 轉成 IGES (mm) 供 ANSYS IGESIN + 使用者檢查.
用 OCP (OpenCascade, cadquery 內含) headless 轉換.
單位坑: STEP 標 INCH -> OCC 讀取套 inch->mm (x25.4) -> ~1587mm. 偵測後 scale 回 mm (x1/25.4).
輸出座標不平移 -> 原點仍 = CAD 全域原點.
"""
import sys
from OCP.STEPControl import STEPControl_Reader
from OCP.IGESControl import IGESControl_Writer
from OCP.IFSelect import IFSelect_RetDone
from OCP.Interface import Interface_Static
from OCP.Bnd import Bnd_Box
from OCP.BRepBndLib import BRepBndLib
from OCP.gp import gp_Trsf
from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform

STEP = r"G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\CAD_model\hung_hexapole\STEP\Full_Assembly.STEP"
OUTS = [
    r"G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\IGES\hung_hexapole\hung_hexapole_full.iges",
    r"G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\model_check\hung_hexapole\hung_hexapole_full.iges",
]

def bbox(shape):
    b = Bnd_Box()
    BRepBndLib.Add_s(shape, b)
    xmin, ymin, zmin, xmax, ymax, zmax = b.Get()
    return (xmin, ymin, zmin, xmax, ymax, zmax)

def main():
    r = STEPControl_Reader()
    if r.ReadFile(STEP) != IFSelect_RetDone:
        print("ERROR: cannot read STEP"); sys.exit(1)
    r.TransferRoots()
    shape = r.OneShape()
    bb = bbox(shape)
    print("bbox after read (OCC applied file unit):")
    print("  X[%.2f, %.2f] Y[%.2f, %.2f] Z[%.2f, %.2f]" % (bb[0],bb[3],bb[1],bb[4],bb[2],bb[5]))
    span = max(bb[3]-bb[0], bb[4]-bb[1])
    # inch->mm 誤放大: 真實應 ~125mm; 若 span > 500 判定被 x25.4, scale 回
    if span > 500:
        t = gp_Trsf(); t.SetScaleFactor(1.0/25.4)
        shape = BRepBuilderAPI_Transform(shape, t, True).Shape()
        bb = bbox(shape)
        print("scaled x1/25.4 -> corrected to mm:")
        print("  X[%.2f, %.2f] Y[%.2f, %.2f] Z[%.2f, %.2f]" % (bb[0],bb[3],bb[1],bb[4],bb[2],bb[5]))
    else:
        print("no scaling (span=%.2f already mm-scale)" % span)
    # write IGES in mm, BRep mode (solids)
    Interface_Static.SetCVal_s("write.iges.unit", "MM")
    Interface_Static.SetIVal_s("write.iges.brep.mode", 1)
    for out in OUTS:
        w = IGESControl_Writer("MM", 1)
        w.AddShape(shape)
        ok = w.Write(out)
        print("wrote", out, "ok=" , bool(ok))

if __name__ == "__main__":
    main()
