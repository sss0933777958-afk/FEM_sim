"""Generate gap200um geom export scripts (mm + metre) from halfcut baselines.

Pipeline mirrors run_geom_export.ps1 conventions:
  - mm script    -> IGES_converted/.../Geom_gap200um.iges (units flag patched 6->2)
  - metre script -> IGES/.../Geom_gap200um.iges            (units flag 6 = ANSYS default, FEM reload source of truth)

This script only WRITES the two APDL files; ANSYS run + units-flag patch is
done by run_geom_export_gap200um.ps1 (one-shot wrapper).

Patch applied: same protrusion CYL4 -> 2-stacked split as
_generate_halfcut_sims_gap.py.

Run from any shell:
    python _generate_geom_export_gap.py
"""
from pathlib import Path

REPO = Path(r"G:\my_workspace\code\FEM_sim")
GEOM_DIR = REPO / "kuo" / "apdl" / "geom" / "long2016_hexapole_halfcut"


# ---------------------------------------------------------------------
# mm version (PROT_H = 7.0 raw mm; AIR_GAP = 0.2 mm)
# ---------------------------------------------------------------------
PROT_FROM_MM = """*do, prot_ind, 1, POLE_NUM,1
    WPLANE,1,0,0,0,,,,,,
    WPcsys,-1,0
    *if, prot_ind,LE,3, then
        temp_prot_h = -PROT_H
        temp_prot_ang = 120*(prot_ind - 1)
    *else
        WPOFFS,0,0,YOKE_H
        temp_prot_h = PROT_H
        temp_prot_ang = 120*(prot_ind - 4) + 60
    *endif
    CYL4,YOKE_MID_R*cos(temp_prot_ang),YOKE_MID_R*sin(temp_prot_ang),PROT_R, ,,,temp_prot_h
*enddo"""

PROT_TO_MM = """! [ADDED gap200um] 200 um air slab (raw mm units; PROT_H = 7.0 mm)
AIR_GAP = 0.2
PROT_STEEL_H = PROT_H - AIR_GAP

! [ADDED gap200um] Each protrusion = 2 stacked CYL4 (yoke-side air slab + steel)
*do, prot_ind, 1, POLE_NUM,1
    WPLANE,1,0,0,0,,,,,,
    WPcsys,-1,0
    *if, prot_ind,LE,3, then
        temp_prot_ang = 120*(prot_ind - 1)
        xc = YOKE_MID_R*cos(temp_prot_ang)
        yc = YOKE_MID_R*sin(temp_prot_ang)
        CYL4, xc, yc, PROT_R, 0, , , -AIR_GAP        ! yoke-side air slab
        WPOFFS, 0, 0, -AIR_GAP
        CYL4, xc, yc, PROT_R, 0, , , -PROT_STEEL_H   ! steel main
    *else
        temp_prot_ang = 120*(prot_ind - 4) + 60
        xc = YOKE_MID_R*cos(temp_prot_ang)
        yc = YOKE_MID_R*sin(temp_prot_ang)
        WPOFFS, 0, 0, YOKE_H
        CYL4, xc, yc, PROT_R, 0, , , AIR_GAP         ! yoke-side air slab
        WPOFFS, 0, 0, AIR_GAP
        CYL4, xc, yc, PROT_R, 0, , , PROT_STEEL_H    ! steel main
    *endif
*enddo"""


# ---------------------------------------------------------------------
# metre version (PROT_H = 7.0e-3 SI; AIR_GAP = 200e-6 = 200 um SI)
# ---------------------------------------------------------------------
PROT_FROM_METRE = PROT_FROM_MM  # source uses same protrusion block structure

PROT_TO_METRE = """! [ADDED gap200um] 200 um air slab (SI metres; PROT_H = 7.0e-3 m)
AIR_GAP = 200e-6
PROT_STEEL_H = PROT_H - AIR_GAP

! [ADDED gap200um] Each protrusion = 2 stacked CYL4 (yoke-side air slab + steel)
*do, prot_ind, 1, POLE_NUM,1
    WPLANE,1,0,0,0,,,,,,
    WPcsys,-1,0
    *if, prot_ind,LE,3, then
        temp_prot_ang = 120*(prot_ind - 1)
        xc = YOKE_MID_R*cos(temp_prot_ang)
        yc = YOKE_MID_R*sin(temp_prot_ang)
        CYL4, xc, yc, PROT_R, 0, , , -AIR_GAP        ! yoke-side air slab
        WPOFFS, 0, 0, -AIR_GAP
        CYL4, xc, yc, PROT_R, 0, , , -PROT_STEEL_H   ! steel main
    *else
        temp_prot_ang = 120*(prot_ind - 4) + 60
        xc = YOKE_MID_R*cos(temp_prot_ang)
        yc = YOKE_MID_R*sin(temp_prot_ang)
        WPOFFS, 0, 0, YOKE_H
        CYL4, xc, yc, PROT_R, 0, , , AIR_GAP         ! yoke-side air slab
        WPOFFS, 0, 0, AIR_GAP
        CYL4, xc, yc, PROT_R, 0, , , PROT_STEEL_H    ! steel main
    *endif
*enddo"""


JOBS = [
    # (src, dst, prot_from, prot_to, cwd_from, cwd_to, iges_from, iges_to)
    (
        "MT_Geom_Export_mm.txt",
        "MT_Geom_Export_mm_gap200um.txt",
        PROT_FROM_MM, PROT_TO_MM,
        r"\geom_export_mm'", r"\geom_export_mm_gap200um'",
        r"\Long2016_HexapoleHalfcut_Geom', 'iges'",
        r"\Long2016_HexapoleHalfcut_Geom_gap200um', 'iges'",
    ),
    (
        "MT_Geom_Export.txt",
        "MT_Geom_Export_gap200um.txt",
        PROT_FROM_METRE, PROT_TO_METRE,
        r"\geom_export_metre'", r"\geom_export_metre_gap200um'",
        r"\Long2016_HexapoleHalfcut_Geom', 'iges'",
        r"\Long2016_HexapoleHalfcut_Geom_gap200um', 'iges'",
    ),
]


def patch(src_name, dst_name, prot_from, prot_to, cwd_from, cwd_to, iges_from, iges_to):
    src = GEOM_DIR / src_name
    dst = GEOM_DIR / dst_name
    if not src.exists():
        raise FileNotFoundError(src)
    text = src.read_text(encoding="utf-8")

    if prot_from not in text:
        raise RuntimeError(f"{src_name}: PROT_FROM block not found")
    text = text.replace(prot_from, prot_to)

    if cwd_from not in text:
        raise RuntimeError(f"{src_name}: /CWD marker {cwd_from!r} not found")
    text = text.replace(cwd_from, cwd_to)

    if iges_from not in text:
        raise RuntimeError(f"{src_name}: IGESOUT marker {iges_from!r} not found")
    text = text.replace(iges_from, iges_to)

    dst.write_text(text, encoding="utf-8")
    print(f"Wrote {dst.name}")


if __name__ == "__main__":
    for j in JOBS:
        patch(*j)
