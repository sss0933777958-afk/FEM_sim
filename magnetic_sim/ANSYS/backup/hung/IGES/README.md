# IGES Part Exports

## Description
Individual part IGES files exported from ANSYS with SolidWorks-compatible units (MM=1/25.4).
SolidWorks reads values as inches x 25.4 = correct mm dimensions.

## Files

| File | Part | Key Dimensions |
|------|------|---------------|
| `Mag_Pole_Bottom.iges` | Magnetic pole (D-shape + cone) | R=3.175mm, cone=15.875mm, total=43mm, flat=28mm |
| `Pole_Block_Top.iges` | Upper block (L-shape) | 25x22x10mm, penetration=7mm |
| `Pole_Block_Bottom.iges` | Lower block (L-shape) | 22x22x10mm, penetration=4.5mm |
| `Mag_Guide_Post.iges` | Guide post cylinder | R=4mm, H=46mm |
| `Coil.iges` | Excitation coil ring | R_in=10mm, R_out=12mm, H=15mm |
| `Upper_Ring.iges` | Yoke (iron ring) | R_in=38mm, R_out=62.5mm, T=2mm |
| `Full_Assembly.iges` | Complete hexapole assembly | All parts assembled, unit flag=1 (inches) |
| `Mag_Pole_Bottom_filleted.iges` | Pole with 40 um diameter tip fillet (single part) | R_tip=20um, semi-angle 11.31°, junction at 15.793mm |
| `Full_Assembly_filleted.iges` | Complete assembly with 40 um tip fillet on all 6 poles | Same as above for all poles |

## STEP Conversion
To convert to STEP for other CAD software:
1. Open IGES in SolidWorks (dimensions will be correct mm)
2. File → Save As → .STEP

`Full_Assembly.iges` has unit flag fixed to 1 (inches) for reliable STEP conversion.

## Regenerating
When part dimensions change:
1. Update values in `magnetic_sim/hung/apdl/geom/export_parts.txt`
2. Run: `MAPDL -b -np 1 -m 2000 -dir "IGES" -j "parts" -i "apdl/geom/export_parts.txt"`
3. For full assembly: run `apdl/variants/MT_Hung_SphereModel.txt`, copy output to `Full_Assembly.iges`, fix unit flag with `sed -i "s/,1.0,6,,/,1.0,1,,/" Full_Assembly.iges`

## Filleted Versions (40 um diameter tip fillet)

The `_filleted.iges` versions add a 40 um diameter (20 um radius) spherical fillet at each pole tip,
replacing the mathematically sharp point with a smooth rounded tip.

**Why fillet?** Real machined parts cannot have R=0 sharp tips. EDM machining typically achieves
20-50 um. Adding a fillet matches reality and improves FEM mesh stability at the tip.

**Geometry (Method A — smooth tangent):**
- Foremost point of fillet at original tip position (preserves R=0.5mm sphere constraint)
- Cone semi-angle unchanged at 11.31° (matches original Hung STEP)
- Cone-cyl junction shifts back from 15.875 -> 15.793 mm (delta = 0.082 mm)
- Total pole length unchanged at 43 mm
- All other features (flat at 28mm, far end at 43mm, POLE_R=3.175mm) unchanged

**Constraints preserved:**
- All 6 fillet foremost points still on R=0.5mm sphere
- Orthogonal pair axes still 90 degrees
- tip-to-tip distance still 1.0 mm
- Magic angle 54.74 deg unchanged

**To regenerate:**
- Single pole: `MAPDL ... -i apdl/geom/export_pole_filleted.txt`
- Full assembly: `MAPDL ... -i apdl/geom/MT_Hung_Assembly_Dfillet.txt`
