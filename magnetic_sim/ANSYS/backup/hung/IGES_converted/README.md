# IGES Converted — Unit Flag Fixed

## Description
Same parts as `magnetic_sim/ANSYS/backup/hung/IGES/`, but with IGES unit flag corrected from 6 (mm) to 1 (inches).
This ensures reliable conversion when saving as STEP in SolidWorks or other CAD software.

## Difference from `magnetic_sim/ANSYS/backup/hung/IGES/`
- `magnetic_sim/ANSYS/backup/hung/IGES/` — raw ANSYS export, unit flag = 6 (mm), values in inches. SolidWorks opens correctly but STEP conversion may fail.
- `magnetic_sim/ANSYS/backup/hung/IGES_converted/` — unit flag fixed to 1 (inches). SolidWorks opens correctly AND STEP conversion preserves dimensions.

## Files

| File | Description |
|------|-------------|
| `Mag_Pole_Bottom.iges` | Magnetic pole (D-shape + cone) |
| `Pole_Block_Top.iges` | Upper block (L-shape) |
| `Pole_Block_Bottom.iges` | Lower block (L-shape) |
| `Mag_Guide_Post.iges` | Guide post (R=4mm, H=46mm) |
| `Coil.iges` | Excitation coil (R_in=10, R_out=12, H=15mm) |
| `Upper_Ring.iges` | Yoke ring (R_in=38, R_out=62.5, T=2mm) |
| `Full_Assembly.iges` | Complete hexapole assembly |
| `Mag_Pole_Bottom_filleted.iges` | Pole with 40 um diameter tip fillet (single part) |
| `Full_Assembly_filleted.iges` | Complete assembly with 40 um tip fillet on all 6 poles |

## Units
All files: values in inches, unit flag = 1 (inches).
SolidWorks reads as inches x 25.4 = correct mm dimensions.

## Regenerating
When parts change, re-export from `magnetic_sim/ANSYS/backup/hung/IGES/` and fix unit flag:
```bash
cp magnetic_sim/ANSYS/backup/hung/IGES/Part.iges magnetic_sim/ANSYS/backup/hung/IGES_converted/Part.iges
sed -i "s/,1.0,6,,/,1.0,1,,/" magnetic_sim/ANSYS/backup/hung/IGES_converted/Part.iges
```
