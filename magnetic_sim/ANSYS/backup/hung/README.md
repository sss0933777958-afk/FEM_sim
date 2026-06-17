# Hung Hexapole Magnetic Tweezers — FEM Simulation

## Overview
ANSYS APDL magnetostatic simulation of the Hung hexapole magnetic tweezers design.
6 tilted poles on a R=0.5mm sphere at magic angle (54.74°), with 1018 steel (mu_r=280).

## Folder Structure

```
magnetic_sim/ANSYS/backup/hung/
├── apdl/              APDL scripts
│   ├── MT_Hung_Simulate_Coil[1-6].txt        Original D-shape poles
│   ├── MT_Hung_Simulate_Coil1_filleted.txt    D-shape + R=40um fillet
│   ├── MT_Hung_Simulate_Coil1_filled.txt      Full-round sharp tip
│   ├── MT_Hung_Simulate_Coil1_round_filleted.txt  Full-round + R=40um fillet
│   ├── MT_Hung_SphereModel*.txt               IGES export scripts
│   ├── mesh_steel_only.txt                    Steel-only mesh (no air domain)
│   └── post_*.txt                             Post-processing scripts
├── analysis/          MATLAB scripts: data import, figure generation
├── comsol/            COMSOL models (gitignored)
├── docs/              Technical documentation
├── figures/
│   ├── analytic/      Analysis figures (sensor, fitting, parameter sweeps)
│   ├── geometry/      3D geometry annotations
│   ├── coil1/         Coil 1 B-field distribution figures
│   │   └── fail/      Rejected figures
│   └── coil2~6/       (empty, for future coil simulations)
│       └── fail/
├── IGES/              IGES exports: raw ANSYS output (unit flag=6)
├── IGES_converted/    IGES exports: unit flag=1 (SolidWorks compatible)
├── results/           ANSYS outputs (gitignored, ~10GB each)
└── scripts/           Utility scripts
```

## Pole Variants

| Script | Pole Shape | Tip | Status |
|--------|-----------|-----|--------|
| `Simulate_Coil1.txt` | D-shape | Sharp | Baseline |
| `Simulate_Coil1_filleted.txt` | D-shape + R=40um fillet | Rounded | Main model |
| `Simulate_Coil1_filled.txt` | Full-round | Sharp | Tested |
| `Simulate_Coil1_round_filleted.txt` | Full-round + R=40um fillet | Rounded | Tested |

## Design Parameters

| Parameter | Value |
|-----------|-------|
| R_sphere | 0.5 mm |
| TILT_UP (P2,P4,P5) | 35° |
| TILT_DN (P1,P3,P6) | 5.71° |
| Pole total length | 43.0 mm |
| Pole shaft radius | 3.175 mm |
| Pole cone length | 15.875 mm |
| Tip fillet radius | 40 um |
| Coil R_in / R_out | 10 / 12 mm |
| Coil height | 15 mm |
| Turns | 70 |
| Steel mu_r | 280 |

## Key Results (Coil1, 70 A-turns, D-shape + fillet)

| Location | |B| |
|----------|-----|
| WP center | 9.66 mT |
| P1 tip surface | 424 mT |
| Iron max (guide post) | 1493 mT |
| P1 surface (s=15.75mm) | 2.82 mT |

## Quick Start

```bash
cd magnetic_sim/ANSYS/backup/hung

# Run simulation (D-shape + fillet, coil 1)
MAPDL -b -np 2 -m 8000 -dir "results/coil1/filleted" -j "coil1" \
  -i "apdl/sim/MT_Hung_Simulate_Coil1_filleted.txt" \
  -o "results/coil1/filleted/solve.out"

# Export B-field data
MAPDL -b -np 1 -m 4000 -dir "results/coil1/filleted" -j "coil1" \
  -i "apdl/postproc/post_export_data.txt"

# Generate figures (MATLAB)
run('analysis/generate_figures.m')
```
