# Hung Hexapole — Figures

## Folder Structure

```
figures/
├── coil1/ ~ coil6/    B-field distribution plots per single-coil excitation
│   └── fail/          Rejected or superseded versions of those plots
├── analytic/          Analysis figures: fitting, geometry, sensor sensitivity
└── README.md          This file
```

**Convention** (per project figure-organization rule):
- Successful distribution plots go in the `coilN/` root and **overwrite the previous file with the same name**.
- Rejected / superseded versions move to `coilN/fail/`.
- Analysis, parameter-sweep, geometry-annotation, and sensor-study figures live in `analytic/` (not coil-specific).

Currently only **coil1/** has populated B-field plots; coil2–6 have empty `fail/` subfolders reserved.

## coil1/ — B-field plots (D-shape + 40 µm fillet, 70 A-turns)

| File | Description |
|------|-------------|
| `Bvector_topview.png` | Top-view B-field vectors over the full model (±80 mm, Tesla) |
| `Bvector_topview_mT.png` | Same top view, rescaled to mT (close-up near working point) |
| `Bvector_3D_iron.png` | 3D B-field arrows inside iron only, showing flux circulation |
| `Bcontour_P1_topview.png` | \|B\| contour at P1 tip top view |
| `Bcontour_xaza_Dfillet.png` | \|B\| contour on x_a-z_a plane (y_a=0), D+fillet 3-pass NREFINE data, 2D slab interpolation. Actuator frame: x_a→P1, y_a→P3, z_a→P5 |
| `Bcontour_xaza_Dfillet_smooth.png` | Same view, upgraded with 3D interpolation + Gaussian smoothing; clearest production version (WP = 8.74 mT, 526,645 nodes) |
| `Bcontour_xaza_Dshape_baseline.png` | Same view, D-shape **no fillet / no NREFINE** baseline (368,686 nodes, WP = 9.31 mT) — shown as convergence contrast to `_smooth` |

## analytic/ — Fitting & analysis figures

### Point-charge model fitting
| File | Description |
|------|-------------|
| `charge_model_cost_landscape.png` | Method [A] cost vs ell (1D search landscape) |
| `charge_positions_J.png` | 6 fitted charge positions from Method [J] ideal K_I |
| `charge_projection_lower_upper.png` | Lower/Upper charge projection schematic (s, θ, l relationship) |
| `fitting_J_idealKI_quiver_error.png` | [J] ideal K_I: FEM vs model quiver + % error scatter |
| `fitting_J_idealKI_RMSE.png` | [J] ideal K_I: FEM vs model quiver + RMSE (mT) |
| `fitting_J_idealKI_50um_RMSE.png` | [J] ±50 µm cube (Long-compatible) + RMSE |
| `fitting_B6x_quiver_error.png` | [B-6x] alternating superposition: FEM vs model quiver + error |

### Direction vectors (d̂)
| File | Description |
|------|-------------|
| `dhat_vs_pole_axis.png` | d̂ (WP→charge) vs pole axis direction, Hung |
| `dhat_vs_pole_axis_Long.png` | Same comparison for Long 2016 (reference) |
| `dhat_3D_Long.png` | 3D view of Long 2016 d̂ vectors |

### Hall-sensor analysis
| File | Description |
|------|-------------|
| `sensor_sensitivity_full.png` | Sensor sensitivity vs position, full range with installable annotation |
| `sensor_sensitivity_vs_position.png` | Same, zoomed to installable range |
| `sensor_ratio_fit.png` | Sensitivity ratio fit curve |
| `sensor_ratio_vs_distance.png` | Sensitivity ratio vs sensor-to-WP distance |
| `sensor_Btop_full_range.png` | B-field magnitude at sensor top surface, full range |
| `sensor_Btop_vs_position.png` | Same, zoomed |

### Geometry annotations
| File | Description |
|------|-------------|
| `hexapole_tip_annotated.png` | 6 pole tips on R=0.5 mm sphere with P1–P6 labels and l annotation |
| `mesh_steel_filleted_3D.png` | 3D mesh visualization of steel body (D-shape + fillet, no yoke) |

### Mesh convergence
| File | Description |
|------|-------------|
| `WP_convergence_Hung.png` | WP \|B\| vs total mesh node count: baseline (no NREFINE, 368k, 9.31 mT) → 3-pass NREFINE (526k, 8.74 mT); Richardson asymptote ≈ 8.60 mT |

## Generating Figures

Paths reflect the current post-reorg structure:

```bash
# Step 1: Run ANSYS simulation (D-shape + fillet, coil1)
cd magnetic_sim/ANSYS/backup/hung/
"C:\Program Files\ANSYS2025R2\v252\ansys\bin\winx64\MAPDL.exe" -b -np 4 -m 8000 \
  -dir "results/coil1/filleted" -j "coil1" \
  -i "apdl/sim/MT_Hung_Simulate_Coil1_filleted.txt" \
  -o "results/coil1/filleted/solve.out"

# Step 2: Export B-field data for MATLAB
"C:\Program Files\ANSYS2025R2\v252\ansys\bin\winx64\MAPDL.exe" -b -np 1 -m 4000 \
  -dir "results/coil1/filleted" -j "coil1" \
  -i "apdl/postproc/post_export_data.txt"

# Step 3: Generate plots in MATLAB
cd analysis/plot
# In MATLAB:
#   run('plot_Bfield_2d.m')      -> figures/ root (move to figures/coil1/)
#   run('plot_Bfield_3d.m')      -> Bvector_3D_iron.png
#   run('plot_J_quiver.m')       -> figures/analytic/fitting_J_idealKI_quiver_error.png
#   run('plot_J_rmse.m')         -> figures/analytic/fitting_J_idealKI_RMSE.png
#   run('plot_charge_proj.m')    -> figures/analytic/charge_projection_lower_upper.png
#   run('plot_Bcontour_xaza.m')          -> figures/coil1/Bcontour_xaza_Dfillet_smooth.png
#   run('plot_Bcontour_xaza_baseline.m') -> figures/coil1/Bcontour_xaza_Dshape_baseline.png
```

Note: `plot_Bfield_2d.m` and `plot_Bfield_3d.m` save to the `figures/` root — manually move the output into the appropriate `coilN/` subfolder after each run.

## Simulation Parameters

| Parameter | Value |
|-----------|-------|
| Excitation | 70 A-turns (TURNS=70, I=1 A) |
| Steel mu_r | 280 (linear, constant) |
| Pole shape | D-shape + 40 µm diameter tip fillet |
| Mesh | tips fine → steel ~ coarse air |
| Solver | magsolv,3 (DSP) |
