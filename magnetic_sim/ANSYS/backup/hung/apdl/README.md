# APDL Scripts

## Directory Structure

```
magnetic_sim/ANSYS/backup/hung/apdl/
├── geom/        Geometry build + IGES export (D-shape + 40 um fillet, main design)
├── sim/         6-coil simulation scripts (main D+fillet design)
├── postproc/    Post-processing: extract/export B-field data, generate figures
└── variants/    Alternative pole designs (baseline D-shape, full-round, round+fillet).
                 Not in the current pipeline, but kept fully functional for future use.
```

## geom/ — Geometry & IGES Export

| File | Description |
|------|-------------|
| `MT_Hung_Assembly_Dfillet.txt` | Build full hexapole (D-shape + 40 um fillet) and export `Full_Assembly_filleted.iges`. **Main geometry script.** Method A: smooth tangent fillet, cone semi-angle 11.31°, junction at 15.793 mm |
| `export_pole_filleted.txt` | Export single D+fillet pole as `Mag_Pole_Bottom_filleted.iges` |
| `export_parts.txt` | Export individual parts as separate IGES files to `magnetic_sim/ANSYS/backup/hung/IGES/`. Dimensions are hardcoded — update values when part sizes change, then re-run |

## sim/ — Main Simulation (D + fillet, 6 coils)

All 6 scripts share identical geometry, materials, mesh, and BC. Only `CURR_ARRAY` and output `/CWD` differ.

| File | Active Coil | CURR_ARRAY |
|------|------------|------------|
| `MT_Hung_Simulate_Coil1_filleted.txt` | P1 | [1,0,0,0,0,0] |
| `MT_Hung_Simulate_Coil2_filleted.txt` | P2 | [0,1,0,0,0,0] |
| `MT_Hung_Simulate_Coil3_filleted.txt` | P3 | [0,0,1,0,0,0] |
| `MT_Hung_Simulate_Coil4_filleted.txt` | P4 | [0,0,0,1,0,0] |
| `MT_Hung_Simulate_Coil5_filleted.txt` | P5 | [0,0,0,0,1,0] |
| `MT_Hung_Simulate_Coil6_filleted.txt` | P6 | [0,0,0,0,0,1] |

### Coil Settings (shared)
- Position: block top + COIL_DZ/2 (touching block upper face)
- Winding: clockwise (N1/N2 swapped, flux toward block)
- R_mean = 11 mm, DY = 2 mm, DZ = 15 mm, TURNS = 70

## postproc/ — Post-Processing

| File | Description | Used by |
|------|-------------|---------|
| `post_export_data.txt` | Export full-model + WP-region coordinate and B-field `.dat` files (Coil1) | `scripts/run/run_coil2to6_refined_v2.sh`, fitting pipeline |
| `post_export_data_coil[2-6].txt` | Same, coil-specific for Coils 2-6 | `scripts/run/run_post2to6.sh`, fitting pipeline |
| `post_extract_wp.txt` | Extract BX/BY/BZ/BSUM at WP (origin) | manual, quick WP measurement |
| `post_plot_geometry.txt` | Generate geometry PNG views (tip region, isometric, closeups) | manual, documentation figures |
| `post_plot_model.txt` | Generate PNG model views (isometric, top, front) | manual |
| `post_trace_tips.txt` | Extract BSUM at all 6 pole tips + WP center | manual, diagnostic |
| `post_trace_circuit.txt` | Extract BSUM at 8 points along P1 magnetic circuit | manual, diagnostic |

## variants/ — Alternative Designs

Not used in the current fitting pipeline, but kept working (filenames unchanged, no internal path issues) so a future study can revive any of them.

| Variant | Files | Status |
|---------|-------|--------|
| **Baseline D-shape (no fillet)** | `MT_Hung_Simulate_Coil[1-6].txt` (6 files) + `MT_Hung_SphereModel.txt` | All 6 coils available, can run a complete baseline sweep |
| **Full-round (no fillet)** | `MT_Hung_Simulate_Coil1_filled.txt` + `MT_Hung_SphereModel_filled.txt` + `export_poles_filled.txt` | Only Coil1 exists; to do a 6-coil study you'd need to create `_filled` variants for Coil2-6 |
| **Full-round + fillet** | `MT_Hung_Simulate_Coil1_round_filleted.txt` | Only Coil1 exists |
| **Steel-only diagnostic mesh** | `mesh_steel_only.txt` | Builds pole mesh without the air domain; useful for mesh QA |

To enable a variant, run the script directly — all use the same post-processing pipeline in `postproc/`.

## Running the Main Pipeline

```bash
# From project root, paths assume magnetic_sim/ANSYS/backup/hung/ as working tree
MAPDL="C:\Program Files\ANSYS2025R2\v252\ansys\bin\winx64\MAPDL.exe"
BASE="$PWD/hung"

# Simulation (example: Coil1)
"$MAPDL" -b -np 4 -m 8000 \
  -dir "$BASE/results/coil1/filleted" -j "coil1" \
  -i "$BASE/apdl/sim/MT_Hung_Simulate_Coil1_filleted.txt" \
  -o "$BASE/results/coil1/filleted/solve.out"

# Post-processing: export .dat for MATLAB fitting
"$MAPDL" -b -np 1 -m 4000 \
  -dir "$BASE/results/coil1/filleted" -j "coil1" \
  -i "$BASE/apdl/postproc/post_export_data.txt" \
  -o "$BASE/results/coil1/filleted/post.out"

# Quick WP measurement only
"$MAPDL" -b -np 1 -m 4000 \
  -dir "$BASE/results/coil1/filleted" -j "coil1" \
  -i "$BASE/apdl/postproc/post_extract_wp.txt" \
  -o "$BASE/results/coil1/filleted/wp_extract.out"
```

For batch-running Coils 2-6, use the shell scripts in `magnetic_sim/ANSYS/backup/hung/scripts/run/`:
- `run_coil2to6_refined_v2.sh` — **latest**, 3-pass NREFINE + POST1 export (recommended)
- `run_post2to6.sh` — re-export data only (reuses existing `sim1.db`)

See `magnetic_sim/ANSYS/backup/hung/scripts/run/README.md` for details and `magnetic_sim/ANSYS/backup/hung/scripts/run/variants/` for earlier versions.

## Legacy → New Path Mapping

For anyone looking up old paths referenced in git history or older docs:

| Old | New |
|-----|-----|
| `apdl/MT_Hung_SphereModel_filleted.txt` | `apdl/geom/MT_Hung_Assembly_Dfillet.txt` (**renamed**) |
| `apdl/MT_Hung_Simulate_Coil[1-6]_filleted.txt` | `apdl/sim/MT_Hung_Simulate_Coil[1-6]_filleted.txt` |
| `apdl/export_pole_filleted.txt` | `apdl/geom/export_pole_filleted.txt` |
| `apdl/export_parts.txt` | `apdl/geom/export_parts.txt` |
| `apdl/post_*.txt` | `apdl/postproc/post_*.txt` |
| `apdl/MT_Hung_Simulate_Coil[1-6].txt` (baseline) | `apdl/variants/…` |
| `apdl/MT_Hung_Simulate_Coil1_filled.txt` | `apdl/variants/…` |
| `apdl/MT_Hung_Simulate_Coil1_round_filleted.txt` | `apdl/variants/…` |
| `apdl/MT_Hung_SphereModel.txt`, `_filled.txt` | `apdl/variants/…` |
| `apdl/export_poles_filled.txt` | `apdl/variants/…` |
| `apdl/mesh_steel_only.txt` | `apdl/variants/…` |
