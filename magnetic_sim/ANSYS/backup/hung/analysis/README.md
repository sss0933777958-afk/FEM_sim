# MATLAB Analysis Scripts

## Prerequisites
- Run ANSYS simulation first (`MT_Hung_Simulate_CoilN_filleted.txt`) for all 6 coils
- Run data export (`post_export_data_coilN.txt`) to produce `.dat` files in `results/coilN/filleted/`

## Directory Structure

```
magnetic_sim/ANSYS/backup/hung/analysis/
├── core/        Shared utilities (constants, data loader)
├── fit/         Charge-model fitting pipeline ([pre-J] → [J] → [B-6x])
├── plot/        Figure generation (B-field, quiver error, RMSE, projection)
├── util/        Small helpers (position printout)
└── variants/    Alternative or previous-version scripts (not actively used but may be revisited)
```

## Scripts

### core/ — shared utilities
| File | Description |
|------|-------------|
| `mt_constants.m` | Hung hexapole design constants: R_norm, tip positions, yoke dimensions, pole angles |
| `import_ansys_data.m` | Read ANSYS-exported coordinate + B-field `.dat` files. Handles MAPDL banner headers and concatenated negative numbers. Returns struct with node_id, x, y, z, bx, by, bz, bsum |

### fit/ — fitting pipeline
Run in this order. Each script is independent but later steps load outputs from earlier ones.

| Order | File | Output | Description |
|------:|------|--------|-------------|
| 1 | `fit_ell_percoil.m` | `single_charge_ell.mat` | Per-coil single-charge ell sweep (no K_I), gives rough |l| per pole |
| 2 | `fit_ell_perlayer.m` | `KI_fit.mat` | Per-layer ell averages + sphere `pos` initial values for [J] |
| 3 | `fit_J.m` | `J_ideal_fit.mat` | **[J] Joint 6-coil fit with ideal K_I** — main result (err ~0.46%). Needed as init by B6x_1C/6C. |
| 3a | `fit_J_50um.m` | `J_idealKI_50um.mat` + figure | [J] restricted to ±50 μm cube (Long-compatible) + RMSE figure |
| 4 | `fit_B6x_1C.m` | `B6x_hung_1C.mat` | [B-6x] 18 pos + 1 shared C (VarPro). Loads `J_ideal_fit.mat` and `KI_fit.mat`. |
| 4a | `fit_B6x_6C.m` | `B6x_hung_6C.mat` | [B-6x] 18 pos + 6 per-pole C (linear LS). Same loads as 4. |
| 4b | `fit_B6x_allcoil.m` | `all6_bias_fit.mat` | [B-6x] all-coil superposition, requires `charge_model_fit.mat` (no Hung producer — see `data/README.md`) |

### plot/ — figures
| File | Output | Description |
|------|--------|-------------|
| `plot_Bfield_2d.m` | `figures/Bvector_topview.png`, `Bfield_xy.png`, `Bfield_xz.png` | Long2016-style vector + contour plots |
| `plot_Bfield_3d.m` | `figures/Bfield_3d.png` | 3D B-field arrow plot (iron flux flow) |
| `plot_J_quiver.m` | `figures/analytic/fitting_J_idealKI_quiver_error.png` | [J] FEM vs Model quiver + % error |
| `plot_J_rmse.m` | `figures/analytic/fitting_J_idealKI_RMSE.png` | [J] FEM vs Model quiver + RMSE (mT) |
| `plot_charge_proj.m` | `figures/analytic/charge_projection_lower_upper.png` | Lower/Upper charge projection schematic |

### util/
| File | Description |
|------|-------------|
| `print_J_positions.m` | Print per-pole positions table from [J] ideal K_I fit |

### variants/ — alternative versions, kept available for future use
These are not part of the current main pipeline, but are kept intact (paths already adjusted) so you can re-run them when needed (e.g., to re-check a previous method against the current main one, or to explore an alternative).

| File | Current main | Context / when to use |
|------|--------------|----------------------|
| `fit_KI_v1.m` | `fit/fit_ell_perlayer.m` | v1 uses single ell from method [A]; v2 uses per-layer ell. Re-run v1 if you want to compare the single-ell path. |
| `fit_J_fittedKI.m` | `fit/fit_J.m` | Uses **fitted** K_I; main pipeline uses **ideal** K_I (0.46% error, confirmed sufficient). Re-run if you want to revisit non-ideal K_I behaviour. |

## Path Mechanism

Each active script adds `core/` to the MATLAB path automatically:

```matlab
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'core'));
```

This lets you `cd` into `fit/` (or `plot/`, `util/`) and run scripts directly without prior setup.

## Usage Example

```matlab
% Generate B-field figures for coil1
cd magnetic_sim/ANSYS/backup/hung/analysis/plot
run('plot_Bfield_2d.m')
run('plot_Bfield_3d.m')

% Run full fitting pipeline
cd ../fit
run('fit_ell_percoil.m')     % produces single_charge_ell.mat
run('fit_ell_perlayer.m')    % produces KI_fit.mat (needed by later steps)
run('fit_J.m')               % produces J_ideal_fit.mat (needed by B6x_1C/6C)
run('fit_B6x_1C.m')          % [B-6x] with shared C
run('fit_B6x_6C.m')          % [B-6x] with per-pole C

% Plot fitting error
cd ../plot
run('plot_J_quiver.m')
run('plot_J_rmse.m')
```

Output `.mat` files go to `magnetic_sim/ANSYS/backup/hung/data/`. Output `.png` files go to `magnetic_sim/ANSYS/backup/hung/figures/` (2D/3D B-field) or `magnetic_sim/ANSYS/backup/hung/figures/analytic/` (fitting error, projection).

## Data Format
The `.dat` files produced by `post_export_data_coilN.txt`:
- `coilN_coord_{all,wp}.dat` — node coordinates (NODE, X, Y, Z)
- `coilN_bfield_{all,wp}.dat` — B-field (NODE, BX, BY, BZ, BSUM)
- WP region = nodes within 2 mm of origin
- Units: SI (meters, Tesla)

## Figure Style (matching Long2016)
- Font: Helvetica, 12pt labels, 13pt titles
- Colormap: turbo(256)
- DPI: 300 (figures), 400 (3D)
- Contour: 20 levels, no line edges

## Legacy → New Filename Mapping

For anyone looking up old filenames referenced in git history or documentation:

| Old | New |
|-----|-----|
| `fit_single_charge.m` | `fit/fit_ell_percoil.m` |
| `fit_KI_v2.m` | `fit/fit_ell_perlayer.m` |
| `test_J_idealKI.m` | `fit/fit_J.m` |
| `fit_J_idealKI_50um.m` | `fit/fit_J_50um.m` |
| `fit_B6x_hung.m` | `fit/fit_B6x_1C.m` |
| `fit_B6x_perC.m` | `fit/fit_B6x_6C.m` |
| `fit_all6_with_bias.m` | `fit/fit_B6x_allcoil.m` |
| `generate_figures.m` | `plot/plot_Bfield_2d.m` |
| `generate_fig_3d.m` | `plot/plot_Bfield_3d.m` |
| `plot_J_idealKI.m` | `plot/plot_J_quiver.m` |
| `plot_J_idealKI_RMSE.m` | `plot/plot_J_rmse.m` |
| `plot_charge_projection.m` | `plot/plot_charge_proj.m` |
| `print_J_idealKI_positions.m` | `util/print_J_positions.m` |
| `fit_KI.m` | `variants/fit_KI_v1.m` |
| `fit_joint_6coil.m` | `variants/fit_J_fittedKI.m` |
| `plot_fitting_J.m` | *(deleted — obsolete figure also removed)* |
