# Fitting Data (`.mat` files)

MATLAB fitting results produced by scripts in `magnetic_sim/hung/analysis/`. Structure mirrors `magnetic_sim/hung/analysis/`:
- Root level → produced by `analysis/fit/` main pipeline
- `variants/` → produced by `analysis/variants/` alternative-version scripts, or legacy orphans

## Active pipeline (`data/` root)

| File | Producer | Consumers | Notes |
|------|----------|-----------|-------|
| `single_charge_ell.mat` | `analysis/fit/fit_ell_percoil.m` | (reference only) | Step 1: per-coil single-charge ell sweep |
| `KI_fit.mat` | `analysis/fit/fit_ell_perlayer.m` | fit_J, fit_J_50um, fit_B6x_1C, fit_B6x_6C, plot_J_quiver, plot_J_rmse, print_J_positions | Step 2: per-layer ell + sphere `pos` initial values. **Core hub — most downstream scripts load this.** |
| `J_ideal_fit.mat` | `analysis/fit/fit_J.m` | fit_B6x_1C, fit_B6x_6C | Step 3: [J] joint 6-coil fit with ideal K_I (main result). Needed as init for B6x. |
| `J_idealKI_50um.mat` | `analysis/fit/fit_J_50um.m` | (reference only) | Step 3a: [J] restricted to ±50 μm cube (Long-compatible) |
| `B6x_hung_1C.mat` | `analysis/fit/fit_B6x_1C.m` | (reference only) | Step 4: [B-6x] with 1 shared C |
| `B6x_hung_6C.mat` | `analysis/fit/fit_B6x_6C.m` | (reference only) | Step 4a: [B-6x] with 6 per-pole C |
| `all6_bias_fit.mat` | `analysis/fit/fit_B6x_allcoil.m` | (reference only) | Step 4b: [B-6x] all-coil superposition |
| `charge_model_fit.mat` | **⚠ No Hung script produces this** | `fit_B6x_allcoil.m`（**required**） | Method [A] baseline ell. Must exist for `fit_B6x_allcoil.m` to run. If regenerating, a Hung-specific method [A] script is needed (not yet written). |

## Pipeline execution order

```matlab
cd magnetic_sim/hung/analysis/fit
run fit_ell_percoil         % → single_charge_ell.mat
run fit_ell_perlayer        % → KI_fit.mat
run fit_J                   % → J_ideal_fit.mat      (needed by B6x_1C/6C)
run fit_B6x_1C              % → B6x_hung_1C.mat
run fit_B6x_6C              % → B6x_hung_6C.mat
run fit_B6x_allcoil         % → all6_bias_fit.mat    (needs charge_model_fit.mat beforehand)
```

## `variants/` — alternative / legacy versions

Not in the active pipeline. Kept available; re-runnable from `magnetic_sim/hung/analysis/variants/` if needed.

| File | Source | Status |
|------|--------|--------|
| `joint_6coil_fit.mat` | `analysis/variants/fit_J_fittedKI.m` | [J] with **fitted** K_I (superseded by fit_J.m which uses ideal K_I) |
| `KI_fit_v1.mat` | `analysis/variants/fit_KI_v1.m` | v1 K_I fit (single ell); superseded by `fit_ell_perlayer.m` (per-layer ell) |
| `B6x_hung_fit.mat` | **orphan** | Pre-split residual from before `fit_B6x_hung.m` was split into `fit_B6x_1C.m` and `fit_B6x_6C.m`. No longer produced or consumed. |
| `wp_fitting_data.mat` | **orphan** | No script in the repo references this file. Origin unknown; kept for safety. |

## Regenerating a file

All active files are regenerable by running the corresponding script in `magnetic_sim/hung/analysis/fit/`, assuming the upstream prerequisites exist (see pipeline order).

`charge_model_fit.mat` is the one exception — there is no Hung script that produces it yet. If it's ever deleted and needs to be rebuilt, a Hung method [A] script must be written (similar to `magnetic_sim/hexapole-long2016/analysis/fit_charge_model.m`).
