# Hexapole Long 2016 Design

6-pole magnetic tweezers simulation based on Fei Long's 2016 PhD dissertation (Ohio State University).

## Design Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| Poles | 6 (3 upper + 3 lower) | Hexapole configuration |
| Pole tip radius | 40 um | Spherical fillet at tip |
| Workspace radius (rho) | 500 um | Physical workspace |
| Yoke inner/outer radius | 42 / 53 mm | Iron yoke ring |
| Coil turns (N_c) | 70 | Per coil |
| Material (murx) | 280 | 1018 steel, linear approx. |
| Mesh | SmartSize 5 | ~2.4M elements, ~494k nodes |
| Solver | DSP (magsolv,3) | Requires MAG=0 boundary |

## Pole Naming Convention

| Paper | APDL Index | Angle | Layer |
|-------|-----------|-------|-------|
| P1 | Coil 1 | 0 deg | Lower |
| P2 | Coil 5 | 180 deg | Upper |
| P3 | Coil 2 | 120 deg | Lower |
| P4 | Coil 6 | 300 deg | Upper |
| P5 | Coil 4 | 60 deg | Upper |
| P6 | Coil 3 | 240 deg | Lower |

## Simulation Approach

Linear superposition: 6 independent unit-excitation solves (1A through one coil, 0A for rest). Any multi-coil excitation is a linear combination of these 6 basis solutions.

## Key Results

- WP center field: ~8.7 mT per coil (unit excitation)
- [B-6x] fitting error: 0.07% (19 parameters)
- R_a: ~1.03e9 A/Wb (vs dissertation 6.3e8, difference due to murx=280)

## File Structure

- `apdl/` — 6 coil scripts + 1 variant (Coil5_sph10) + 7 POST1 extraction scripts
- `analysis/` — MATLAB: utilities, fitting ([A]/[J]/[B-6x]), verification, figures
- `data/` — Fitting results (.mat)
- `figures/` — Publication figures (.png)
- `docs/` — Technical documentation (fitting methods, notation, troubleshooting, etc.)
