# Troubleshooting

## Negative Pivot at Node 316044
- **Symptom:** `*** WARNING: Negative pivot at node 316044 MAG`
- **Cause:** Missing or insufficient boundary conditions for DSP method
- **Fix:** Apply `D,ALL,MAG,0` on all outer air cylinder boundary surfaces (lateral + top/bottom caps)
- **Notes:** The DSP method (magsolv,3) requires explicit boundary conditions to obtain a unique solution

## Coefficient Ratio > 1e8
- **Symptom:** `Coefficient ratio exceeds 1.0e8`
- **Cause:** Extreme element size ratio: pole tip (40 um) vs outer air (80 mm) = 2000:1
- **Fix:** This is expected for this geometry. The SmartSize algorithm handles the transition. Monitor solution convergence but no action needed.

## SmartSizing Small Angle Warning
- **Symptom:** `Small angle detected in area near pole tip, mesh coarsened`
- **Cause:** Sharp cone geometry at pole tips (40 um radius fillet on 3 mm base cone)
- **Fix:** Expected behavior. ANSYS auto-coarsens to avoid degenerate elements. Does not affect solution accuracy significantly.

## "No Constraints" Warning
- **Symptom:** `No constraints have been defined` during /SOLU
- **Cause:** Forgot to apply MAG=0 boundary condition before entering /SOLU
- **Fix:** Exit /SOLU, apply boundary conditions in /PREP7, then re-enter /SOLU
- **Prevention:** Always verify the `[ADDED] Apply far-field boundary condition` block exists before magsolv

## Memory Exceeded
- **Symptom:** `Insufficient memory` or ANSYS crashes
- **Cause:** Large mesh in outer air sphere (~2.4M elements at SMRT=5)
- **Fix:** Use `-m 24000` flag (24 GB) when launching batch mode. If still failing, increase SMRT level (coarser mesh) or reduce AIR_CYL dimensions.
