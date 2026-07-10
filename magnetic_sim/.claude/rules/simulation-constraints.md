---
globs: ["**/apdl/*.txt", "**/apdl/**"]
---

# Hexapole Simulation Constraints (Auto-loaded for APDL files)

## Before ANY geometry edit, verify these constraints:

1. **alpha = 54.74 degrees is FIXED** — arctan(sqrt(2)), derived from orthogonality + sphere + 60-degree offset. NEVER change this angle.
2. **Tip position formulas are locked:**
   - `R_norm_xy = R_norm * sqrt(2.0/3.0)`
   - `R_norm_z  = R_norm / sqrt(3.0)`
   - These follow directly from alpha. Only R_norm is a free parameter.
3. **Lower poles: 0, 120, 240 degrees. Upper poles: 60, 180, 300 degrees.** NEVER change this arrangement.
4. **3 pair axes must be orthogonal.** If you add or move poles, verify dot products = 0.

## Before ANY material/element edit:

5. **murx must be a constant** (not B-H curve) if linear superposition is used.
6. **Element types: SOLID96 (volumes) + SOURC36 (coils).** Do not change without user approval.

## Before ANY solver/BC edit:

7. **D,ALL,MAG,0 must exist** on all outer air domain surfaces before /SOLU.
8. **magsolv,3** (DSP method) requires the MAG=0 boundary for a unique solution.

## SOURC36 coil edits:

9. **Each coil = 3 nodes + Real Constants.** Do not create volume geometry for coils.
10. **Only CURR_ARRAY differs between coil scripts.** All other content must be synchronized.

## General:

- NEVER modify geometric parameters without explicit user approval
- Mark all changes with `[ADDED]` or `[MODIFIED]` comments
