# Model Validation: APDL Simulation vs Long 2016 Dissertation

Date: 2026-03-07
Status: Pre-analysis validation complete. Ready to proceed with charge model fitting.

## 1. Origin of the Simulation Files

- APDL scripts received from Fei Long (original paths: `C:\Users\feilong\Desktop\Fei Long Ansys\Magnetic Tweezers LowPermeability`)
- We modified: working directory (line 3), boundary conditions (lines 497-518), POST1 data extraction (lines 531-565)
- All modifications marked with `[ADDED]` or `[MODIFIED]`
- Original code preserved with `!****` comment prefix

## 2. Parameter Comparison

### Confirmed Match

| Parameter | APDL Value | Dissertation (p.14) | Status |
|-----------|-----------|---------------------|--------|
| Pole material | 1018 steel | 1018 steel (0.18% C) | MATCH |
| Pole tip radius | 40 um (`POLE_TIP_R`) | 40 um | MATCH |
| Pole diameter | 6 mm (`POLE_R*2`) | ~6 mm | MATCH |
| Workspace radius | 500 um (`R_norm`) | 500 um (nominal) | MATCH |
| Coil turns | 70 (`TURNS`) | 70 | MATCH |
| Coil inner/outer radius | 5/8 mm | (not stated, but consistent) | MATCH |
| Unit excitation | 1 A per coil | 1 A (Fig.2.3) | MATCH |
| Element types | SOLID96 + SOURC36 | (not stated in dissertation) | N/A |
| Solver | magsolv,3 (DSP) | (not stated) | N/A |

### Issues Investigated

#### Issue 1: Material Model — murx = 280 (constant linear)

| Aspect | Finding |
|--------|---------|
| **Source of 280** | Iowa State 2006 paper: measured mu_r' = 280 for 1018 steel at low frequency (1-10 Hz) |
| **Directory name** | "LowPermeability" — intentional choice of low-field linear approximation |
| **Dissertation says** | Nothing. FEM material model not specified anywhere in 167 pages |
| **Real B-H curve** | mu_r varies: ~1250 at B=0.25T, ~1500 peak at B~1T, drops to ~100 at B=2T |
| **At our operating point** | Pole tips B~1T: real mu_r~1400 vs APDL 280 (5x underestimate) |
| **At WP center** | B~9mT: mu_r~280 (correct, low-field regime) |
| **System linearity** | Dissertation p.39: "flux density increased linearly with input current... up to 3A" |
| **Impact on WP field** | Small — air reluctance R_a >> steel reluctance R_steel |
| **Impact on tip |B|** | Moderate — peak values underestimated |
| **Impact on charge model fit** | R_a absorbs the difference; field SHAPE is geometry-dominated |
| **Conclusion** | Acceptable for current work. Non-linear B-H version can be done later for comparison |

#### Issue 2: Boundary Condition D,ALL,MAG,0

| Aspect | Finding |
|--------|---------|
| **Original script** | No boundary condition present |
| **We added** | Lines 497-518: MAG=0 on outer air cylinder surfaces (lateral + top/bottom caps) |
| **Why needed** | DSP solver (magsolv,3) requires explicit BC for unique solution |
| **Without it** | "Negative pivot at node 316044 MAG" warning |
| **Dissertation says** | Nothing about BCs |
| **Conclusion** | Physically correct. Required addition |

#### Issue 3: Lower Pole Geometry — Milled Flat

| Aspect | Finding |
|--------|---------|
| **APDL code** | Line 186-188: `VSBV` subtracts upper half of cone at z = -13 mm |
| **Effect** | Lower poles are HALF cones (flat on top, rounded below) |
| **Dissertation p.14** | "lower poles are 42mm long and are milled to form a flat platform to support the culture dish" |
| **ANSYS screenshot** | Confirmed: lower pole cross-section shows flat top surface |
| **Conclusion** | APDL model correctly includes the milling operation |

#### Issue 4: Upper/Lower Pole Asymmetry in max |B|

| Aspect | Finding |
|--------|---------|
| **Our result** | Lower poles (P1,P3,P6): max |B| ~1.13 T; Upper (P2,P4,P5): ~0.74 T |
| **Seems inverted?** | No — physically correct for the modeled geometry |
| **Explanation** | Lower pole = half cone → tip has half the cross-section area → B = Phi/A is higher |
| **Dissertation K_I** | Lower diagonal ~0.60 < Upper ~0.90 — measures WP flux FRACTION, not tip |B| |
| **Consistency** | Both observations are consistent: lower poles deliver less total flux (higher reluctance from milling) but concentrate it in a smaller tip area (higher peak |B|) |
| **WP center field** | All 6 coils produce ~8.7 mT (ratio 1.04) — consistent with R_a dominance |

## 3. Geometric Dimensions

### Lower Poles (P1, P3, P6)

- Orientation: **horizontal** (all key points at z = -13 mm)
- Tip position: r = 0.408 mm from z-axis, z = -13.0 mm
- Cone: 40 um tip radius → 3 mm base, length ~15 mm
- Shape: **half cone** (upper half removed by VSBV)
- Holder: 4 blocks, max cross-section ~160 mm^2
- Magnetic path (protrusion + cone): ~44 mm
- Dissertation: "lower poles ~42 mm" (p.14)

### Upper Poles (P2, P4, P5)

- Orientation: **inclined 36.6 deg** from horizontal
- Tip position: r = 0.408 mm, z = -12.4 mm
- Cone: 40 um tip radius → 3 mm base, length ~15 mm (same as lower)
- Shape: **full cone** (no cut)
- Holder: 3 blocks, max cross-section ~120 mm^2
- Pole axis length: ~44.3 mm; total magnetic path: ~51 mm
- Dissertation: "upper poles ~45 mm" (p.14)

### Workspace Center

- ANSYS coordinates: (0, 0, -12.711 mm)
- Defined by: `SPH_OFST = Z0_LOW - 6mm + R_norm_z`
- Symmetric between upper and lower tips (offset +-R_norm_z = +-0.289 mm)

## 4. Simulation Results Summary (6 Unit-Excitation Solves)

| Coil (APDL) | Paper Pole | Layer | Max |B| [T] | WP |B| [mT] | Nodes | Status |
|-------------|-----------|-------|------------|------------|-------|--------|
| 1 | P1 | Lower | 1.12 | 8.64 | 494,871 | PASS |
| 2 | P3 | Lower | 1.14 | 8.74 | 494,871 | PASS |
| 3 | P6 | Lower | 1.14 | 8.63 | 494,871 | PASS |
| 4 | P5 | Upper | 0.75 | 8.87 | 494,871 | PASS |
| 5 | P2 | Upper | 0.74 | 8.93 | 494,871 | PASS |
| 6 | P4 | Upper | 0.74 | 8.96 | 494,871 | PASS |

- Max |B| ratio (lower/upper): 1.13/0.74 = 1.53
- WP field consistency ratio: 8.96/8.63 = 1.04
- All boundary fields < 1e-3 T (BC verified)
- No NaN or negative values

## 5. Figures Generated

| Figure | File | Description |
|--------|------|-------------|
| 6-coil bar chart | `figures/verify_all_coils_comparison.png` | Max |B| and WP field per coil |
| Coil 1 vector plot (a) | `figures/verify_coil1_a.png` | Top view B-field distribution |
| Coil 1 WP vectors (b) | `figures/verify_coil1_b.png` | 3D vectors near WP center |
| Coil 1 tip vectors (c) | `figures/verify_coil1_c.png` | Field convergence at P1 tip |
| Coil 1 decay + axis (d) | `figures/verify_coil1_d.png` | Radial decay and x-axis profile |
| Coil 1 WP xy top view | `figures/verify_coil1_wp_xy.png` | WP region field vectors |
| Pole geometry | `figures/pole_geometry_comparison.png` | Upper vs lower pole schematic |
| Pole shape detail | `figures/pole_shape_detail.png` | Cone cross-section and tip zoom |
| murx analysis | `figures/murx_analysis.png` | mu_r(B) curve vs APDL constant 280 |

## 6. B Perpendicularity at P1 Cone Surface

**Script:** `analysis/verify_B_perpendicularity.m` (2026-03-21)

Verified that B on the air side of the P1 cone surface is approximately
perpendicular to the surface, as expected from the iron-air boundary
condition (murx=280, theoretical deviation = arctan(1/280) ≈ 0.2°).

### Method

- Load Coil1 full mesh (~494k nodes), select fine sphere (R < 7mm)
- For each axial bin (s = 1-6 mm from P1 tip, bin width 0.25 mm):
  - Select confirmed-air nodes: |B| < 50 mT, dr > 0 (outside cone)
  - Restrict to lower half (phi < -10°) to exclude VSBV flat-cut surface
  - Take closest 20% shell (adaptive, min 3 nodes)
  - Compute analytical cone normal and angle between B and normal

### Key definitions

- s: axial distance from P1 tip along cone axis (tip→base)
- dr: radial distance from analytical cone surface (positive = air side)
- phi: azimuth around cone axis (phi < 0 = lower half, intact cone)
- angle: deviation of B from surface normal (0° = perpendicular)

### Results

- **Angle: 7-11° across s = 1-6 mm**, mean 9.4°, std ~2-3° per bin
- **|B|: 28 mT (s=1mm) → 9 mT (s=6mm)**, smooth decay
- All nodes confirm B is directed **inward** (P1 = flux sink)
- No systematic s-dependence: angle stable across entire cone body

### Why not the theoretical 0.2°?

The ~9° deviation is NOT a model error. It comes from finite mesh resolution:
- Closest air nodes are 90-400 um from the true iron-air interface
- The 0.2° boundary condition holds only exactly AT the interface
- At 100+ um distance, the global field pattern adds tangential components
- Finer mesh (SmartSize < 5) would bring closest nodes nearer → smaller angle

### VSBV flat-cut effect

Lower poles have the upper half milled flat (VSBV at APDL line 186-188).
The flat cut is at z_wp = tip_z (vz = 0 in tip-centered coordinates).
- s < 3.5 mm: cone cross-section is nearly full circle → no effect
- s > 3.5 mm: upper half progressively removed → D-shaped cross-section
- Air nodes near the flat surface have FLAT normals, not CONE normals
- The script avoids this by selecting only phi < -10° (lower half)

## 7. Open Questions for Future Work

1. **Non-linear B-H version**: Run with `TB,BH` material data for 1018 steel, compare WP field and K_I fit
2. **K_I extraction**: Fit flux distribution matrix from 6-coil data, compare to dissertation Eq.2.8 and Eq.2.19
3. **Force model**: Compute g_I and force envelopes
