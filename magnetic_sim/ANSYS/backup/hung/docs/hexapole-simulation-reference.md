# Hexapole Magnetic Tweezers — Simulation Reference

> Cross-design reference for all hexapole magnetic tweezers ANSYS simulations.
> This document defines mandatory constraints and modeling procedures.
> Any new hexapole design MUST satisfy all constraints listed here.

---

## Part A: Mandatory Geometric Constraints

### Constraint 1: Orthogonal Pair Axes

The 6 poles form 3 opposing pairs. The line connecting each pair of tips defines an **actuation axis**. The 3 actuation axes MUST be mutually orthogonal (pairwise dot product = 0).

```
Pair 1: P1 <--> P2    (actuation axis x_a)
Pair 2: P3 <--> P4    (actuation axis y_a)
Pair 3: P5 <--> P6    (actuation axis z_a)

x_a . y_a = 0
y_a . z_a = 0
z_a . x_a = 0
```

This is the fundamental requirement for 3D force control in arbitrary directions.

### Constraint 2: Tips on a Common Sphere

All 6 pole tips MUST lie on a sphere of radius `R_norm` centered at the working point (WP):

```
|tip_i - WP_center| = R_norm    for all i = 1..6
```

`R_norm` is a free design parameter (e.g., 500 um in the Long 2016 design).

### Constraint 3: 60-degree Azimuthal Offset

The 3 Upper poles are offset by 60 degrees in the azimuthal direction relative to the 3 Lower poles:

```
Lower poles (P1, P3, P6):  theta = 0, 120, 240 degrees
Upper poles (P2, P4, P5):  theta = 180, 300, 60 degrees
                           (= Lower angles + 180, each shifted by the pair relationship)
```

Equivalently: within each layer (upper or lower), the 3 poles are spaced 120 degrees apart, and the two layers are rotated 60 degrees relative to each other.

### Constraint 4: Fixed Polar Angle alpha = 54.74 degrees

**Derived from Constraints 1 + 2 + 3.** This is NOT a free parameter.

```
alpha = arctan(sqrt(2)) = 54.7356... degrees
```

This uniquely determines the XY and Z components of each tip position:

```
R_norm_xy = R_norm * sqrt(2/3)     (XY-plane projection)
R_norm_z  = R_norm / sqrt(3)       (Z component)
```

**Derivation sketch:**

A Lower tip at azimuthal angle theta has position:
  (R_norm_xy * cos(theta), R_norm_xy * sin(theta), -R_norm_z)

Its opposing Upper tip (offset by 180 degrees in theta, with +z) has position:
  (R_norm_xy * cos(theta+180), R_norm_xy * sin(theta+180), +R_norm_z)

The pair axis vector is:
  (-2 * R_norm_xy * cos(theta), -2 * R_norm_xy * sin(theta), 2 * R_norm_z)

Setting the dot product of two pair axes to zero and solving:

  sin^2(alpha) * cos(delta_theta) + cos^2(alpha) = 0

where delta_theta is the angle between two pair axes projected onto XY.
For the hexapole geometry (120 degree separation with 60 degree offset):

  tan^2(alpha) = 2
  alpha = arctan(sqrt(2)) = 54.74 degrees

Verification:
  R_norm_xy^2 + R_norm_z^2 = R_norm^2 * (2/3 + 1/3) = R_norm^2  (on sphere)

---

## Part B: ANSYS Modeling Procedure (8 Steps)

Extracted from the hexapole-long2016 Coil 1 APDL script. The procedure is general for any hexapole design.

### Step 1: Environment Setup

```apdl
/CLEAR,NOSTART
/PREP7
*afun, DEG
BTOL, 1e-6
/UNITS, MKS
```

Define all geometric parameters. The tip positions are fully determined by R_norm and alpha = 54.74 degrees:

```apdl
R_norm    = <design value>
R_norm_xy = R_norm * sqrt(2.0/3.0)
R_norm_z  = R_norm / sqrt(3.0)
```

### Step 2: Yoke + Protrusion

Build the yoke (hollow disk) and 6 cylindrical protrusions:
- Lower protrusions (1-3): extend downward from yoke bottom, at 0/120/240 degrees
- Upper protrusions (4-6): extend upward from yoke top, at 60/180/300 degrees
- Boolean union (VADD) all 7 volumes into one solid

### Step 3: Lower Poles (Loop, 3 poles)

For each lower pole:
1. Define keypoints for the axial cross-section (tip fillet arc + cone + base)
2. Create area from lines (AL)
3. Rotate area 360 degrees to get solid of revolution (VROTAT)
4. Add holder blocks, boolean union (VADD)
5. Subtract excess volume (VSBV)

### Step 4: Upper Poles (Loop, 3 poles)

Similar to lower poles but:
- Poles originate from yoke top and angle downward toward WP center
- Incline angle computed from tip position to protrusion center
- Full cone geometry (not truncated like lower)

### Hung-Specific Differences (vs Long2016)

The general procedure above follows Long2016. The Hung design differs in these steps:

**Pole construction (Steps 3-4):**
- Hung uses **two-part VROTAT** (front 180° D-shape + back 360° full round), NOT VSBV
- D-shape flat cut is achieved by rotating only 180° for the front section (0~28mm)
- Upper poles: `-side_perp` offset → flat face up
- Lower poles: `+side_perp` offset → flat face down

**Block positioning:**
- Hung blocks are **centered on pole end** (`fc_endz ± BLK_T/2`), NOT entry-point based
- This method works for any TILT_UP angle (the original `BLK_T/sin(TILT_UP)` formula fails for TILT_UP < 28°)

**Coil placement:**
- Hung coils sit on guide posts / upper cores (above blocks), NOT on protrusions
- Coil position Z = `block_top + COIL_DZ/2` (touching block upper face)
- Clockwise winding (N1/N2 swapped) → flux toward block (-Z direction)

**Guide posts + Upper cores:**
- Hung has 3 lower guide posts (P1, P3, P6) + 3 upper cores (P2, P4, P5)
- Both are vertical cylinders (R=4mm) connecting blocks to yoke
- Long2016 uses "protrusions" extending directly from the yoke

### Step 5: Air Domains

1. **Fine mesh sphere** (radius = SPH_FINE_R): centered at WP center, provides mesh refinement near tips
2. **Outer air cylinder** (radius = AIR_CYL_R, height = AIR_CYL_H): encloses entire device

### Step 6: Boolean Overlap

```apdl
VOVLAP,ALL       ! overlap all volumes
NUMCMP,VOLU      ! compress volume numbering
VGLUE, <iron volumes>   ! glue iron parts for mesh continuity
```

### Step 7: Material + Element + Coil + Mesh

**Materials:**
- MAT 1: Air (murx = 1) — center sphere
- MAT 2: Steel (murx = design value, must be constant for linear superposition) — yoke + poles
- MAT 3: Air (murx = 1) — outer cylinder

**Elements:**
- ET 1: SOLID96 — 3D magnetic scalar potential (all volumes)
- ET 2: SOURC36 — current source (coil primitives only)

**Coil definition:** See Part C below.

**Mesh:**
```apdl
smrt, <level>
mshape, 1, 3D       ! tetrahedral
mshkey, 0            ! free mesh
vmesh, all
```

### Step 8: Boundary Condition + Solve

```apdl
! Select all outer air cylinder surfaces (lateral + top + bottom)
! Apply MAG=0 on all selected nodes
D, ALL, MAG, 0

/SOLU
magsolv, 3, , , , , 1    ! DSP solver
```

---

## Part C: SOURC36 Coil Definition

SOURC36 coil primitives do NOT require volume geometry or mesh. Each coil is defined by:

### Real Constants

```apdl
R, <set_id>, 1, NI, COIL_THK, COIL_H
```

| Field | Meaning |
|-------|---------|
| R1 = 1 | Coil primitive type flag |
| R2 = NI | Ampere-turns (TURNS * current) |
| R3 = COIL_THK | Radial thickness of coil cross-section |
| R4 = COIL_H | Axial height of coil cross-section |

### 3-Node Definition

Each coil requires exactly 3 nodes:

```
node3 = coil center point
node1 = point on coil plane, offset by COIL_R in one direction (defines radius + plane)
node2 = point on coil plane, offset by COIL_R in perpendicular direction
```

```apdl
N, node1, X0 + COIL_R, Y0,          Z0
N, node2, X0,          Y0 + COIL_R, Z0
N, node3, X0,          Y0,          Z0       ! center
TYPE, ET_curr
REAL, <set_id>
E, node1, node2, node3
```

ANSYS internally computes a uniform current density J = NI / (COIL_THK * COIL_H) distributed over the toroidal cross-section, and uses Biot-Savart to inject the source into the field equations.

### Coil Placement

- Lower coils: centered on protrusion, at Z = protrusion_bottom + COIL_H/2
- Upper coils: centered on protrusion, at Z = protrusion_top - COIL_H/2
- COIL_R = (COIL_IN_R + COIL_OUT_R) / 2 (average radius of coil cross-section)

---

## Part D: Excitation and Superposition Strategy

### Unit Excitation Principle

Because the material model is **linear** (constant murx), the magnetic field equation is linear. This enables superposition:

```
B_total(I1, I2, ..., I6) = I1 * B_coil1 + I2 * B_coil2 + ... + I6 * B_coil6
```

### Procedure

1. Run N independent simulations (N = number of coils)
2. Each simulation: set CURR_ARRAY so only one coil has current = 1A, rest = 0
3. The N solutions form a complete basis
4. Any arbitrary current combination is a linear combination of these basis solutions

### CURR_ARRAY for Each Simulation

For a 6-pole design:

| Simulation | CURR_ARRAY |
|-----------|------------|
| Coil 1 | [1, 0, 0, 0, 0, 0] |
| Coil 2 | [0, 1, 0, 0, 0, 0] |
| Coil 3 | [0, 0, 1, 0, 0, 0] |
| Coil 4 | [0, 0, 0, 1, 0, 0] |
| Coil 5 | [0, 0, 0, 0, 1, 0] |
| Coil 6 | [0, 0, 0, 0, 0, 1] |

All scripts share identical geometry, materials, mesh, and boundary conditions. Only CURR_ARRAY differs.

### Critical Requirement

**murx MUST be a constant** (not a B-H curve) for superposition to be valid. If nonlinear materials are used, each current combination must be solved independently, and the number of simulations grows from N to the full space of interest.

---

## Part E: New Design Checklist

Before creating APDL scripts for a new hexapole design, verify:

### Geometry
- [ ] R_norm defined (all tips equidistant from WP center)
- [ ] R_norm_xy = R_norm * sqrt(2/3)
- [ ] R_norm_z = R_norm / sqrt(3)
- [ ] Alpha = 54.74 degrees (NOT adjustable)
- [ ] Lower poles at 0, 120, 240 degrees
- [ ] Upper poles at 60, 180, 300 degrees
- [ ] 3 pair axes are orthogonal (verify dot products = 0)

### Materials
- [ ] murx is a constant (required for linear superposition)
- [ ] Air regions have murx = 1

### Elements
- [ ] SOLID96 for all volumes
- [ ] SOURC36 for coil primitives

### Coils
- [ ] Each coil defined by 3 nodes + Real Constants
- [ ] NI = TURNS * CURR_ARRAY(i)
- [ ] Coil center on protrusion axis at correct Z height

### Boundary & Solver
- [ ] MAG = 0 on all outer air domain surfaces (lateral + top + bottom)
- [ ] magsolv,3 (DSP method)
- [ ] Outer air domain large enough for field decay

### Excitation
- [ ] N scripts with identical geometry, only CURR_ARRAY differs
- [ ] Each script excites exactly one coil with 1A
