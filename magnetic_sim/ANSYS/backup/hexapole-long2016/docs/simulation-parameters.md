# Simulation Parameters

## Geometric Parameters
| Parameter | Value | Description |
|-----------|-------|-------------|
| R_norm | 500 um | Working radius (center of pole tips) |
| YOKE_IN_R | 42 mm | Inner radius of yoke ring |
| YOKE_OUT_R | 53 mm | Outer radius of yoke ring |
| YOKE_H | 2 mm | Yoke thickness |
| PROT_R | 5 mm | Protrusion radius |
| PROT_H | 7 mm | Protrusion height |
| POLE_R | 3 mm | Pole base radius |
| POLE_TIP_R | 40 um | Pole tip radius |
| POLE_CONE_LEN | 15 mm | Pole cone length |
| SPH_FINE_R | 7 mm | Fine mesh sphere radius (around tips) |
| AIR_CYL_R | 80 mm | Outer air cylinder radius |
| AIR_CYL_H | 70 mm | Outer air cylinder height |

## Pole Arrangement
- 6 poles total: 3 lower + 3 upper
- Lower poles (1-3): 120 degrees apart, below yoke
- Upper poles (4-6): 120 degrees apart, offset 60 degrees, above yoke
- Poles converge toward R_norm center point

## Material Properties
| MAT ID | Material | murx | Usage |
|--------|----------|------|-------|
| 1 | Air | 1 | Center sphere region |
| 2 | 1018 Steel | 280 | Yoke + poles |
| 3 | Air | 1 | Outer cylinder region |

## Element Types
| ET ID | Type | Usage |
|-------|------|-------|
| 1 | SOLID96 | 3D magnetic scalar potential (volumes) |
| 2 | SOURC36 | Current source (coil primitives) |

## Coil Specifications
- **Turns:** 70
- **Inner radius:** 5 mm
- **Outer radius:** 8 mm
- **Height:** 7 mm (= PROT_H)
- **Current:** 1 A per active coil (unit excitation for superposition)

## Mesh Settings
- **SmartSize level:** 5 (SMRT,5)
- **Shape:** Free tetrahedral (mshape,1,3D / mshkey,0)
- **Total elements:** ~2.4 million
- **Fine sphere:** separate refinement zone around pole tips

## Solver
- **Method:** magsolv,3 (Differential Scalar Potential, DSP)
- **Boundary condition:** MAG=0 on all outer air cylinder surfaces
- **Post-processing:** /POST1 with NLIST and PRNSOL,B
