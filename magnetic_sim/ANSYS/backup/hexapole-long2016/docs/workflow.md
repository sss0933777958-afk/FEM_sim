# Workflow: Simulation to Publication

## Overview
```
Stage 1          Stage 2           Stage 3              Stage 4
APDL scripts --> ANSYS batch   --> MATLAB import     --> MATLAB plots
(6 scripts)      (6 result sets)   (superposition)      (figures)
```

## Stage 1: APDL Simulation
- Run each CoilN.txt independently (only that coil's current = 1A)
- Each produces a full FEA solution in results/coilN/
- Scripts are identical except for CURR_ARRAY values
- **Input:** MT_Modeling_Geometry_Meshing_Solving_CoilN.txt
- **Output:** results/coilN/*.rst, *.db, *.full, etc.

## Stage 2: Data Extraction (POST1)
- Enter /POST1 after solve completes
- Extract node coordinates: NLIST,ALL,,,COORD
- Extract B field components: PRNSOL,B
- Redirect output to text files using /OUTPUT
- **Input:** ANSYS result files
- **Output:** Text files with node coordinates and B-field values

## Stage 3: MATLAB Linear Superposition
- Read 6 text files (one per coil unit response)
- For arbitrary current vector [I1, I2, I3, I4, I5, I6]:
  - B_total = I1*B_coil1 + I2*B_coil2 + ... + I6*B_coil6
- This works because the simulation is linear (constant permeability)
- **Input:** 6 exported text files
- **Output:** data/*.mat, data/*.csv

## Stage 4: Publication Figures
- MATLAB plotting scripts in analysis/
- Generate vector field plots, contour maps, force calculations
- Export to figures/ as .png and .eps
- **Input:** data/*.mat
- **Output:** figures/*.png, figures/*.eps

## Data Flow Diagram
```
CoilN.txt --> results/coilN/ --> (POST1 export) --> data/coilN_nodes.txt
                                                    data/coilN_Bfield.txt
                                                         |
                                              MATLAB superposition
                                                         |
                                                    data/B_combined.mat
                                                         |
                                                  figures/*.png, *.eps
```
