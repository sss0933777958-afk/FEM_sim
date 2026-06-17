# FEM Simulation Workspace

`FEM_sim/` 是通用 FEM 模擬容器。磁學模擬（磁鑷／hexapole／quadrupole）全部收在 `magnetic_sim/` 類別下；未來其他 FEM 類別會與之並列。

ANSYS APDL magnetostatic finite element simulation of hexapole magnetic tweezers for biophysics research.

## Overview

This project implements FEM simulation of a 6-pole magnetic tweezers device, computing the magnetic field (B) under unit-excitation for each coil. Results feed into a point-charge model fitting pipeline for force calibration.

Based on: **Fei Long, "Design, Fabrication, and Calibration of a Hexapole Magnetic Tweezers," PhD Dissertation, Ohio State University, 2016.**

## Repository Structure

```
FEM_sim/                     通用 FEM 模擬容器
├── magnetic_sim/            磁學模擬類別
│   └── ANSYS/               ANSYS 求解器子層（未來可並列 COMSOL/ 等）
│       ├── main/            ★ 活躍設計：4-pole MEMS Quadrupole (Harrison-style；原 kuo/)
│       │   ├── apdl/        ANSYS APDL simulation & extraction scripts
│       │   ├── matlab/      MATLAB analysis (含 common/ resolver)
│       │   ├── ANSYS_data/  FEM .dat/.db (gitignored)
│       │   ├── MATLAB_data/ MATLAB outputs (.mat/.csv)
│       │   ├── figures/     Figures (.png)
│       │   └── ...          CAD/ IGES/ comsol/ mph/ doc/ ...
│       └── backup/          歸檔（非活躍設計）
│           ├── hexapole-long2016/  Long 2016 dissertation hexapole design
│           └── hung/        Hung hexapole design
│
└── (future: electric_sim/ 等其他 FEM 類別，與 magnetic_sim/ 並列)
```

## Prerequisites

- **ANSYS MAPDL** 2025 R2 (or compatible version)
- **MATLAB** R2024b+ (Optimization Toolbox for fitting)
- ~60 GB disk space for full simulation results (6 coils)

## Quick Start

Run a single coil simulation (batch mode):
```bash
cd magnetic_sim/hexapole-long2016
"C:\Program Files\ANSYS2025R2\v252\ansys\bin\winx64\MAPDL.exe" -b -np 4 -m 24000 \
  -dir "results/coil1" -j "coil1" \
  -i "$(pwd)/apdl/MT_Modeling_Geometry_Meshing_Solving_Coil1.txt" \
  -o "results/coil1/solve.out"
```

Process results in MATLAB:
```matlab
cd magnetic_sim/hexapole-long2016/analysis
fit_charge_model        % [A] baseline fit
fit_all6_with_bias      % [B-6x] final 19-parameter fit
```

## Key Results (hexapole-long2016)

| Method | Parameters | Error | R_a (A/Wb) |
|--------|-----------|-------|------------|
| [A] Baseline | ell=835 um | 4.94% | 9.21e8 |
| [J] Joint 6-coil | ell=766-818 um | 1.11% | ~1.01e9 |
| [B-6x] Final | 19 params | 0.07% | 1.03e9 |

## References

- Long, F. (2016). PhD Dissertation, Ohio State University.
- Zhang, Z. & Menq, C.H. (2011). IEEE/ASME Trans. Mechatronics.
- Long, F., Matsuura, T. & Bhatt, D. (2016). Actively Controlled Hexapole Electromagnetic Actuating System.
