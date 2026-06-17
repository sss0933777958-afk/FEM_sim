# Notation Glossary

Unified notation mapping between Fei Long's 2016 dissertation, our APDL scripts, and MATLAB analysis code. **This is the single source of truth for all symbol/term usage in this project.**

## Physical Quantities

| Dissertation Symbol | Name (EN) | Name (ZH) | APDL Variable | MATLAB Variable | Units | Dissertation Ref |
|---|---|---|---|---|---|---|
| **B** | Magnetic flux density | 磁通密度 | (solver output) | `bx, by, bz` | T (Tesla) | Eq. 2.1-2.3 |
| \|**B**\| | Flux density magnitude | 磁通密度大小 | BSUM | `bsum` | T | Fig. 2.3-2.4 |
| **Phi** (Phi_i) | Magnetic flux through pole i | 第 i 極磁通量 | (implicit) | (not yet computed) | Wb (Weber) | Eq. 2.1 |
| **q_i** | Magnetic charge at pole i | 第 i 極等效磁荷 | (implicit) | (not yet computed) | A*m | Eq. 2.1: q_i = Phi_i / mu_0 |
| **Q** | Charge vector [q_1...q_6]^T | 磁荷向量 | (implicit) | (not yet computed) | A*m | Eq. 2.3 |
| **I** | Current vector [I_1...I_6]^T | 電流向量 | `CURR_ARRAY` | (not yet computed) | A (Ampere) | Eq. 2.4 |
| **I_hat** | Normalized current I/I_max | 歸一化電流 | — | — | dimensionless | Eq. 2.10 |
| **F** | Magnetic force on bead | 磁珠受力 | — | (not yet computed) | N or pN | Eq. 2.5-2.7 |
| **F_hat** | Normalized force F/F_N | 歸一化力 | — | — | dimensionless | Eq. 2.10 |
| **m** | Bead magnetic moment | 磁珠磁矩 | — | — | A*m^2 | Eq. 4.12 |
| **V_H** | Hall sensor voltage vector | 霍爾電壓向量 | — | — | V | Eq. 4.1-4.4 |

## Key Matrices

| Dissertation Symbol | Name (EN) | Name (ZH) | Dimension | Definition | Dissertation Ref |
|---|---|---|---|---|---|
| **K_I** | Flux distribution matrix | 磁通分佈矩陣 | 6x6 | Phi = K_I * F_mmf / R_a (Eq. 2.4) | Eq. 2.8 (nominal), Eq. 2.19 (calibrated) |
| **R_hat** | Charge-bead distribution matrix | 磁荷-磁珠分佈矩陣 | 3x6 | B = (k_m/rho^2) * R_hat * Q | Eq. 2.3, p.18 |
| **L_i** | Charge-bead gradient matrix (axis i) | 磁荷-磁珠梯度矩陣 | 6x6 | F_i = f_Phi * Phi^T * L_i * Phi / mu_0^2 | Eq. 2.6, p.22 |
| **D_H** | Hall flux-gain matrix | 霍爾磁通增益矩陣 | 6x6 diag | Phi = D_H * V_H | Eq. 4.2, p.67 |

## Model Parameters

| Dissertation Symbol | Name (EN) | Name (ZH) | Value | APDL/MATLAB Variable | Dissertation Ref |
|---|---|---|---|---|---|
| **rho** (physical) | Physical workspace radius | 物理工作距離 | 500 um | `R_norm` | p.14 |
| **rho** (fitted) | Effective charge location | 等效磁荷位置 | 900 um | (post-fit parameter) | p.20 |
| **R_a** | Lumped air reluctance | 集總空氣磁阻 | 6.3e8 A/Wb | (post-fit parameter) | p.20 |
| **N_c** | Coil turns per pole | 每極線圈匝數 | 70 | `TURNS` | p.14 |
| **I_max** | Maximum operating current | 最大操作電流 | 3 A | — | p.39 |
| **k_m** | mu_0 / (4*pi) | 磁常數 | 1.0e-7 N/A^2 | — | p.18 |
| **mu_0** | Vacuum permeability | 真空磁導率 | 4*pi*1e-7 H/m | — | — |
| **mu_r** | Relative permeability (1018 steel) | 相對磁導率 | 280 (APDL linear) | `murx` (MAT_MT) | p.14 |
| **g_I** | Current-based force gain | 電流力增益 | [7.56, 8.55, 7.62] pN (calibrated option 3) | (not yet computed) | Eq. 2.7, p.38 |
| **f_Phi** | Flux-based force gain | 磁通力增益 | — | — | Eq. 2.6, p.22 |
| **f_hat** | Hall-sensor force gain | 霍爾力增益 | 24.939 pN | — | p.105 |
| **m_s** | Bead magnetic saturation | 磁珠飽和磁矩 | 1.40e-12 A*m^2 (recal.) | — | p.105 |
| **a** | Langevin shape parameter | 朗之萬形狀參數 | 244.2 (recal.) | — | p.105 |
| **gamma** | Viscous drag coefficient | 黏滯阻力係數 | 8.5e-6 N*s/m (glycerol) | — | p.34 |

## Geometric Parameters

| Dissertation Name | APDL Variable | Value | Description (ZH) | Dissertation Ref |
|---|---|---|---|---|
| Pole tip radius | `POLE_TIP_R` | 40 um | 極尖半徑 | p.14 |
| Pole diameter | `POLE_R` * 2 | 6 mm | 極桿直徑 | p.14 |
| Cone length | `POLE_CONE_LEN` | 15 mm | 錐體長度 | (not stated) |
| Yoke inner radius | `YOKE_IN_R` | 42 mm | 軛環內徑 | (not stated) |
| Yoke outer radius | `YOKE_OUT_R` | 53 mm | 軛環外徑 | (not stated) |
| Yoke thickness | `YOKE_H` | 2 mm | 軛環厚度 | (not stated) |
| Protrusion radius | `PROT_R` | 5 mm | 突起半徑 | (not stated) |
| Protrusion height | `PROT_H` | 7 mm | 突起高度 | (not stated) |
| Coil inner radius | `COIL_IN_R` | 5 mm | 線圈內徑 | p.14 |
| Coil outer radius | `COIL_OUT_R` | 8 mm | 線圈外徑 | p.14 |
| Coil height | `COIL_H` = `PROT_H` | 7 mm | 線圈高度 | p.14 |
| Fine mesh sphere radius | `SPH_FINE_R` | 7 mm | 細網格球半徑 | (ANSYS only) |
| Outer air cylinder radius | `AIR_CYL_R` | 80 mm | 外部空氣圓柱半徑 | (ANSYS only) |
| Outer air cylinder height | `AIR_CYL_H` | 70 mm | 外部空氣圓柱高度 | (ANSYS only) |

## Coordinate Systems

| Dissertation Name | Symbol | Description | Dissertation Ref |
|---|---|---|---|
| Measurement coordinate | {O; x_m, y_m, z_m} | Aligned with microscope/camera; our ANSYS global CS | p.17, Fig. 2.5(a) |
| Actuation coordinate | {O; x_a, y_a, z_a} | Aligned with pole-pair axes; rotated from measurement | p.17, Fig. 2.5(b) |
| Rotation matrix | R_ma | Measurement-to-actuation transform | p.17 |
| WP center (ANSYS coords) | (0, 0, SPH_OFST) | SPH_OFST = -PROT_H - 6mm + R_norm_z ~ -12.711 mm | APDL line 62, 316 |

**Note**: ANSYS global coordinate system corresponds to the **measurement coordinate** in the dissertation. The actuation coordinate is a rotated frame used only in the force/control model, not in FEM.

## Pole Naming Convention

### Paper Poles (Dissertation Fig. 2.1, p.12)

| Paper Name | Angle (deg) | Layer | Opposite Pole | Actuation Axis |
|---|---|---|---|---|
| **P1** | 0 | Lower | P2 | +x_a |
| **P2** | 180 | Upper | P1 | +x_a |
| **P3** | 120 | Lower | P4 | +y_a |
| **P4** | 300 | Upper | P3 | +y_a |
| **P5** | 60 | Upper | P6 | +z_a |
| **P6** | 240 | Lower | P5 | +z_a |

### APDL-to-Paper Mapping

| APDL Coil Index | APDL Angle (deg) | Layer | Paper Pole |
|---|---|---|---|
| 1 | 0 | Lower | **P1** |
| 2 | 120 | Lower | **P3** |
| 3 | 240 | Lower | **P6** |
| 4 | 60 | Upper | **P5** |
| 5 | 180 | Upper | **P2** |
| 6 | 300 | Upper | **P4** |

### In Code

```matlab
% MATLAB mapping arrays (from mt_constants.m)
apdl_to_paper = {'P1','P3','P6','P5','P2','P4'};  % apdl_to_paper{i} = paper name of APDL coil i
pole_labels   = {'P1','P2','P3','P4','P5','P6'};   % paper order
pole_angles   = [0, 180, 120, 300, 60, 240];        % degrees, paper pole order P1-P6
```

## Nominal K_I Matrix (Dissertation Eq. 2.8)

Row/column order: [P1, P2, P3, P4, P5, P6] (paper order)

```
K_I = [ 5/6  -1/6  -1/6  -1/6  -1/6  -1/6 ]
      [-1/6   5/6  -1/6  -1/6  -1/6  -1/6 ]
      [-1/6  -1/6   5/6  -1/6  -1/6  -1/6 ]
      [-1/6  -1/6  -1/6   5/6  -1/6  -1/6 ]
      [-1/6  -1/6  -1/6  -1/6   5/6  -1/6 ]
      [-1/6  -1/6  -1/6  -1/6  -1/6   5/6 ]
```

Physical meaning: when coil j is excited with MMF = N_c * I_j, pole j gets 5/6 of the total flux, each of the other 5 poles gets -1/6 (return path).

## Units Convention

| Quantity | ANSYS Output | MATLAB Analysis | Dissertation Figures | Conversion |
|---|---|---|---|---|
| Coordinates | m (MKS) | um (for WP plots), mm (for device plots) | um / mm | 1 m = 1e6 um |
| B-field | T (Tesla) | mT (for WP plots) | Gauss (Fig. 2.4) | 1 T = 1e3 mT = 1e4 Gauss |
| Force | — | pN | pN | — |
| Current | A | A | A | — |
| Flux | — | Wb | Wb | — |

**Important**: Dissertation Fig. 2.4 contour plots use **Gauss**. When comparing our mT plots: **1 mT = 10 Gauss**.

## Model / Fitting Variables

Mapping between dissertation notation, MATLAB variable names, and their roles in the fitting pipeline.

| Concept | Paper Symbol | MATLAB Variable | Dimension | Notes |
|---|---|---|---|---|
| Charge position | c_i | `pos(:,i)` | 3x6 | c_i = ell * d_hat_i + delta_i |
| Baseline sphere radius | ell (script L) | `ell` or `ell_fixed` | scalar | Shared across all methods |
| Bias/offset vector | b_i | `delta(:,i)` | 3x6 | Code uses `delta` to avoid confusion with `b_fem` |
| Direction unit vector | d_hat_i | `d_hat(:,i)` | 3x6 | WP center -> pole tip direction, fixed by geometry |
| Air reluctance | R_a | `R_a` or `R_a_k` | scalar or 6x1 | Per-coil when using C_k; shared when using C |
| Lumped amplitude | C = N_c/(mu_0*R_a) | `C` or `C_k` | scalar or 6x1 | VarPro solves analytically |
| Magnetic charge vector | Q | `-(N_c/(mu_0*R_a))*K_I*I_vec` | 6x1 | Implicit in code |
| Flux distribution matrix | K_I | `K_I` | 6x6 | = I_6 - ones(6)/6 (nominal, Eq. 2.8) |
| Current vector (model) | I_diss | `I_vec` | 6x1 | Includes coil_sign correction |
| Weight vector | w | `w` = K_I * I_vec | 6x1 | Determines per-pole contribution |
| Coil sign correction | — | `coil_sign` | scalar | +1 (lower), -1 (upper); see coil-winding-sign-convention.md |

## Fitting Methods

| Tag | Description | Data | Nonlinear DOF | Analytic DOF | Script |
|---|---|---|---|---|---|
| [A] | Shared ell, b=0 | 1 coil | 1 (ell) | 1 (C -> R_a) | `fit_charge_model.m` |
| [a] | Per-pole ell_i, b=0 | 1 coil | 6 (ell_i) | 1 (C -> R_a) | (inline in test script) |
| [J] | Free 3D positions | 6 coils joint | 18 (3x6 coords) | 6 (C_k) | (test script) |
| [D] | ell(fixed) + delta(3D) | 6 coils joint | 18 (delta) | 6 (C_k) | `fit_ell_delta_6coil.m` |
| [B-sc] | ell(fixed) + delta(3D), shared C | 1 coil | 18 (delta) | 1 (C) | `fit_single_coil_with_bias.m` |
| [B-6x] | ell(fixed) + delta(3D), shared C | all6 | 18 (delta) | 1 (C) | `fit_all6_with_bias.m` |

## ANSYS-Specific Terms (Not in Dissertation)

| Term | Description (ZH) | Why Needed |
|---|---|---|
| `SOLID96` | 3D 磁標量位元素 | Volume meshing for magnetostatic |
| `SOURC36` | 電流源元素 | Coil primitive definition |
| `magsolv,3` | DSP (微分標量位) 求解器 | Magnetostatic solver method |
| `D,ALL,MAG,0` | 遠場邊界條件 (MAG=0) | Required by DSP for unique solution |
| `SmartSize 5` | 自適應網格等級 5 | Mesh refinement level |
| `BSUM` | \|B\| = sqrt(BX^2+BY^2+BZ^2) | ANSYS post-processing quantity name |
| `PRNSOL,B,COMP` | 列印 B 分量 (BX,BY,BZ,BSUM) | POST1 data extraction command |
