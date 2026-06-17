# Coil Winding Direction & Sign Convention

How the APDL SOURC36 coil winding creates opposite flux polarity at upper vs lower pole tips, and how the point-charge model corrects for this.

## 1. SOURC36 Three-Node Definition

All 6 coils use the **same** node pattern (APDL lines 411-449):

```
node1: (center_x + COIL_R,  center_y,          z_mid)   <- +x direction
node2: (center_x,           center_y + COIL_R, z_mid)   <- +y direction
node3: (center_x,           center_y,          z_mid)   <- coil axis point
```

Right-hand rule: fingers curl from (node1-node3) toward (node2-node3), thumb = **+z**.

**Result: all 6 coils drive flux in +z (upward) inside the iron core.**

## 2. Geometry: Why Upper and Lower Tips Have Opposite Polarity

```
                z (APDL coords)
                ^
   z=9mm  ------| <- upper protrusion top -> cone -> down to tip at z~-12.4mm
                |    +===========+
   z=5.5mm     |    | UPPER COIL|  B_iron = +z (upward)
                |    +===========+
   z=2mm  ------| <- yoke top
   z=0mm  ------| <- yoke bottom
                |    +===========+
   z=-3.5mm    |    | LOWER COIL|  B_iron = +z (upward)
                |    +===========+
   z=-7mm ------| <- lower protrusion bottom -> cone -> down to tip at z~-13mm
                |
   z~-12.7mm ---| <- WP center (SPH_OFST)
```

### Lower poles (P1, P3, P6)
- Coil at z = -3.5 mm
- **Tip is BELOW the coil** (z ~ -13 mm), yoke is ABOVE (z = 0)
- +z flux direction = from tip upward through iron to yoke
- Flux enters tip from air -> **tip = SINK (flux converges in)**

### Upper poles (P2, P4, P5)
- Coil at z = 5.5 mm
- **Yoke is BELOW the coil** (z = 2 mm), protrusion top is ABOVE (z = 9 mm)
- +z flux direction = from yoke upward through iron to protrusion top -> cone -> tip
- Flux exits tip into air -> **tip = SOURCE (flux diverges out)**

### Summary Table

| Pole Layer | Coil z  | Tip relative to coil | +z flux path (iron) | Tip polarity |
|------------|---------|---------------------|---------------------|-------------|
| Lower      | -3.5 mm | Below (z ~ -13 mm)  | tip -> coil -> yoke | **SINK**    |
| Upper      | +5.5 mm | Beyond top (via cone)| yoke -> coil -> tip | **SOURCE**  |

## 3. Sign Correction in Point-Charge Model

The dissertation model (Eq. 2.2-2.4) assumes the excited pole is always a **sink**:

```
Q = -(N_c / (mu_0 * R_a)) * K_I * I_eff
```

Since APDL gives opposite tip polarity for upper poles, we apply `coil_sign`:

```matlab
coil_sign = +1   % lower poles (P1, P3, P6) - naturally sink
coil_sign = -1   % upper poles (P2, P4, P5) - flip to make sink
I_eff = coil_sign * I_nominal
```

Implementation in `fit_charge_model_6coil.m` (line 62):
```matlab
paper_idx = apdl_to_paper_idx(k);
coil_sign = 1 - 2 * (1 - c.pole_is_lower(paper_idx));  % +1 lower, -1 upper
I_vec(paper_idx) = coil_sign;
```

## 4. Common Pitfalls

### Forgetting coil_sign when fitting
- **Symptom**: Upper coil fits give wrong R_a sign or physically nonsensical parameters
- **Fix**: Always apply `coil_sign` to `I_vec` before calling `point_charge_model()`

### Single-coil fitting (Coil 1 = P1)
- P1 is a **lower pole**, so `coil_sign = +1` and no correction needed
- This is why `fit_charge_model.m` (single-coil) works without sign correction
- The issue only appears when extending to upper pole coils (Coil 4/P5, Coil 5/P2, Coil 6/P4)

### Recovering R_a from single-charge fit
In `fit_3d_vector_all_coils.m`, recovering R_a from the fitted Q requires the sign:
```matlab
R_a = -(k_m * N_c * 5 * coil_sign) / (6 * mu_0 * Q_opt)
```
Missing `coil_sign` here gives R_a < 0 for upper poles.

### Multi-coil joint fitting
When stacking all 6 coils for joint fitting (`helper_cost_joint6.m`), each coil's `KI_w` must already include the `coil_sign` correction. The sign is baked into `I_vec` before computing `KI_w = K_I * I_vec`.

## 5. Verification

The sign convention was verified empirically by checking the B-field direction at WP center for each coil excitation:
- Lower coils (1,2,3): B at center points TOWARD the excited pole tip (sink)
- Upper coils (4,5,6): B at center points AWAY from the excited pole tip (source)
- After applying `coil_sign`, all 6 coils show B pointing toward the excited pole (unified sink convention)

See `fit_charge_model_6coil.m` header comments (lines 7-11) for the original documentation.
