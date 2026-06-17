"""Build 6 gap200um_mu_eq sim variants — μ_r equivalent for 200 µm gap.

Strategy:
  - Keep Long2016 source geometry COMPLETELY intact (no geom changes, no
    V_ID issues, no mesh failures).
  - Add a new material MAT_PROT with μ_r = 31, which gives equivalent
    reluctance for 7 mm post + 200 µm air gap in series with steel
    μ_r = 280.
  - After meshing, reassign elements in the 6 protrusion regions from
    MAT_MT (μ_r=280) to MAT_PROT (μ_r=31).
  - Restrict to elements with current MAT = MAT_MT so air elements in
    the protrusion bbox aren't accidentally turned into steel.

Equivalent calculation:
    R_steel_post = h/(μ_r·μ₀·A) = 7e-3 / (280·4πe-7·π·25e-6) ≈ 2.54e5 A/Wb
    R_gap_200µm  = l_gap/(μ₀·A) = 200e-6/(4πe-7·π·25e-6)    ≈ 2.03e6 A/Wb
    R_total      = R_steel + R_gap ≈ 2.28e6 A/Wb
    μ_r_eff      = h/(R_total·μ₀·A) ≈ 31

Patches applied:
  1. Add MAT_PROT = 4 + mp,murx,MAT_PROT,31 (next to existing MP commands)
  2. After vmesh,all — reassign protrusion elements from MAT_MT to MAT_PROT
  3. /CWD → coilN_gap200um_mueq
  4. Current stays at 1.0 A (matches Long2016 verbatim baseline)

Usage:
    python _build_gap200um_mu_eq.py
"""
from pathlib import Path

SRC_DIR = Path(r"G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\apdl\long2016_hexapole_halfcut\sim")

MU_R_EFF = 31    # equivalent μ_r for 7 mm post + 200 µm gap (μ_r_steel = 280)

# ---------------------------------------------------------------------
# Patch 1: Add MAT_PROT after MAT_AIR2 definition + corresponding mp
# ---------------------------------------------------------------------
MAT_FROM = """MAT_AIR1 = 1   !----MAT number 1: Center air
MAT_MT   = 2   !----MAT number 2: magnetic tweezers
MAT_AIR2 = 3   !----MAT number 3: outer air
mp , murx , MAT_AIR1 , 1     !----MAT NUM 1 : center air
mp , murx , MAT_MT   , 280   !----MAT NUM 2 : 1018 steel
mp , murx , MAT_AIR2 , 1     !----MAT NUM 3 : outer air"""

MAT_TO = f"""MAT_AIR1 = 1   !----MAT number 1: Center air
MAT_MT   = 2   !----MAT number 2: magnetic tweezers
MAT_AIR2 = 3   !----MAT number 3: outer air
MAT_PROT = 4   !----MAT number 4: protrusion (μ_r equivalent for 200µm gap)
mp , murx , MAT_AIR1 , 1     !----MAT NUM 1 : center air
mp , murx , MAT_MT   , 280   !----MAT NUM 2 : 1018 steel (full)
mp , murx , MAT_AIR2 , 1     !----MAT NUM 3 : outer air
mp , murx , MAT_PROT , {MU_R_EFF}    !----MAT NUM 4 : prot equivalent (7mm steel + 200µm gap series)"""


# ---------------------------------------------------------------------
# Patch 2: After vmesh,all — reassign protrusion elements to MAT_PROT
# Anchor on the unique " vmesh, all" line. Insert reassignment block AFTER.
# ---------------------------------------------------------------------
VMESH_FROM = "vmesh,all"

VMESH_TO = """vmesh,all

! [gap200um μ_r-equivalent] Reassign 6 protrusion regions from MAT_MT to MAT_PROT
! to mimic 7 mm steel post + 200 µm yoke-pole air gap in series.
! ESEL by bbox at each protrusion location, restrict to MAT_MT elements to avoid
! accidentally modifying air elements in the bbox region.
*do, prot_ind, 1, POLE_NUM, 1
    *if, prot_ind, LE, 3, THEN
        ang_prot = 120*(prot_ind - 1)
        zlo = -PROT_H
        zhi = 0
    *ELSE
        ang_prot = 120*(prot_ind - 4) + 60
        zlo = YOKE_H
        zhi = YOKE_H + PROT_H
    *ENDIF
    xc_prot = YOKE_MID_R*cos(ang_prot)
    yc_prot = YOKE_MID_R*sin(ang_prot)
    tol_xy  = PROT_R*1.1

    ESEL, S, CENT, X, xc_prot-tol_xy, xc_prot+tol_xy
    ESEL, R, CENT, Y, yc_prot-tol_xy, yc_prot+tol_xy
    ESEL, R, CENT, Z, zlo, zhi
    ESEL, R, MAT,  , MAT_MT
    EMODIF, ALL, MAT, MAT_PROT
    ALLSEL, ALL
*ENDDO

! Verify: count elements reassigned
ESEL, S, MAT, , MAT_PROT
*get, n_prot, ELEM, , COUNT
ALLSEL, ALL
/COM, ====================================================================
/COM, μ_r-equivalent reassignment: %n_prot% elements set to MAT_PROT (μ_r=31)
/COM, ===================================================================="""


def patch_mu_eq(pole_n: int) -> None:
    src = SRC_DIR / f"MT_Sim_P{pole_n}.txt"
    dst = SRC_DIR / f"MT_Sim_P{pole_n}_gap200um_mueq.txt"

    text = src.read_text(encoding="utf-8")

    # Patch 1: material definitions
    if MAT_FROM not in text:
        raise RuntimeError(f"MAT_FROM marker not found in {src.name}")
    text = text.replace(MAT_FROM, MAT_TO, 1)

    # Patch 2: vmesh,all → vmesh,all + EMODIF block
    if VMESH_FROM not in text:
        raise RuntimeError(f"VMESH_FROM marker not found in {src.name}")
    text = text.replace(VMESH_FROM, VMESH_TO, 1)

    # Patch 3: /CWD coilN → coilN_gap200um_mueq
    cwd_marker = f"coil{pole_n}'"
    if cwd_marker not in text:
        raise RuntimeError(f"/CWD marker not found in {src.name}")
    text = text.replace(cwd_marker, f"coil{pole_n}_gap200um_mueq'", 1)

    dst.write_text(text, encoding="utf-8")
    print(f"Wrote {dst.name}")


if __name__ == "__main__":
    for n in range(1, 7):
        patch_mu_eq(n)
    print(f"\nDone: 6 gap200um_mueq sims built (μ_r_prot = {MU_R_EFF}).")
