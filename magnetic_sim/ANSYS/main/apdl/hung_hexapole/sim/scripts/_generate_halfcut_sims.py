"""Generate 6 halfcut sim scripts from the full-hexapole templates.

Differences vs. full-hexapole simulations:
  1. /CWD points to magnetic_sim/ANSYS/main/ANSYS_data/long2016_hexapole_halfcut/coilN
  2. Each lower pole gets a VSBV cut after its VROTAT (removes upper half
     of the cone segments, leaving a D-shaped half-cone).

Run from any shell:
    python _generate_halfcut_sims.py
"""
from pathlib import Path
import re

REPO = Path(r"G:\my_workspace\code\FEM_sim")
SRC_DIR = REPO / "kuo" / "apdl" / "sim" / "long2016_hexapole_full"
DST_DIR = REPO / "kuo" / "apdl" / "sim" / "long2016_hexapole_halfcut"
DST_DIR.mkdir(parents=True, exist_ok=True)

# The FILLED-cone marker block in the source (will be replaced with halfcut VSBV)
FILLED_MARKER = """    AL, LI_START, LI_START + 1, LI_START + 2, LI_START + 3, LI_START + 4
    VROTAT, AR_START, ,,,,,KP_START,KP_START + 4,,4

    ! [NOTE] NO VSBV cut here — keep cone FILLED (bidirectional)
*enddo"""

# The halfcut replacement: VSBV the upper half off the cone segments.
HALFCUT_REPLACEMENT = """    AL, LI_START, LI_START + 1, LI_START + 2, LI_START + 3, LI_START + 4
    VROTAT, AR_START, ,,,,,KP_START,KP_START + 4,,4

    ! [ADDED] Half-cut lower pole: subtract upper-half block from cone segments
    ! [ADDED] Pre-boolean: capture cone segment numbers before VSBV scrambles them
    vl_cone_first = VL_START + 4
    *get, vl_cone_max, VOLU, , NUM, MAX
    vl_cut = vl_cone_max + 1
    numstr, VOLU, vl_cut
    BLOCK, 0, X0_LOW - 10e-3, Y0_LOW - 10e-3, Y0_LOW + 10e-3, Z0_LOW - 6e-3, Z0_LOW
    FLST, 2, (vl_cone_max - vl_cone_first + 1), 6, ORDE, 2
    FITEM, 2, vl_cone_first
    FITEM, 2, -vl_cone_max
    VSBV, P51X, vl_cut
*enddo"""

# Header banner — update the comment that says "no VSBV cut" to "half-cut"
HEADER_FROM = "Long Fei FULL HEXAPOLE"
HEADER_TO   = "Long Fei HALF-CUT HEXAPOLE (lower poles milled)"
GEOM_FROM   = "filled cones, no VSBV cut"
GEOM_TO     = "lower poles milled half (D-shape), upper poles filled"
GEOM_NOTE_FROM = "3 lower FILLED cones (P1/P3/P6)"
GEOM_NOTE_TO   = "3 lower MILLED-HALF cones (P1/P3/P6)"

# Hung-style per-volume ESIZE mesh (replaces global SmartSize).
# Classification reuses the existing VATT material assignment + a VSIZE filter
# to separate small steel cones (V1-V6, ~1e-7 m^3) from the larger yoke body
# (~1e-5 m^3).  Per-volume ESIZE then VMESH that single selection, looped.
MESH_FROM = """ALLSEL,ALL
smrt, MESH_SIZE, , , , , , , , OFF,
mshape, 1, 3D
mshkey, 0
vmesh, all"""

MESH_TO = """ALLSEL,ALL
smrt, MESH_SIZE, , , , , , , , OFF,
mshape, 1, 3D
mshkey, 0
vmesh, all"""

# [REVERTED 2026-05-28] Coil winding stays CCW from above (original source).
# After user clarification, the project convention is CCW current → moment +z →
# flux exits the tip toward WP (post Bz UP, cone Bx toward tip).  Earlier
# attempt to swap to CW was wrong direction; reverted to the CCW source.

count = 0
for pole_n in range(1, 7):
    src = SRC_DIR / f"MT_Sim_P{pole_n}.txt"
    dst = DST_DIR / f"MT_Sim_P{pole_n}.txt"
    if not src.exists():
        raise FileNotFoundError(src)
    text = src.read_text(encoding="utf-8")

    # 1) /CWD path
    text = text.replace("long2016_hexapole_full", "long2016_hexapole_halfcut")

    # 2) Header / banner
    text = text.replace(HEADER_FROM, HEADER_TO)
    text = text.replace(GEOM_FROM, GEOM_TO)
    text = text.replace(GEOM_NOTE_FROM, GEOM_NOTE_TO)

    # 3) VSBV insert after lower-pole VROTAT
    if FILLED_MARKER not in text:
        raise RuntimeError(f"FILLED_MARKER not found in {src.name}; cannot patch")
    text = text.replace(FILLED_MARKER, HALFCUT_REPLACEMENT)

    # 4) Replace global SmartSize with Hung-style per-volume ESIZE
    if MESH_FROM not in text:
        raise RuntimeError(f"MESH_FROM not found in {src.name}; cannot patch mesh")
    text = text.replace(MESH_FROM, MESH_TO)

    # 5) Coil winding stays CCW from above — no patch (reverted 2026-05-28)

    dst.write_text(text, encoding="utf-8")
    count += 1
    print(f"Wrote {dst}")

print(f"\nDone: {count} files generated under {DST_DIR}")
