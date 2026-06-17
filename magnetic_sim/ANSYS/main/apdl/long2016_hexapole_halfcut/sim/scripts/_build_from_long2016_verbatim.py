"""Build 6 halfcut sims by copying Long2016 source verbatim + minimal kuo tweaks.

Long2016 source (magnetic_sim/hexapole-long2016/apdl/MT_Modeling_Geometry_Meshing_Solving_Coil1.txt)
already IS the halfcut version with VADD merges (yoke+prots, cone+yoke per pole) which
the kuo source `long2016_hexapole_full/MT_Sim_P*.txt` had stripped out — that strip
caused the topology break that produced 30-70% weaker flux to other poles vs Long2016.

This script bypasses both the kuo "full" source and the halfcut generator by writing
6 sim files directly from Long2016 source with 3 minimal modifications per file:

1. /CWD line points to magnetic_sim/ANSYS/main/ANSYS_data/long2016_hexapole_halfcut/coil{N}/
2. CURR_ARRAY: only coil{N} active for sim N (=1), others 0
3. POST1 export at end replaced with working kuo POST1 (nlist/prnsol → .dat files)

After running this, no need to run _generate_halfcut_sims.py for baseline.
gap200um generator (_generate_halfcut_sims_gap.py) still reads from these halfcut sims.

Usage:
    python _build_from_long2016_verbatim.py
"""
from pathlib import Path

LONG2016_SRC = Path(r"G:\my_workspace\code\FEM_sim\magnetic_sim\hexapole-long2016\apdl\MT_Modeling_Geometry_Meshing_Solving_Coil1.txt")
KUO_HALFCUT_DIR = Path(r"G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\apdl\long2016_hexapole_halfcut\sim")

# Kuo POST1 export block (verbatim from current MT_Sim_P1.txt end)
POST1_KUO = """
! ====================================================================
! POST1: export 4 .dat files (kuo standard, matches gen_Vout_Vin_4p572.m)
! ====================================================================
/POST1
SET,LAST
ALLSEL,ALL

/HEADER,OFF,OFF,OFF,OFF,OFF,OFF
/FORMAT,10,,20,12,,
/PAGE,1000000,,1000000

/OUTPUT,'coil{n}_coord_all','dat'
NLIST,ALL, , ,COORD
/OUTPUT

/OUTPUT,'coil{n}_bfield_all','dat'
PRNSOL,B,COMP
/OUTPUT

LOCAL,12,2, 0, 0, SPH_OFST
CSYS,12
NSEL,S,LOC,X, 0, 2e-3
CSYS,0

/OUTPUT,'coil{n}_coord_wp','dat'
NLIST,ALL, , ,COORD
/OUTPUT

/OUTPUT,'coil{n}_bfield_wp','dat'
PRNSOL,B,COMP
/OUTPUT

ALLSEL,ALL
FINISH
"""

# Long2016 source's terminal /POST1 block (commented PRNSOL,B) — replaced wholesale
LONG2016_POST1_TAIL = """/POST1
/FORMAT,7,,17,9,1000000,
/PAGE,1000000,,1000000
!NLIST,ALL,,,COORD
!PRNSOL,B
"""

def build(pole_n: int) -> None:
    text = LONG2016_SRC.read_text(encoding="utf-8")

    # 1. /CWD point to kuo halfcut result dir (line 3 has pmero path)
    cwd_pattern_pmero = "/CWD,'C:\\Users\\pmero\\Documents\\Lab406\\FEM_sim\\hexapole-long2016\\results\\coil1'   ! [MODIFIED] output to magnetic_sim/hexapole-long2016/results/coil1"
    cwd_kuo = f"/CWD,'G:\\my_workspace\\code\\FEM_sim\\kuo\\results\\long2016_hexapole_halfcut\\coil{pole_n}'   ! [KUO] halfcut sim output"
    if cwd_pattern_pmero not in text:
        raise RuntimeError("/CWD pmero path marker not found")
    text = text.replace(cwd_pattern_pmero, cwd_kuo)

    # 2. CURR_ARRAY: only coil{pole_n} = 1, others = 0
    #    Long2016 source has CURR_ARRAY(1) = 1, (2..6) = 0
    for k in range(1, 7):
        old_val = "1" if k == 1 else "0"
        new_val = "1" if k == pole_n else "0"
        old_line = f"CURR_ARRAY({k}) = {old_val}"
        new_line = f"CURR_ARRAY({k}) = {new_val}"
        if old_line not in text:
            raise RuntimeError(f"CURR_ARRAY({k}) = {old_val} marker not found")
        text = text.replace(old_line, new_line, 1)

    # 3. Replace terminal POST1 block with kuo POST1 export
    if LONG2016_POST1_TAIL not in text:
        raise RuntimeError("Long2016 terminal POST1 marker not found")
    text = text.replace(LONG2016_POST1_TAIL, POST1_KUO.format(n=pole_n))

    # 4. Sanity-check no kuo-source-specific artefact left (mat_assign was kuo only)
    assert "mat_assign" not in text, "mat_assign should not appear (Long2016 uses explicit VATT)"

    dst = KUO_HALFCUT_DIR / f"MT_Sim_P{pole_n}.txt"
    dst.write_text(text, encoding="utf-8")
    print(f"Wrote {dst}")


if __name__ == "__main__":
    for n in range(1, 7):
        build(n)
    print(f"\nDone: 6 halfcut sims built from Long2016 verbatim source.")
