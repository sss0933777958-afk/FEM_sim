"""Build 6 gap100um_mueq sim decks (coil1..6) from the proven gap200um_mueq decks.

  gap = 50 um each end (top+bottom) = 100 um total in series with 7 mm steel post.
  mu_eff = h / [ (h-g)/mu_steel + g ] = 7 / [ (7-0.1)/280 + 0.1 ] = 56.

  Changes vs gap200um deck (per coil N=1..6):
    1. mp,murx,MAT_PROT  31 -> 56
    2. comments 200um/mu_r=31 -> 100um(50um each end)/mu_r=56
    3. /CWD  ...\coilN\gap200um_mueq'  ->  ...\db\gap100um\coilN'
    4. insert a "! [MESH] baseline smrt5" header after the active /CWD line
  Mesh = baseline smrt5 (unchanged); EMODIF 6 protrusions -> MAT_PROT (unchanged); per-coil CURR_ARRAY kept.

Run:  python _build_gap100um_mu_eq.py
"""
import pathlib

ROOT = pathlib.Path(r"G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\apdl\long2016_hexapole_halfcut\sim")
SRC  = ROOT / "gap200um_mueq"
DST  = ROOT / "gap100um_mueq"
DST.mkdir(exist_ok=True)

MESH_HDR = "! [MESH] baseline - smrt 5 free-tet SOLID96 (mu_r-equivalent gap, geometry/mesh unchanged)"

for N in range(1, 7):
    src = SRC / f"MT_Sim_P{N}_gap200um_mueq.txt"
    t = src.read_text(encoding="utf-8")

    # 1. mu_r value
    t = t.replace("mp , murx , MAT_PROT , 31", "mp , murx , MAT_PROT , 56")
    # 2. comments
    t = t.replace("for 200µm gap", "for 100µm gap (50µm each end)")
    t = t.replace("7mm steel + 200µm gap series", "7mm post + 100µm gap (50µm each end) series")
    t = t.replace("(μ_r=31)", "(μ_r=56)")
    t = t.replace("(μ_r equivalent for 200", "(μ_r equivalent for 100")
    # 3. CWD path  \coilN\gap200um_mueq'  ->  \db\gap100um\coilN'
    t = t.replace(f"\\coil{N}\\gap200um_mueq'", f"\\db\\gap100um\\coil{N}'")

    # 4. insert MESH header after the active /CWD line (not the commented !****/CWD ones)
    out = []
    for line in t.split("\n"):
        out.append(line)
        if line.startswith("/CWD,"):
            out.append(MESH_HDR)
    t = "\n".join(out)

    dst = DST / f"MT_Sim_P{N}_gap100um_mueq.txt"
    dst.write_text(t, encoding="utf-8")
    has56 = "MAT_PROT , 56" in t
    hascwd = f"\\db\\gap100um\\coil{N}'" in t
    print(f"P{N}: wrote {dst.name}  mu_r=56:{has56}  CWD-ok:{hascwd}")

print("done.")
