#!/usr/bin/env bash
# run_ws_sim.sh — 對 hung 工作空間變體 mesh 跑 6-coil 磁靜 FEM 求解（交付 sim/data/db/sim）
# ---------------------------------------------------------------------------
# 從 canonical solve deck `../gap_200um/MT_Sim_P1.txt` 生每 (變體,coil) 暫存 deck：
#   RESUME → db/mesh/mesh_graded_R<tag>（我的變體 mesh）；激發 coil N（R,N=TURNS*1）；
#   D,ALL,MAG,0 @R180 → magsolv,3 DSP → 抽 coilN_{coord,bfield}_{all,wp}.dat →
#   SAVE sim_coilN.db。輸出 db/sim/R<tag>/ + data/R<tag>/coilN/。跑完 half-clean。
# 用法：  bash run_ws_sim.sh R700 R300        # 每變體 6 coil
# ---------------------------------------------------------------------------
set -u
MAPDL="G:\\ANSYS Inc\\v252\\ansys\\bin\\winx64\\MAPDL.exe"
ROOT="G:/my_workspace/code/FEM_sim/magnetic_sim/ANSYS/main"
DECK="$ROOT/apdl/hung_hexapole/sim/gap_200um/MT_Sim_P1.txt"
SIMDIR="$ROOT/apdl/hung_hexapole/sim"

for tag in "$@"; do
  meshdb_w="G:\\my_workspace\\code\\FEM_sim\\magnetic_sim\\ANSYS\\main\\ANSYS_data\\hung_hexapole\\db\\mesh\\mesh_graded_${tag}"
  simcwd_u="$ROOT/ANSYS_data/hung_hexapole/db/sim/${tag}"
  simcwd_w="G:\\my_workspace\\code\\FEM_sim\\magnetic_sim\\ANSYS\\main\\ANSYS_data\\hung_hexapole\\db\\sim\\${tag}"
  datadir="$ROOT/ANSYS_data/hung_hexapole/data/${tag}"
  mkdir -p "$simcwd_u"
  for N in 1 2 3 4 5 6; do
    mkdir -p "$datadir/coil${N}"
    tmp="$SIMDIR/_ws_${tag}_coil${N}.txt"
    python - "$DECK" "$tmp" "$tag" "$N" "$meshdb_w" "$simcwd_w" <<'PY'
import sys
deck,out,tag,N,meshdb_w,simcwd_w = sys.argv[1:7]
t=open(deck,encoding='utf-8').read()
# CWD -> db/sim/<tag>
t=t.replace(r"\db\sim\gap_200um_solve", r"\db\sim\%s"%tag)
# RESUME -> my variant mesh
t=t.replace(r"db\geom\mesh_graded_basegap", r"db\mesh\mesh_graded_%s"%tag)
# energise coil N (canonical energises coil1 at 'R,1,1,TURNS*1')
t=t.replace("R,1,1,TURNS*1,COIL_THK,COIL_H", "R,%s,1,TURNS*1,COIL_THK,COIL_H"%N)
# output file/dir coil1 -> coilN ; data\gap_200um\coil1 -> data\<tag>\coilN
t=t.replace(r"\data\hung_hexapole\data\gap_200um\coil1", r"\data\hung_hexapole\data\%s\coil%s"%(tag,N))
for f in ("coord_all","bfield_all","coord_wp","bfield_wp"):
    t=t.replace("'coil1_%s'"%f, "'coil%s_%s'"%(N,f))
t=t.replace("coil 1", "coil %s"%N)
# SAVE solved db before exit (deliver db/sim)
t=t.replace("/EXIT,NOSAVE", "SAVE,'sim_coil%s','db'\n/EXIT,NOSAVE"%N)
open(out,'w',encoding='utf-8').write(t)
PY
    echo "[$tag coil$N] solving ..."
    ( cd "$simcwd_u" && "$MAPDL" -b -np 4 -m 24000 -dir "$simcwd_w" -j "sim_coil${N}" -i "$tmp" -o "$simcwd_u/solve_coil${N}.out" >/dev/null 2>&1 )
    # half-clean this coil (keep .db + 主 .rmg + solve log; drop worker/esav/etc)
    find "$simcwd_u" -maxdepth 1 -type f \( -name "sim_coil${N}_*.rmg" -o -name "sim_coil${N}*.esav" -o -name "sim_coil${N}*.full" -o -name "sim_coil${N}*.DSP*" -o -name "sim_coil${N}*.page*" -o -name "sim_coil${N}*.err" -o -name "sim_coil${N}*.stat" -o -name "sim_coil${N}*.ldhi" \) -delete 2>/dev/null
    rm -f "$tmp"
    nn=$(grep -ic "coil${N}_bfield_wp" "$simcwd_u/solve_coil${N}.out" 2>/dev/null)
    err=$(grep -c "ERROR" "$simcwd_u/solve_coil${N}.out" 2>/dev/null)
    echo "[$tag coil$N] done (db=$([ -f "$simcwd_u/sim_coil${N}.db" ] && echo Y || echo N), errors=$err)"
  done
done
