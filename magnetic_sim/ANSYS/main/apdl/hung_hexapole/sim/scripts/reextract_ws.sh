#!/usr/bin/env bash
# reextract_ws.sh — 從已解好的 db/sim/R<tag>/sim_coilN.{db,rmg} RESUME + POST1 重抽 4 .dat
#   到正確的 data/R<tag>/coilN/（修正 run_ws_sim.sh 的路徑 bug：原本誤寫進 data/gap_200um/）。
#   不重解（SET,LAST 讀 .rmg 的場）。WP=origin、抽 all + WP(2mm) 區 coord/bfield。
#   用法：  bash reextract_ws.sh R700 R300          # 每變體 6 coil
#          bash reextract_ws.sh R700 1              # 單一變體單一 coil（測試）
set -u
MAPDL="G:\\ANSYS Inc\\v252\\ansys\\bin\\winx64\\MAPDL.exe"
ROOT="G:/my_workspace/code/FEM_sim/magnetic_sim/ANSYS/main"
DROOT_W="G:\\my_workspace\\code\\FEM_sim\\magnetic_sim\\ANSYS\\main\\ANSYS_data\\hung_hexapole"

# 若第二個參數是單一數字 → 單 coil 測試
tags=(); coils=(1 2 3 4 5 6)
for a in "$@"; do
  if [[ "$a" =~ ^[1-6]$ ]]; then coils=("$a"); else tags+=("$a"); fi
done

for tag in "${tags[@]}"; do
  simcwd_u="$ROOT/ANSYS_data/hung_hexapole/db/sim/${tag}"
  simcwd_w="${DROOT_W}\\db\\sim\\${tag}"
  for N in "${coils[@]}"; do
    datadir_w="${DROOT_W}\\data\\${tag}\\coil${N}"
    mkdir -p "$ROOT/ANSYS_data/hung_hexapole/data/${tag}/coil${N}"
    tmp="$ROOT/apdl/hung_hexapole/sim/_reext_${tag}_c${N}.txt"
    cat > "$tmp" <<EOF
/CWD,'${simcwd_w}'
/CLEAR,NOSTART
RESUME,'sim_coil${N}','db','${simcwd_w}'
/POST1
SET,LAST
ALLSEL,ALL
/HEADER,OFF,OFF,OFF,OFF,OFF,OFF
/FORMAT,10,,20,12,,
/PAGE,1000000,,1000000
/OUTPUT,'coil${N}_coord_all','dat','${datadir_w}'
NLIST,ALL, , ,COORD
/OUTPUT
/OUTPUT,'coil${N}_bfield_all','dat','${datadir_w}'
PRNSOL,B,COMP
/OUTPUT
LOCAL,12,2,0,0,0
CSYS,12
NSEL,S,LOC,X,0,2e-3
CSYS,0
/OUTPUT,'coil${N}_coord_wp','dat','${datadir_w}'
NLIST,ALL, , ,COORD
/OUTPUT
/OUTPUT,'coil${N}_bfield_wp','dat','${datadir_w}'
PRNSOL,B,COMP
/OUTPUT
ALLSEL,ALL
CSYS,0
FINISH
/EXIT,NOSAVE
EOF
    echo "[$tag coil$N] re-extracting ..."
    ( cd "$simcwd_u" && "$MAPDL" -b -np 4 -m 24000 -dir "$simcwd_w" -j "sim_coil${N}" -i "$tmp" -o "$simcwd_u/reext_coil${N}.out" >/dev/null 2>&1 )
    rm -f "$tmp"
    nb=$(wc -l < "$ROOT/ANSYS_data/hung_hexapole/data/${tag}/coil${N}/coil${N}_bfield_all.dat" 2>/dev/null || echo 0)
    echo "[$tag coil$N] done (bfield_all lines=$nb)"
  done
done
