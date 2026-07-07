#!/bin/bash
# run_sweep.sh — single-pole TIPCUT geometry sweep driver (cut-plane x from 500..550 um).
#   Each case: sed __XCUT__/__OUTDIR__ into template -> one MAPDL run (geom+mesh+solve+extract)
#   -> verify 4 .dat -> FULL-CLEAN scratch (rm swp.*; only .dat in data dir + run log survive).
# Usage:  ./run_sweep.sh [um ...]     (default 500 501 ... 550 ; pilot e.g. ./run_sweep.sh 500 525 550)
set -u
ANSYS="G:/ANSYS Inc/v252/ansys/bin/winx64/MAPDL.exe"
BASE="G:/my_workspace/code/FEM_sim/magnetic_sim/ANSYS/main"
TEMPLATE="$BASE/apdl/long2016_hexapole_halfcut/sim/singlepole_tipcut_sweep/MT_Sweep_Tipcut_template.txt"
SCRATCH="$BASE/ANSYS_data/long2016_hexapole_halfcut/db/singlepole/tipcut_sweep"
DATA_UNIX="$BASE/ANSYS_data/long2016_hexapole_halfcut/data/coil1/singlepole/tipcut_sweep"
DATA_WIN='G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data\coil1\singlepole\tipcut_sweep'
DECKS="$SCRATCH/_decks"
mkdir -p "$SCRATCH" "$DECKS" "$DATA_UNIX"

UMS="$@"
if [ -z "$UMS" ]; then UMS=$(seq 500 10 550); fi

for um in $UMS; do
  xcut="0.$(printf '%03d' "$um")e-3"                       # 500 -> 0.500e-3 (=500um)
  outwin="${DATA_WIN}\\tipcut_x${um}"
  outunix="${DATA_UNIX}/tipcut_x${um}"
  mkdir -p "$outunix"
  outesc=$(printf '%s' "$outwin" | sed 's/\\/\\\\/g')      # escape backslashes for sed replacement
  deck="$DECKS/tipcut_x${um}.txt"
  sed -e "s|__XCUT__|${xcut}|" -e "s|__OUTDIR__|${outesc}|g" "$TEMPLATE" > "$deck"

  echo "=== x=${um}um  XCUT=${xcut}  $(date +%H:%M:%S) ==="
  "$ANSYS" -b -np 4 -m 24000 -dir "$SCRATCH" -j swp -i "$deck" -o "$SCRATCH/run_x${um}.out"

  n=$(ls "$outunix"/coord_all.dat "$outunix"/bfield_all.dat "$outunix"/coord_tip.dat "$outunix"/bfield_tip.dat 2>/dev/null | wc -l)
  nn=$(grep -o "MESH DONE.*elems" "$SCRATCH/run_x${um}.out" 2>/dev/null | tail -1)
  err=$(grep -c "ERROR" "$SCRATCH/run_x${um}.out" 2>/dev/null)
  echo "  -> ${n}/4 .dat | ${nn} | ERROR lines=${err}"

  rm -f "$SCRATCH"/swp* "$SCRATCH"/magsolv.out "$SCRATCH"/scratch 2>/dev/null   # FULL-CLEAN: swp.db + swp0.* DMP workers (.rmg/.esav/.full)
done
echo "DONE sweep: $UMS  $(date +%H:%M:%S)"
