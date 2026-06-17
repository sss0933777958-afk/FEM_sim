#!/bin/bash
# Run POST1 data extraction for Coil2-6 (resume db, export .dat files)

MAPDL="C:/Program Files/ANSYS2025R2/v252/ansys/bin/winx64/MAPDL.exe"
BASE="C:/Users/pmero/Documents/Lab406/FEM_sim/magnetic_sim/hung"
LOG="$BASE/results/run_post2to6.log"

echo "" > "$LOG"

for i in 2 3 4 5 6; do
  echo "=== POST Coil${i} === $(date)" | tee -a "$LOG"

  "$MAPDL" -b -m 4000 \
    -dir "$BASE/results/coil${i}/filleted" \
    -j "coil1" \
    -i "$BASE/apdl/postproc/post_export_data_coil${i}.txt" \
    -o "$BASE/results/coil${i}/filleted/post.out" \
    >> "$LOG" 2>&1

  echo "=== POST Coil${i} done === $(date)" | tee -a "$LOG"
done

echo "=== ALL POST DONE === $(date)" | tee -a "$LOG"
