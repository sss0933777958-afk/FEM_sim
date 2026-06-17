#!/bin/bash
# Re-run Coil2-6 with 3-pass NREFINE + POST1 export (more memory)
MAPDL="C:/Program Files/ANSYS2025R2/v252/ansys/bin/winx64/MAPDL.exe"
BASE="C:/Users/pmero/Documents/Lab406/FEM_sim/magnetic_sim/hung"
LOG="$BASE/results/run_coil2to6_refined_v2.log"

echo "" > "$LOG"

for i in 2 3 4 5 6; do
  echo "=== Coil${i} solve === $(date)" | tee -a "$LOG"

  # Clean old result files to avoid confusion
  rm -f "$BASE/results/coil${i}/filleted/coil1_"*.rmg
  rm -f "$BASE/results/coil${i}/filleted/coil1.rmg"

  "$MAPDL" -b -m 8000 \
    -dir "$BASE/results/coil${i}/filleted" \
    -j "coil1" \
    -i "$BASE/apdl/sim/MT_Hung_Simulate_Coil${i}_filleted.txt" \
    -o "$BASE/results/coil${i}/filleted/solve.out" \
    >> "$LOG" 2>&1

  echo "=== Coil${i} POST1 === $(date)" | tee -a "$LOG"

  "$MAPDL" -b -m 4000 \
    -dir "$BASE/results/coil${i}/filleted" \
    -j "coil1" \
    -i "$BASE/apdl/postproc/post_export_data_coil${i}.txt" \
    -o "$BASE/results/coil${i}/filleted/post.out" \
    >> "$LOG" 2>&1

  echo "=== Coil${i} done === $(date)" | tee -a "$LOG"
done

echo "=== ALL DONE === $(date)" | tee -a "$LOG"
