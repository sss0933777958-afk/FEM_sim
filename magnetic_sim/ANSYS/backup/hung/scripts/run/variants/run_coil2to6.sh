#!/bin/bash
# Run Coil2-6 simulations sequentially (single process, no -np)
# Same parameters as successful Coil1 run

MAPDL="C:/Program Files/ANSYS2025R2/v252/ansys/bin/winx64/MAPDL.exe"
BASE="C:/Users/pmero/Documents/Lab406/FEM_sim/magnetic_sim/hung"
LOG="$BASE/results/run_coil2to6.log"

echo "" > "$LOG"

for i in 2 3 4 5 6; do
  echo "=== Starting Coil${i} === $(date)" | tee -a "$LOG"

  "$MAPDL" -b -m 8000 \
    -dir "$BASE/results/coil${i}/filleted" \
    -j "coil1" \
    -i "$BASE/apdl/sim/MT_Hung_Simulate_Coil${i}_filleted.txt" \
    -o "$BASE/results/coil${i}/filleted/solve.out" \
    >> "$LOG" 2>&1

  echo "=== Coil${i} done === $(date)" | tee -a "$LOG"
done

echo "=== ALL DONE === $(date)" | tee -a "$LOG"
