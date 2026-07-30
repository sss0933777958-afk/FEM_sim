#!/usr/bin/env bash
# run_graded_all.sh — 重解 long2016 graded / tip 變體 指定 coil，統一補上全域 all 場。
#   用法： bash run_graded_all.sh 1 2 3 4 5 6                    # baseline graded（VARIANT 預設 graded）
#          VARIANT=tip400um bash run_graded_all.sh 1 2 3 4 5 6   # tip-radius 變體（RESUME db/mesh/tip400um）
#   每顆：sed XVARIANTX→$VARIANT + XCOILX→N → 命中檢查 → MAPDL -b RESUME mesh_$VARIANT.db + 施電流
#         + magsolv + 抽 all/wp .dat + SAVE solve db → half-clean（保 sim_coilN.db + 主 .rmg + solve_coilN.out）。
#   coil↔pole = long2016 [1,3,6,5,2,4]；電流 coilN=+1（照 basegap）；all-source 在 MATLAB 端。
#   ⚠ VARIANT=graded 時輸出路徑/RESUME 與原 graded 完全等價（baseline 不變）。
set -u
VARIANT="${VARIANT:-graded}"
MAPDL="G:/ANSYS Inc/v252/ansys/bin/winx64/MAPDL.exe"
ROOT="G:/my_workspace/code/FEM_sim/magnetic_sim/ANSYS/main"
DECK="$ROOT/apdl/long2016_hexapole_halfcut/sim/graded/MT_Sim_graded_all.txt"
SIMCWD="$ROOT/ANSYS_data/long2016_hexapole_halfcut/db/sim/$VARIANT"
DATADIR="$ROOT/ANSYS_data/long2016_hexapole_halfcut/data/$VARIANT"
MESHDB="$ROOT/ANSYS_data/long2016_hexapole_halfcut/db/mesh/$VARIANT/mesh_$VARIANT.db"
mkdir -p "$SIMCWD"
echo "VARIANT=$VARIANT | RESUME mesh=$MESHDB"
if [ ! -f "$MESHDB" ]; then echo "[ABORT] 找不到 mesh db：$MESHDB（先跑 mesh）"; exit 1; fi

for N in "$@"; do
  echo "==================== coil $N ===================="
  mkdir -p "$DATADIR/coil$N"
  tmp="$SIMCWD/_deck_coil$N.txt"
  sed "s/XVARIANTX/$VARIANT/g; s/XCOILX/$N/g" "$DECK" > "$tmp"

  # --- 命中檢查：確認 sed 有換到「激發行」「輸出檔名」「variant 路徑」，否則中止（防呆）---
  if ! grep -q "^R,$N,1,TURNS\*1," "$tmp" || ! grep -q "coil${N}_bfield_all" "$tmp"; then
    echo "[coil$N] ABORT: sed 未正確替換 XCOILX（R,$N,1,TURNS*1 或 coil${N}_bfield_all 未命中）"; rm -f "$tmp"; continue
  fi
  if ! grep -q "db\\\\mesh\\\\$VARIANT" "$tmp" || ! grep -q "data\\\\$VARIANT\\\\coil$N" "$tmp"; then
    echo "[coil$N] ABORT: sed 未正確替換 XVARIANTX（db\\mesh\\$VARIANT 或 data\\$VARIANT\\coil$N 未命中）"; rm -f "$tmp"; continue
  fi
  # 防呆：殘留的 XCOILX / XVARIANTX 表示替換不完全
  if grep -qE "XCOILX|XVARIANTX" "$tmp"; then
    echo "[coil$N] ABORT: 仍有殘留 XCOILX/XVARIANTX"; rm -f "$tmp"; continue
  fi

  echo "[coil$N] solving ..."
  "$MAPDL" -b -np 4 -m 24000 -dir "$SIMCWD" -j "sim_coil$N" -i "$tmp" -o "$SIMCWD/solve_coil$N.out"

  # --- 驗證：0 error + 四個 .dat 落點 ---
  nerr=$(grep -c "\*\*\* ERROR" "$SIMCWD/solve_coil$N.out" 2>/dev/null || echo "?")
  ndat=$(ls "$DATADIR/coil$N/coil${N}_"{coord,bfield}_{all,wp}.dat 2>/dev/null | wc -l)
  ball=$(wc -l < "$DATADIR/coil$N/coil${N}_bfield_all.dat" 2>/dev/null || echo 0)
  bwp=$(wc -l < "$DATADIR/coil$N/coil${N}_bfield_wp.dat" 2>/dev/null || echo 0)
  db=$( [ -f "$SIMCWD/sim_coil$N.db" ] && echo Y || echo N )
  echo "[coil$N] done | ERROR=$nerr | .dat=$ndat/4 | bfield_all=$ball 行 | bfield_wp=$bwp 行 | solve.db=$db"

  # --- half-clean：只刪該 coil 的 worker/中間檔，保 .db + 主 .rmg + solve log ---
  rm -f "$SIMCWD/sim_coil${N}_"*.rmg "$SIMCWD/sim_coil$N"*.esav "$SIMCWD/sim_coil$N"*.full \
        "$SIMCWD/sim_coil$N"*.DSP* "$SIMCWD/sim_coil$N"*.page* "$SIMCWD/sim_coil$N"*.stat \
        "$SIMCWD/sim_coil$N"*.err "$SIMCWD/sim_coil$N"*.ldhi "$SIMCWD/sim_coil$N"*.mlv \
        "$SIMCWD/sim_coil$N"*.BCS "$SIMCWD/_deck_coil$N.txt" 2>/dev/null
done
echo "==================== 全部完成 ===================="