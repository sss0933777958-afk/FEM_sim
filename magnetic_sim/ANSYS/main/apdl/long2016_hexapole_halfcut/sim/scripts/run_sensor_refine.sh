#!/usr/bin/env bash
# run_sensor_refine.sh — sensor 加密收斂研究：每級 mesh(6 幾何球內細 ESIZE) + 解 coil1。
#   用法： bash run_sensor_refine.sh 1 2 3          # mesh + solve
#          MESH_ONLY=1 bash run_sensor_refine.sh 1  # 只 mesh（驗 sliver，先測 lv1）
#   每級（LV=1/2/3 → 球內 ESIZE = 0.10/0.05/0.03mm）：
#     ① mesh deck sed XLVX→LV、XESX→ESIZE → 產 db/mesh/sensor_refine/mesh_lvLV.db
#     ② solve deck sed XLVX→LV → RESUME mesh_lvLV + coil1=+1 + magsolv + 抽 all/wp
#        → data/sensor_refine/lvLV/coil1 + db/sim/sensor_refine/coil1_lvLV.db
#     ③ half-clean。coil↔電流照 long2016（coil1=+1、map [1,3,6,5,2,4]）。
#   幾何球路線（EREFINE 死路後改）：6 R0.3 純空氣球建進幾何、球內 ESIZE 生細、其餘 smrt5。
set -u
MESH_ONLY="${MESH_ONLY:-0}"
MAPDL="G:/ANSYS Inc/v252/ansys/bin/winx64/MAPDL.exe"
ROOT="G:/my_workspace/code/FEM_sim/magnetic_sim/ANSYS/main"
MDECK="$ROOT/apdl/long2016_hexapole_halfcut/mesh/MT_Mesh_SensorRefine.txt"
SDECK="$ROOT/apdl/long2016_hexapole_halfcut/sim/sensor_refine/MT_Sim_SensorRefine.txt"
MCWD="$ROOT/ANSYS_data/long2016_hexapole_halfcut/db/mesh/sensor_refine"
SCWD="$ROOT/ANSYS_data/long2016_hexapole_halfcut/db/sim/sensor_refine"
DATA="$ROOT/ANSYS_data/long2016_hexapole_halfcut/data/sensor_refine"
mkdir -p "$MCWD" "$SCWD"

# LV → 球內 ESIZE（公尺）
esize_for() { case "$1" in 1) echo "0.10e-3";; 2) echo "0.05e-3";; 3) echo "0.03e-3";; *) echo "";; esac; }

for LV in "$@"; do
  echo "==================== 級數 lv$LV ===================="
  ES=$(esize_for "$LV")
  if [ -z "$ES" ]; then echo "[lv$LV] ABORT: 未定義的級數（只支援 1/2/3）"; continue; fi
  # ---- ① mesh（若 mesh_lvLV.db 已存在則重用，除非 REMESH=1）----
  if [ -f "$MCWD/mesh_lv$LV.db" ] && [ "${REMESH:-0}" != "1" ]; then
    echo "[lv$LV] 重用既有 mesh_lv$LV.db（REMESH=1 可強制重建）"
  else
    tmpm="$MCWD/_mesh_lv$LV.txt"
    sed -e "s/XLVX/$LV/g" -e "s/XESX/$ES/g" "$MDECK" > "$tmpm"
    if ! grep -q "ESIZE_SENS = $ES" "$tmpm" || ! grep -q "FILNAME,mesh_lv$LV" "$tmpm"; then
      echo "[lv$LV] ABORT mesh: sed XLVX/XESX 未命中"; rm -f "$tmpm"; continue
    fi
    echo "[lv$LV] meshing (smrt5 + cone EREFINE + 6 sensor 球內 ESIZE=$ES) ..."
    "$MAPDL" -b -np 4 -m 24000 -dir "$MCWD" -j "mesh_lv$LV" -i "$tmpm" -o "$MCWD/mesh_lv$LV.out"
    nn=$(grep -iE "TOTAL  nodes" "$MCWD/mesh_lv$LV.out" 2>/dev/null | tail -1)
    nsens=$(grep -iE "SENSOR 6-sphere" "$MCWD/mesh_lv$LV.out" 2>/dev/null | tail -1)
    nerrm=$(grep -c "\*\*\* ERROR" "$MCWD/mesh_lv$LV.out" 2>/dev/null || echo "?")
    echo "[lv$LV] mesh done | ERROR=$nerrm | $nn | $nsens"
    rm -f "$MCWD/mesh_lv${LV}_"*.rmg "$MCWD/mesh_lv$LV"*.esav "$MCWD/mesh_lv$LV"*.page* "$MCWD/mesh_lv$LV"*.err "$MCWD/mesh_lv$LV"*.stat "$MCWD/mesh_lv${LV}_"*.out "$MCWD/mesh_lv${LV}_"*.log "$MCWD/_mesh_lv$LV.txt" 2>/dev/null
  fi
  if [ ! -f "$MCWD/mesh_lv$LV.db" ]; then echo "[lv$LV] ABORT: mesh_lv$LV.db 沒產出"; continue; fi
  # mesh-only 清中間檔即結束（驗 sliver 用）
  if [ "$MESH_ONLY" = "1" ]; then
    rm -f "$MCWD/mesh_lv${LV}_"*.rmg "$MCWD/mesh_lv$LV"*.esav "$MCWD/mesh_lv$LV"*.page* "$MCWD/mesh_lv$LV"*.err "$MCWD/mesh_lv$LV"*.stat "$MCWD/_mesh_lv$LV.txt" 2>/dev/null
    echo "[lv$LV] MESH_ONLY：略過 solve"; continue
  fi

  # ---- ② solve coil1 ----
  mkdir -p "$DATA/lv$LV/coil1"
  tmps="$SCWD/_solve_lv$LV.txt"
  sed "s/XLVX/$LV/g" "$SDECK" > "$tmps"
  if ! grep -q "mesh_lv$LV" "$tmps" || ! grep -q "sensor_refine\\\\lv$LV\\\\coil1" "$tmps"; then
    echo "[lv$LV] ABORT solve: sed XLVX 未命中"; rm -f "$tmpm" "$tmps"; continue
  fi
  echo "[lv$LV] solving coil1 ..."
  "$MAPDL" -b -np 4 -m 24000 -dir "$SCWD" -j "coil1_lv$LV" -i "$tmps" -o "$SCWD/solve_lv$LV.out"
  nerr=$(grep -c "\*\*\* ERROR" "$SCWD/solve_lv$LV.out" 2>/dev/null || echo "?")
  ball=$(wc -l < "$DATA/lv$LV/coil1/coil1_bfield_all.dat" 2>/dev/null || echo 0)
  simdb=$( [ -f "$SCWD/coil1_lv$LV.db" ] && echo Y || echo N )
  echo "[lv$LV] solve done | ERROR=$nerr | bfield_all=$ball 行 | sim.db=$simdb"

  # ---- ③ half-clean（保 .db+主.rmg+主 log；刪 worker/中間檔）----
  rm -f "$MCWD/mesh_lv${LV}_"*.rmg "$MCWD/mesh_lv$LV"*.esav "$MCWD/mesh_lv$LV"*.page* "$MCWD/mesh_lv$LV"*.err "$MCWD/mesh_lv$LV"*.stat "$MCWD/_mesh_lv$LV.txt" 2>/dev/null
  rm -f "$SCWD/coil1_lv${LV}_"*.rmg "$SCWD/coil1_lv$LV"*.esav "$SCWD/coil1_lv$LV"*.full "$SCWD/coil1_lv$LV"*.DSP* "$SCWD/coil1_lv$LV"*.page* "$SCWD/coil1_lv$LV"*.stat "$SCWD/coil1_lv$LV"*.err "$SCWD/coil1_lv$LV"*.ldhi "$SCWD/coil1_lv$LV"*.BCS "$SCWD/coil1_lv${LV}_"*.out "$SCWD/coil1_lv${LV}_"*.log "$SCWD/_solve_lv$LV.txt" 2>/dev/null
done
echo "==================== 全部完成 ===================="