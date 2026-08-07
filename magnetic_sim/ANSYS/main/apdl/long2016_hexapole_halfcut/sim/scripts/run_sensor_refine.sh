#!/usr/bin/env bash
# run_sensor_refine.sh — sensor 加密收斂研究：每級 mesh(幾何球內細 ESIZE) + 解 $COILS 各顆 coil。
#   用法： bash run_sensor_refine.sh 1 2 3          # mesh + solve
#          MESH_ONLY=1 bash run_sensor_refine.sh 1  # 只 mesh（驗 sliver，先測 lv1）
#   每級（LV=1/2/3 → 球內 ESIZE = 0.10/0.05/0.03mm）：
#     ① mesh deck sed XLVX→LV、XESX→ESIZE → 產 db/mesh/$CASE/mesh_lvLV.db
#     ② solve deck sed XLVX→LV、XCOILX→C → RESUME mesh_lvLV + coilC=+1 + magsolv + 抽 all/wp
#        → data/$CASE/lvLV/coilC + db/sim/$CASE/coilC_lvLV.db
#     ③ half-clean。coil↔電流照 long2016（激發顆=+1、map [1,3,6,5,2,4]）。
#   幾何球路線（EREFINE 死路後改）：R0.3 純空氣球建進幾何、球內 ESIZE 生細、其餘 smrt5。
#   球組/球心由 mesh deck 決定（現行 3 顆 @ SOFF 4.5）；case 與 coil 由 CASE / COILS 控制。
set -u
MESH_ONLY="${MESH_ONLY:-0}"
MAPDL="G:/ANSYS Inc/v252/ansys/bin/winx64/MAPDL.exe"
ROOT="G:/my_workspace/code/FEM_sim/magnetic_sim/ANSYS/main"
MDECK="$ROOT/apdl/long2016_hexapole_halfcut/mesh/MT_Mesh_SensorRefine.txt"
SDECK="$ROOT/apdl/long2016_hexapole_halfcut/sim/sensor_refine/MT_Sim_SensorRefine.txt"
# [ADDED 2026-08-05] CASE：加密球組的 case tag。
#   預設 sensor_refine_p1p2 = 3 顆球 @ SOFF 4.5（P1 削平面 / P1 底錐面 / P2 上錐面），
#   對應 mesh/MT_Mesh_SensorRefine.txt 現行內容（其 /CWD 也指向同一個 case）。
#   舊 case sensor_refine（6 顆 @ 4.572）已封存，要重跑得同時把 deck 的球心與 /CWD 改回去。
CASE="${CASE:-sensor_refine_p1p2}"
# [ADDED 2026-08-05] COILS：要解哪幾顆 coil（APDL 索引）。
#   coil1 = P1 自激、coil5 = P2 自激（long2016 map [1,3,6,5,2,4]；⚠ coil2 是 P3，不是 P2）。
#   例： COILS="1" bash run_sensor_refine.sh 3      # 只解 lv3 的 coil1
COILS="${COILS:-1 5}"
MCWD="$ROOT/ANSYS_data/long2016_hexapole_halfcut/db/mesh/$CASE"
SCWD="$ROOT/ANSYS_data/long2016_hexapole_halfcut/db/sim/$CASE"
DATA="$ROOT/ANSYS_data/long2016_hexapole_halfcut/data/$CASE"
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
    echo "[lv$LV] meshing (smrt5 + cone EREFINE + sensor 球內 ESIZE=$ES) ..."
    "$MAPDL" -b -np 4 -m 24000 -dir "$MCWD" -j "mesh_lv$LV" -i "$tmpm" -o "$MCWD/mesh_lv$LV.out"
    nn=$(grep -iE "TOTAL  nodes" "$MCWD/mesh_lv$LV.out" 2>/dev/null | tail -1)
    nsens=$(grep -iE "SENSOR sphere" "$MCWD/mesh_lv$LV.out" 2>/dev/null | tail -1)
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

  # ---- ② solve：對 $COILS 每一顆 coil 各解一次 ----
  # [MODIFIED 2026-08-05] 原本寫死 coil1，改成內層迴圈（deck 的 XCOILX 由此 sed）。
  for C in $COILS; do
    mkdir -p "$DATA/lv$LV/coil$C"
    tmps="$SCWD/_solve_lv${LV}_c$C.txt"
    sed -e "s/XLVX/$LV/g" -e "s/XCOILX/$C/g" "$SDECK" > "$tmps"
    # 守衛：確認 sed 真的命中（網格名 + 本 case 的輸出路徑 + 激發那一行）。
    #   若 solve deck 的路徑還指向舊 case，或 XCOILX 沒被換掉 → ABORT，
    #   而不是把新 case 的結果默默寫進舊資料夾。
    # ⚠ 路徑含 Windows 反斜線 → **必須用 grep -F**（固定字串）。
    #   原本的 grep -q "…\\\\lv…\\\\coil…" 在此 bash/grep 組合下**永遠不匹配**（實測），
    #   會讓每一次 solve 都誤判成「sed 未命中」而 ABORT。
    if ! grep -q "mesh_lv$LV" "$tmps" \
       || ! grep -qF "$CASE\\lv$LV\\coil$C" "$tmps" \
       || ! grep -q "^R,$C,1,TURNS\*1," "$tmps"; then
      echo "[lv$LV coil$C] ABORT solve: sed XLVX/XCOILX 未命中"; rm -f "$tmps"; continue
    fi
    echo "[lv$LV coil$C] solving ..."
    "$MAPDL" -b -np 4 -m 24000 -dir "$SCWD" -j "coil${C}_lv$LV" -i "$tmps" -o "$SCWD/solve_lv${LV}_c$C.out"
    nerr=$(grep -c "\*\*\* ERROR" "$SCWD/solve_lv${LV}_c$C.out" 2>/dev/null || echo "?")
    ball=$(wc -l < "$DATA/lv$LV/coil$C/coil${C}_bfield_all.dat" 2>/dev/null || echo 0)
    simdb=$( [ -f "$SCWD/coil${C}_lv$LV.db" ] && echo Y || echo N )
    echo "[lv$LV coil$C] solve done | ERROR=$nerr | bfield_all=$ball 行 | sim.db=$simdb"

    # ---- ③ half-clean（保 .db+主.rmg+主 log；刪 worker/中間檔）----
    rm -f "$SCWD/coil${C}_lv${LV}_"*.rmg "$SCWD/coil${C}_lv$LV"*.esav "$SCWD/coil${C}_lv$LV"*.full "$SCWD/coil${C}_lv$LV"*.DSP* "$SCWD/coil${C}_lv$LV"*.page* "$SCWD/coil${C}_lv$LV"*.stat "$SCWD/coil${C}_lv$LV"*.err "$SCWD/coil${C}_lv$LV"*.ldhi "$SCWD/coil${C}_lv$LV"*.BCS "$SCWD/coil${C}_lv${LV}_"*.out "$SCWD/coil${C}_lv${LV}_"*.log "$tmps" 2>/dev/null
  done
  rm -f "$MCWD/mesh_lv${LV}_"*.rmg "$MCWD/mesh_lv$LV"*.esav "$MCWD/mesh_lv$LV"*.page* "$MCWD/mesh_lv$LV"*.err "$MCWD/mesh_lv$LV"*.stat "$MCWD/_mesh_lv$LV.txt" 2>/dev/null
done
echo "==================== 全部完成 ===================="