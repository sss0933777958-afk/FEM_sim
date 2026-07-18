#!/usr/bin/env bash
# ============================================================================
# run_full_sim.sh — NTU_hexapole full_assembly_sleeve：依序跑 coil 1..6 + half-clean
#
#   對 canonical deck sim/full_assembly/MT_Sim_FullAssembly.txt 的 temp copy 只 sed 兩行
#   （COIL_ID / CN），其餘（線圈參數、BC、magsolv、POST1 路徑）一字不動。
#
#   用法：
#     bash run_full_sim.sh            # 跑 coil 1..6
#     bash run_full_sim.sh 1          # 只跑 coil 1（放行前先驗證用）
#     bash run_full_sim.sh 2 3 4      # 跑指定幾顆
#
#   ⚠ sed 取代字串必須逐字對上 deck 現況（hung run_ws_sim.sh 就是這裡出包：
#     把 \data\hung_hexapole\... 當真實路徑、實際是 \ANSYS_data\... → .dat 全誤寫、
#     覆蓋損毀 baseline）。本腳本只 sed 兩行純參數，路徑一律由 deck 自己寫死，
#     從結構上避開該坑；跑完仍會核 .dat 落點。
# ============================================================================
set -u

MAPDL="G:\\ANSYS Inc\\v252\\ansys\\bin\\winx64\\MAPDL.exe"
ROOT_U="/g/my_workspace/code/FEM_sim/magnetic_sim/ANSYS/main"
ROOT_W="G:\\my_workspace\\code\\FEM_sim\\magnetic_sim\\ANSYS\\main"

DECK="$ROOT_U/apdl/NTU_hexapole/sim/full_assembly/MT_Sim_FullAssembly.txt"
SIMCWD_U="$ROOT_U/ANSYS_data/NTU_hexapole/db/sim/full_assembly_sleeve"
SIMCWD_W="$ROOT_W\\ANSYS_data\\NTU_hexapole\\db\\sim\\full_assembly_sleeve"
DATA_U="$ROOT_U/ANSYS_data/NTU_hexapole/data/full_assembly_sleeve"

COILS=( "$@" )
if [ ${#COILS[@]} -eq 0 ]; then COILS=(1 2 3 4 5 6); fi

mkdir -p "$SIMCWD_U"

for N in "${COILS[@]}"; do
    mkdir -p "$DATA_U/coil${N}"
    tmp="$SIMCWD_U/_sim_coil${N}.txt"

    # 只換 active coil 兩行（deck 內 'COIL_ID = <n>' 與 "CN = '<n>'"；兩行後面都帶行尾註解，不可用 $ 錨定）
    sed -e "s/^COIL_ID = [0-9]\+/COIL_ID = ${N}/" \
        -e "s/^CN = '[0-9]'/CN = '${N}'/" \
        "$DECK" > "$tmp"

    # 驗證 sed 真的有換到（避免靜默跑成 coil1）
    if ! grep -q "^COIL_ID = ${N}\([^0-9]\|\$\)" "$tmp" || ! grep -q "^CN = '${N}'" "$tmp"; then
        echo "[coil${N}] ABORT: sed 未命中 COIL_ID / CN 行（deck 格式已變？）"
        rm -f "$tmp"; exit 1
    fi

    echo "[coil${N}] solving ..."
    ( cd "$SIMCWD_U" && "$MAPDL" -b -np 4 -m 32000 \
        -dir "$SIMCWD_W" -j "sim_coil${N}" -i "$tmp" \
        -o "$SIMCWD_U/solve_coil${N}.out" > /dev/null 2>&1 )

    rm -f "$tmp"

    # half-clean（per .claude/rules/db-folder-retention.md）：
    #   留 sim_coilN.db + 主 sim_coilN.rmg + solve_coilN.out
    #   ⚠ 絕不可寫 rm -f sim_coilN*（會連主 .rmg + .db 一起殺）
    rm -f "$SIMCWD_U/sim_coil${N}_"*.rmg \
          "$SIMCWD_U/sim_coil${N}"*.esav \
          "$SIMCWD_U/sim_coil${N}"*.full \
          "$SIMCWD_U/sim_coil${N}"*.DSP* \
          "$SIMCWD_U/sim_coil${N}"*.page* \
          "$SIMCWD_U/sim_coil${N}_"*.err \
          "$SIMCWD_U/sim_coil${N}_"*.out \
          "$SIMCWD_U/sim_coil${N}_"*.log \
          "$SIMCWD_U/sim_coil${N}_"*.stat \
          "$SIMCWD_U/sim_coil${N}"*.ldhi 2>/dev/null

    db="N"; [ -f "$SIMCWD_U/sim_coil${N}.db" ] && db="Y"
    nb=0; [ -f "$DATA_U/coil${N}/full_assembly_sleeve_bfield_all.dat" ] && \
        nb=$(wc -l < "$DATA_U/coil${N}/full_assembly_sleeve_bfield_all.dat")
    echo "[coil${N}] done (db=${db}, bfield_all lines=${nb})"
done
