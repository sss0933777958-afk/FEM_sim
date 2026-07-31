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

def sub(t, old, new, label):
    """[ADDED] 替換必須命中，否則 abort —— 比照 long2016 run_graded_all.sh 的 sed 防呆。
       原版無此檢查，導致下面 output-dir 那條靜默未命中、六顆 coil 全倒進 data/gap_200um/coil1/。"""
    if old not in t:
        sys.exit("ABORT [%s]: 替換字串未命中: %r" % (label, old))
    return t.replace(old, new)

# CWD -> db/sim/<tag>
t=sub(t, r"\db\sim\gap_200um_solve", r"\db\sim\%s"%tag, "CWD")
# RESUME -> my variant mesh
t=sub(t, r"db\geom\mesh_graded_basegap", r"db\mesh\mesh_graded_%s"%tag, "RESUME-mesh")
# energise coil N (canonical energises coil1 at 'R,1,1,TURNS*1')
#   號誌維持 TURNS*1（正）→ raw FEM 六極全 sink，與 R300/R700 一致；
#   翻號由 MATLAB config s_source=[-1×6] 全域處理，deck 端不可自行改號誌。
t=sub(t, "R,1,1,TURNS*1,COIL_THK,COIL_H", "R,%s,1,TURNS*1,COIL_THK,COIL_H"%N, "excite-coil")
# [FIXED] output dir: 真實路徑是 ...\ANSYS_data\hung_hexapole\data\gap_200um\coil1，
#   原版找 r"\data\hung_hexapole\data\gap_200um\coil1" —— "ANSYS_data" 的 data 前是底線非反斜線
#   → 永不命中。改用無歧義的短片段。
t=sub(t, r"\data\gap_200um\coil1", r"\data\%s\coil%s"%(tag,N), "output-dir")
t=sub(t, "data/gap_200um/coil1/", "data/%s/coil%s/"%(tag,N), "output-dir-comment")
for f in ("coord_all","bfield_all","coord_wp","bfield_wp"):
    t=sub(t, "'coil1_%s'"%f, "'coil%s_%s'"%(N,f), "outfile-"+f)
t=t.replace("coil 1", "coil %s"%N)
# SAVE solved db before exit (deliver db/sim)
t=sub(t, "/EXIT,NOSAVE", "SAVE,'sim_coil%s','db'\n/EXIT,NOSAVE"%N, "SAVE")

# [ADDED] 殘留佔位符檢查：兩個 canonical tag 都不該再出現
for bad in ("gap_200um", "mesh_graded_basegap"):
    if bad in t:
        sys.exit("ABORT: 產生的 deck 仍殘留 %r（替換不完整）" % bad)
open(out,'w',encoding='utf-8').write(t)
PY
    [ -f "$tmp" ] || { echo "[$tag coil$N] deck 生成失敗，中止"; exit 1; }
    echo "[$tag coil$N] solving ..."
    ( cd "$simcwd_u" && "$MAPDL" -b -np 4 -m 24000 -dir "$simcwd_w" -j "sim_coil${N}" -i "$tmp" -o "$simcwd_u/solve_coil${N}.out" >/dev/null 2>&1 )
    # half-clean this coil (keep .db + 主 .rmg + solve log; drop worker/esav/etc)
    # [MODIFIED] 補刪 per-worker .out/.log（db-folder-retention L32 明列該刪；原版漏掉，
    #   造成 db/sim/R700、R300 各殘留 36 個 sim_coilN_{0..3}.{err,out,log}）
    find "$simcwd_u" -maxdepth 1 -type f \( -name "sim_coil${N}_*.rmg" -o -name "sim_coil${N}*.esav" -o -name "sim_coil${N}*.full" -o -name "sim_coil${N}*.DSP*" -o -name "sim_coil${N}*.page*" -o -name "sim_coil${N}*.err" -o -name "sim_coil${N}*.stat" -o -name "sim_coil${N}*.ldhi" -o -name "sim_coil${N}_[0-9].out" -o -name "sim_coil${N}_[0-9].log" \) -delete 2>/dev/null
    rm -rf "$simcwd_u/scratch" 2>/dev/null
    rm -f "$tmp"
    # [MODIFIED] 驗收：4 個 .dat 齊全 + 行數 + 真正的 *** ERROR（原本 grep "ERROR" 會誤中一般字樣）
    err=$(grep -c '\*\*\* ERROR' "$simcwd_u/solve_coil${N}.out" 2>/dev/null)
    ndat=$(ls "$datadir/coil${N}"/coil${N}_{coord,bfield}_{all,wp}.dat 2>/dev/null | wc -l)
    nrow=$(wc -l < "$datadir/coil${N}/coil${N}_bfield_all.dat" 2>/dev/null || echo 0)
    echo "[$tag coil$N] done (db=$([ -f "$simcwd_u/sim_coil${N}.db" ] && echo Y || echo N), dat=$ndat/4, bfield_all=$nrow 行, *** ERROR=$err)"
    if [ "$ndat" -ne 4 ] || [ "$err" -ne 0 ]; then
        echo "[$tag coil$N] ⚠ 驗收未過 —— 中止（檢查 $simcwd_u/solve_coil${N}.out）"; exit 1
    fi
  done
done
