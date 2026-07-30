#!/usr/bin/env bash
# run_mesh_tip_variant.sh — 由 canonical mesh deck MT_Mesh_Graded.txt 產「tip-radius 變體」graded mesh。
#   用法： bash run_mesh_tip_variant.sh <tag> <POLE_TIP_R值> [WP_REF_LVL]
#     例： bash run_mesh_tip_variant.sh tip400um 400.0e-6      # CNC 400 um 鈍尖（WP_REF_LVL 預設 1）
#          bash run_mesh_tip_variant.sh tip20um  20.0e-6       # EDM 20 um 尖尖
#          bash run_mesh_tip_variant.sh tip400um 400.0e-6 2    # WP 內球 refine +2（更密、solve 更重）
#   WP_REF_LVL：鈍尖無奇異點 → smrt5 把 WP 鋪很粗（R<=150um ~11 節點 < 18 參數）；+1 補到 ~78（>18）。
#     局部 WP-core refine 在此幾何一定 fail → 只能整顆內球 uniform refine（deck 內 *if 控制）。
#   做法：把 canonical deck 複製到 scratch，literal-replace 5 個變體旋鈕（POLE_TIP_R + 輸出路徑/檔名
#         + EREFINE 內圈界放寬），再 MAPDL -b 建幾何(scenario-A 錐固定尖端後退)+graded mesh+SAVE。
#   ★ canonical deck 永久保持 baseline（POLE_TIP_R=40e-6、scenario-A 參數在 R=40e-6 時塌回 baseline）；
#     本 driver 只動 scratch copy，絕不改 canonical。
#   EREFINE 放寬：400um 鈍尖後退到 dwp~1.9mm 卡在原 2e-3 內圈界外 → 放寬 1.5e-3 讓鈍尖進加密帶
#     （ESEL,S,MAT,,MAT_MT 先選 → 只多拉鋼、不拉空氣）；Z 界 -11.5→-12.0 使上極後退尖仍排除於 Step2。
set -u
MAPDL="G:/ANSYS Inc/v252/ansys/bin/winx64/MAPDL.exe"
ROOT="G:/my_workspace/code/FEM_sim/magnetic_sim/ANSYS/main"
DECK="$ROOT/apdl/long2016_hexapole_halfcut/mesh/MT_Mesh_Graded.txt"

TAG="${1:?用法: bash run_mesh_tip_variant.sh <tag> <POLE_TIP_R值> [WP_REF_LVL]  例 tip400um 400.0e-6}"
RVAL="${2:?缺 POLE_TIP_R 值（例 400.0e-6）}"
WPLVL="${3:-1}"
MESHCWD="$ROOT/ANSYS_data/long2016_hexapole_halfcut/db/mesh/$TAG"
JOB="mesh_$TAG"
mkdir -p "$MESHCWD"
tmp="$MESHCWD/_deck_$TAG.txt"

# --- literal-replace 5 個變體旋鈕（python 免 sed backslash 地獄）---
DECK="$DECK" TMP="$tmp" TAG="$TAG" RVAL="$RVAL" WPLVL="$WPLVL" python -c '
import os
src=open(os.environ["DECK"],encoding="utf-8",errors="replace").read()
tag=os.environ["TAG"]; rval=os.environ["RVAL"]; wplvl=os.environ["WPLVL"]
reps=[
 ("POLE_TIP_R    = 40.0e-6",            "POLE_TIP_R    = "+rval),          # 1) 變體尖端倒圓半徑
 ("db\\geom\\mesh_graded",              "db\\mesh\\"+tag),                 # 2) /CWD 輸出夾
 ("/FILNAME,mesh_graded,1",             "/FILNAME,mesh_"+tag+",1"),        # 3) jobname/存檔名
 ("CENT,X, 2e-3, 18e-3",                "CENT,X, 1.5e-3, 18e-3"),          # 4) EREFINE 內圈放寬（兩處）
 ("CENT,Z, -1.0, -11.5e-3",             "CENT,Z, -1.0, -12.0e-3"),         # 5) Step2 下極 Z 界
 ("WP_REF_LVL = 0",                     "WP_REF_LVL = "+wplvl),            # 6) WP 內球 uniform refine level
]
cnt={}
for a,b in reps:
    cnt[a]=src.count(a); src=src.replace(a,b)
open(os.environ["TMP"],"w",encoding="utf-8").write(src)
# 命中報告（給 shell 檢查）
print("HIT POLE_TIP_R="+str(cnt["POLE_TIP_R    = 40.0e-6"]))
print("HIT CWD="+str(cnt["db\\geom\\mesh_graded"]))
print("HIT FILNAME="+str(cnt["/FILNAME,mesh_graded,1"]))
print("HIT CENTX="+str(cnt["CENT,X, 2e-3, 18e-3"]))
print("HIT CENTZ="+str(cnt["CENT,Z, -1.0, -11.5e-3"]))
print("HIT WPREF="+str(cnt["WP_REF_LVL = 0"]))
' | tee "$MESHCWD/_subst_report.txt"

# --- 命中檢查（防呆）：每個旋鈕都要有替換到，殘留 mesh_graded 表示漏換 ---
ok=1
grep -q "^HIT POLE_TIP_R=1$" "$MESHCWD/_subst_report.txt" || { echo "[ABORT] POLE_TIP_R 未命中(需=1)"; ok=0; }
grep -q "^HIT CWD=1$"        "$MESHCWD/_subst_report.txt" || { echo "[ABORT] /CWD 路徑未命中(需=1)"; ok=0; }
grep -q "^HIT FILNAME=1$"    "$MESHCWD/_subst_report.txt" || { echo "[ABORT] /FILNAME 未命中(需=1)"; ok=0; }
grep -q "^HIT CENTX=2$"      "$MESHCWD/_subst_report.txt" || { echo "[ABORT] EREFINE CENT,X 需命中2處"; ok=0; }
grep -q "^HIT CENTZ=1$"      "$MESHCWD/_subst_report.txt" || { echo "[ABORT] EREFINE CENT,Z 需命中1處"; ok=0; }
grep -q "^HIT WPREF=1$"      "$MESHCWD/_subst_report.txt" || { echo "[ABORT] WP_REF_LVL 未命中(需=1)"; ok=0; }
# 註：deck 註解行仍會提及 mesh_graded.db（無害），故不做 blanket 殘留檢查；上面 HIT 計數才是權威。
if [ "$ok" != "1" ]; then rm -f "$tmp"; exit 1; fi
grep -n "POLE_TIP_R    =" "$tmp" | head -1
grep -n "/CWD" "$tmp" | head -1
grep -n "/FILNAME" "$tmp" | head -1

# --- 建幾何 + mesh + SAVE ---
echo "[$TAG] meshing (POLE_TIP_R=$RVAL) ..."
"$MAPDL" -b -np 4 -m 48000 -dir "$MESHCWD" -j "$JOB" -i "$tmp" -o "$MESHCWD/$JOB.out"

# --- 驗證：0 error + MESH SUMMARY + .db 落點 ---
nerr=$(grep -c "\*\*\* ERROR" "$MESHCWD/$JOB.out" 2>/dev/null || echo "?")
db=$( [ -f "$MESHCWD/$JOB.db" ] && echo Y || echo N )
echo "[$TAG] done | ERROR=$nerr | $JOB.db=$db"
echo "----- MESH SUMMARY -----"
grep -A6 "MESH SUMMARY" "$MESHCWD/$JOB.out" 2>/dev/null | head -8

# --- half-clean：保 $JOB.db + $JOB.out，清中間檔 + worker log + scratch deck ---
rm -f "$MESHCWD/$JOB"*.esav "$MESHCWD/$JOB"*.full "$MESHCWD/$JOB"*.page* \
      "$MESHCWD/$JOB"*.stat "$MESHCWD/$JOB"*.err "$MESHCWD/$JOB"*.ldhi "$MESHCWD/$JOB"*.mlv \
      "$MESHCWD/$JOB"*.BCS "$MESHCWD/$JOB"*.DSP* "$MESHCWD/$JOB"*.lock "$MESHCWD/$JOB"*.mntr \
      "$MESHCWD/$JOB"[0-9]*.out "$MESHCWD/$JOB"[0-9]*.err "$MESHCWD/$JOB"[0-9]*.log \
      "$MESHCWD/$JOB"*.rmg "$tmp" "$MESHCWD/_subst_report.txt" 2>/dev/null
echo "[$TAG] 保留：$JOB.db + $JOB.out（中間檔已清）"
