#!/usr/bin/env bash
# build_ws_variants.sh — 產生 Hung hexapole 不同「工作空間半徑」(R_sphere) 的幾何變體
# ---------------------------------------------------------------------------
# 用途：從 canonical 幾何 deck `../export/MT_Geom_FullAssembly.txt`（R_sphere=0.5mm）
#   sed 出指定 R 的暫存 deck → 跑 MAPDL 建幾何 → 存 geom .db + IGES → IGES 修 unit flag
#   (6→1, 因 deck 用英寸 MM=1/25.4 建、SolidWorks ×25.4 讀回正確 mm) → 複製到 model_check/。
#   canonical deck 不改；每個 R 產物落各自 db/geom/geom_full_R<tag>/ + model_check/。
#
# R_sphere = tip 到 WP 球心距離（工作空間半徑）；6 tip 依 magic-angle 線性隨 R 縮放，
#   外結構(yoke/blocks/coils/pole 長度+傾角)固定。tip-to-tip = 2R。
#
# 用法：  bash build_ws_variants.sh 0.7 0.3        # 產 R700、R300
#         bash build_ws_variants.sh 0.5            # 重建 baseline
# ---------------------------------------------------------------------------
set -u
MAPDL="G:\\ANSYS Inc\\v252\\ansys\\bin\\winx64\\MAPDL.exe"
ROOT="G:/my_workspace/code/FEM_sim/magnetic_sim/ANSYS/main"
DECK="$ROOT/apdl/hung_hexapole/geom/export/MT_Geom_FullAssembly.txt"
GEOM="$ROOT/ANSYS_data/hung_hexapole/db/geom"
MC="$ROOT/model_check/hung_hexapole"
EXP="$ROOT/apdl/hung_hexapole/geom/export"

for R in "$@"; do
    tag="R$(awk -v r="$R" 'BEGIN{printf "%d", r*1000}')"      # 0.7 -> R700
    cwd_u="$GEOM/geom_full_$tag"; mkdir -p "$cwd_u" "$MC"
    cwd_w="$(echo "$cwd_u" | sed 's|/g/|G:\\|; s|/|\\|g')"
    tmp="$EXP/_ws_$tag.txt"
    sed -e "s|R_sphere = 0.5\\*MM|R_sphere = ${R}*MM|" \
        -e "s|geom_full'|geom_full_${tag}'|" \
        -e "s|!\\*\\*\\*\\* SAVE, 'Full_Assembly_filleted', 'db'|SAVE, 'Full_Assembly_${tag}', 'db'|" \
        -e "s|IGESOUT, 'hung_full_assembly', 'iges'|IGESOUT, 'hung_full_assembly_${tag}', 'iges'|" \
        "$DECK" > "$tmp"
    echo "[$tag] R_sphere=${R}mm  building geometry ..."
    ( cd "$cwd_u" && "$MAPDL" -b -np 2 -m 8000 -dir "$cwd_w" -j "hung_ws_$tag" -i "$tmp" -o "$cwd_u/build.out" >/dev/null 2>&1 )
    # ANSYS IGES(db/geom 中繼) unit flag 6->1，OCC 才讀回正確 mm
    sed -i "s/,1.0,6,,/,1.0,1,,/" "$cwd_u/hung_full_assembly_${tag}.iges"
    # 交付檢查 = STEP (deliver-step-for-check)：OCC IGES->mm STEP -> model_check/，並讀回驗 tip=R
    python "$EXP/../scripts/iges_to_step.py" "$cwd_u/hung_full_assembly_${tag}.iges" "$MC/hung_full_assembly_${tag}.step" "$R"
    # keep only .db + .iges(中繼) + build.out in the geom folder (db-folder-retention)
    find "$cwd_u" -type f ! -name "*.db" ! -name "*.iges" ! -name "build.out" -delete 2>/dev/null
    rm -f "$tmp"
    echo "[$tag] done -> $cwd_u/Full_Assembly_${tag}.db  +  $MC/hung_full_assembly_${tag}.step"
done
