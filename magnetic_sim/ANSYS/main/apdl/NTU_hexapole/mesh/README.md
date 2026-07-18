# apdl/NTU_hexapole/mesh/ — 網格 deck

NTU_hexapole 的 mesh-only deck（建幾何 + 網格 + SAVE，不 solve）。三支：單極 `MT_Mesh_Graded.txt`、上層完整 `MT_Mesh_Graded_Upper.txt`、**全組合 `MT_Mesh_Graded_Full.txt`（6 極 full_assembly）**。

## `MT_Mesh_Graded_Full.txt`（full_assembly：6 磁極 + 6 導柱 + yoke）★保留 5µm 圓角
完整 6 極（世界座標＝總組合 frame）。**保留尖端 5µm 圓角（不簡化直線）** + 含空氣全模型。

- **圓角保留法（mm 建 → VLSCALE 縮公尺）**：5µm 圓角公尺尺度 `LARC` 建不起（絕對 5e-6 near-tol），mm 尺度建得起。故 steel 幾何用 **mm-magnitude 建（SCB=1.0，圓角 LARC OK）→ `VLSCALE,ALL,,,1e-3,1e-3,1e-3,,0,1` 對原點縮成公尺**（只縮 keypoint、不重跑 LARC）→ 圓角變 5e-6m 小弧、單位正確。mesh 時圓角弧 1 元素 + ANSYS 局部漸變（已驗證成功、非簡化）。
- **幾何**：下極 z[8.8,9.05]@0/120/240°、上極 z[9.55,9.8]@60/180/300°（板厚0.25、9 段 profile 含圓角）+ 6 導柱 full r2.5（下極柱 z[9.05,26.8]、上極柱 z[9.80,26.8]）+ yoke（rounded-rect r5 外盤 − cutout r5、**無孔**、z[20.8,26.3]）。**VADD 導柱+yoke=連通 return body、6 薄板獨立**。WP=(0,0,9.3mm)。
- **空氣三層**（公尺，中心 WP）：內球 R7 + 中球 R48 + 外圓柱 R100×H110（內球+外圓柱概念照龍飛；中球/R100×H110 沿用 upper）。
- **網格**：6 薄板 `ESIZE 0.1mm`（≥2 層穿厚）先 mesh → return body+空氣 `SMRTSIZE,5` → 尖端近 WP `EREFINE +1`。
- **產物**：`ANSYS_data/NTU_hexapole/db/**mesh**/full_assembly/full_assembly.db`（**新 db/mesh/ 層、不覆蓋 db/geom**）。**1,787,195 節點 / 10,754,865 元素**（STEEL 2.08M / AIR1 999k / MID 7.67M / AIR2 1911；STEEL 13 vol、體積 1.5795e-5 m³ 與幾何吻合）。`-np 4 -m 32000`。**mesh-only（無 /SOLU/coil/BC）**。
- **檢視**：`../gui/MT_Render_SteelMesh.txt`（batch PNG，可靠）→ `figures/NTU_hexapole/full_assembly_mesh/steel_mesh_{iso,front_xz,top_xy,pole_plates,tip_fillet_zone}.png`。或 `MT_View_SteelMesh.txt`（互動 GUI；4.37GB db 需 `-m 24000`+單一乾淨啟動）。

## `MT_Mesh_Graded.txt`（單磁極薄板 + 短導柱）
比照 `long2016_hexapole_halfcut/mesh/MT_Mesh_Graded.txt` 的 region-graded 做法（純 SMRTSIZE,5 + EREFINE，不用硬 ESIZE）。**建幾何 + 網格 + SAVE，不 solve、不加 BC/coil**。

- **幾何（metre/MKS, SC=1e-3）**：磁極薄板（pole.STEP 9 段 profile + VEXT 0.25mm；⚠ 5µm 尖端 fillet 在 metre scale 無法建 → 改兩段直線、尖端位置不變）+ 短導柱（實心階梯圓柱）VADD 成鋼件。沿用 `../geom/export/MT_Geom_PoleAssembly.txt` 的 profile。
- **空氣三層**（中心=原點≈極尖，比照龍飛外邊界）：內球 7mm + **中球 32mm**（須 >鋼件最遠 ~29mm 且 <外圓柱半高 35mm，避免與端蓋相切→VOVLAP 退化丟外圓柱）+ 外圓柱 R80×H70mm。
- **Sensor 加密球**：2 顆 R0.15mm air 球 @ (13.15,−0.99,−0.41)/(13.15,−0.99,+0.66) mm。**先 mesh 兩球（同 ESIZE 0.025mm→兩球拓樸幾乎相同）再 mesh 鋼件+空氣（conform）**，求兩球節點數相等；實得 S1=1421/S2=1419（free-tet 對嵌入球無法逐顆完全相等，差 2 顆 0.14%，使用者接受）。⚠ n+ S1→−z / S2→+z 供下游 B·n+ 符號；不建 disc 實體。
- **材料**：MAT_MT=2(μr=280 鋼)、AIR1/MID/AIR2/SENS(μr=1)；VOVLAP 後 by location+size 指派（不硬編 vol 號）。
- **網格**：純 smrt5（薄板本質 → 鋼件 ~0.12mm/2 層、~1.37M 顆；試過 AESIZE 封頂會觸發 gradient blowup，故回純 smrt5）+ 極尖 cone EREFINE +1 + sensor 球 +1。
- **產物**：`ANSYS_data/NTU_hexapole/db/mesh_graded/mesh_graded.db`（~279k 節點 / ~1.67M 元素；-np 1，EREFINE 避 DMP 崩）。

## `MT_Mesh_Graded_Upper.txt`（上層完整組合：3 磁極 + 6 導柱 + yoke）
單極版擴成完整上層（供 `../sim/upper_assembly/` solve）。沿用 `../geom/export/MT_Geom_UpperAssembly.txt` 尺寸、metre/MKS。

- **幾何**：3 極 WPROTA 0/120/240（tip 5µm fillet→直線）+ 6 導柱 full r2.5 + yoke（外盤−中心 cutout，**無導柱孔**）。
  - **鋼件融合**：只 `VADD` 導柱+yoke 成 return-path body；**3 薄板保持獨立 volume**。導柱穿過 solid yoke（overlap→VADD 融），極-導柱在 z0.25 面接觸交 VOVLAP conform（已驗 pole0-post0 界面 0 重合節點＝無裂縫、單一連通鋼件）。
- **空氣三層**（enlarged 自單極以含 yoke，yoke 角落距原點 ~43.6mm）：內球 7 / **中球 48**（>43.6 全包、<外圓柱半高 55 避相切）/ 外圓柱 **R100×H110（z±55）**。
- **⚠ 薄板網格關鍵坑**：`VADD` 全鋼件（含板）→ 板被併進大 yoke body、smrt5 全域變粗 → 板只 **1 層 tet slivers**（穿厚僅 z0/z0.25 兩層節點）；**EREFINE 也救不了**（conforming-air 夾住的單層板，13k 元素只加 2k）。
  **解**：對「薄板鋼件（`VSEL` MAT_MT + **LOC**（非 CENT！VSEL 用 LOC）z<0.3mm）」下 **ESIZE 0.1mm 獨立先 mesh**（0.25mm 板→**3 層穿厚**，優於單極 baseline 2 層）；return body（導柱+yoke）+ 空氣後續 smrt5 自動 conform/grade 到細板面。+ 極尖 EREFINE +1。
- **產物**：`ANSYS_data/NTU_hexapole/db/mesh_graded_upper/mesh_graded_upper.db`（**710k 節點 / 4.27M 元素**；STEEL 1.05M / AIR1 340k / MID 2.87M / AIR2 1989 / SENS 13.6k；sensor 球 1421/1419；-np 1 -m 40000，~9min）。
- **跑法**：同下（`-j mesh_graded_upper -m 40000`）。

## 跑法
```bash
cd ANSYS_data/NTU_hexapole/db/mesh_graded
"G:/ANSYS Inc/v252/ansys/bin/winx64/MAPDL.exe" -b -np 1 -m 24000 \
  -dir <該夾> -j mesh_graded -i apdl/NTU_hexapole/mesh/MT_Mesh_Graded.txt -o run.out
```

## 檢查鋼件（GUI）
複製 `mesh_graded.db`→`view.db`（避免 GUI RESUME 蓋 canonical），用 `../gui/MT_View_SteelMesh.txt`：
```bash
MAPDL.exe -g -np 1 -dir <db夾> -j view -i apdl/NTU_hexapole/gui/MT_View_SteelMesh.txt
```
GUI 內 `ESEL,S,MAT,,2 $ EPLOT`＝鋼件；`ESEL,S,MAT,,5 $ EPLOT`＝sensor 加密球。

## 後續（未做）
solve deck（RESUME mesh_graded.db → 加 coil + D,ALL,MAG,0 + magsolv）、sensor 場取樣 postproc（1000 點面積平均）。
