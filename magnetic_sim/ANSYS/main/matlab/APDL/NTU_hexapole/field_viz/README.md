# matlab/NTU_hexapole/field_viz/ — 場視覺化 + sensor 平均

NTU FEM 場的視覺化與 sensor 平均（讀 `ANSYS_data/NTU_hexapole/data/{singlepole,upper_assembly}/coil1/*.dat`，真實 FEM 節點）。
`data/` = 本組 `.mat` 成果（per matlab-output-layout）；`results/` = 只留 `.pdf`。

## `code/plot/plot_sensor_Bcircuit.m` → `figures/sensor_Bcircuit_xz.png`
sensor 附近磁路箭頭圖（**x–z @ y=−0.99**，一張含 S1+S2），比照長飛 `plot_P2sensor_Braw_P1exc.m`：
- per-arrow **|B| log 色階(turbo，mT)**、箭頭=單位方向；y 薄片 + 視窗盒選真實節點 + 格點抽稀（留 |B| 最大真實節點、非內插）。
- 座標 **mm**、|B| 色條 **mT**（unit-reference）。薄板截面**只黑框輪廓、不填色**（磁極內磁路箭頭看得見）。
- S1/S2 綠 disc + 紅 n+；翻電流後 → 板內 in-plane 磁通朝 **−x（極尖）＝往磁極跑**。
- **gap-fill 重心內插**：薄板板內（穿厚只 1 層 tet、無真實節點）+ 板↔sensor 球間粗空氣（節點稀）都是空洞 → 對「附近 0.10mm 內無真實顯示節點」的格點，用單一局部 3D Delaunay(全節點) pointLocation+重心內插補箭頭（板內落鋼件 tet、空氣空洞落跨界 tet 過渡場）。板外密處仍真實節點、不重複。圖上紅字標明「gap-fill ... barycentric interpolation」（依 plot-real-nodes 規則）。

## `code/scripts/sensor_average.m`（單極）→ `results/sensor_average_singlepole.pdf`
單極 sensor 有效範圍（取樣圓柱 R=0.15mm、沿 n+ 高 H=0.1mm）內 **1000 點面積平均（重心內插）**：
- 面積均勻取樣（`r=R√rand`）+ 局部真實節點 `delaunayTriangulation`→`pointLocation`+重心內插（=FEM 線性 shape fn，非 scatteredInterpolant）。
- 報 ⟨B·n+⟩、⟨|B|⟩（**mT**）每 sensor + 1000 點 std；xelatex 出 PDF（`results/` 只留 .pdf）。**重用** `import_ansys_data.m`。

### 單極結果（1A×50T，翻電流後；1000 點內插）
| sensor | ⟨B·n+⟩ | ⟨|B|⟩ | 加密球節點 |
|---|---|---|---|
| S1(n+ −z) | +2.1457 mT (±0.017) | 2.2711 mT | 1421 |
| S2(n+ +z) | +1.8533 mT (±0.026) | 2.0409 mT | 1419 |

⚠ 兩加密球 FEM 節點 1421/1419（先 mesh 兩球同 ESIZE、free-tet 對嵌入球無法逐顆完全相等，差 2 顆 0.14%）；下游 1000 點內插平均兩顆等量、不受影響。B·n+ 兩顆皆正＝場沿 n+（離板）；板內磁通 −x 朝極尖。

## `code/scripts/sensor_realnode_average.m`（獨立）→ `results/sensor_realnode_average.pdf` + `data/sensor_realnode_average.mat`
**真實節點平均**（不內插；per plot-real-nodes）：直接對 sensor 加密球的**實際 FEM 節點**取無權重 mean+std。獨立腳本、不動上面 1000-點 `sensor_average.m`。
- 讀 `import_ansys_data(dir,'sensor_S1'|'S2', 'upper_assembly'|'singlepole')`（sim deck 已抽的球全節點）。⟨|B|⟩=mean(bsum)、⟨B·n+⟩=mean(sgn·bz)（S1 −z 板下 / S2 +z 板上，正=離板），mT。
- 出 PDF（upper vs 單極 baseline 同法 + 倍率）+ 存 `.mat` 到 `data/`。

### 上層完整組合結果（upper_assembly，1A×50T post0；真實節點平均）
| sensor | ⟨\|B\|⟩ upper | ⟨B·n+⟩ upper | ⟨\|B\|⟩ 單極 baseline | 倍率 |
|---|---|---|---|---|
| S1（板下 z−0.41, n+ −z） | **5.3860 mT** (±0.038) | +5.1090 mT (±0.040) | 2.2749 mT | **2.37×** |
| S2（板上 z+0.66, n+ +z） | **4.3160 mT** (±0.050) | +3.8659 mT (±0.058) | 2.0363 mT | **2.12×** |

yoke + 6 導柱 + 鄰極低磁阻 return path → pole0 板磁通大增 → sensor 場 ~2.1–2.4× 於單極。B·n+ 皆正（沿 n+ 離板）。單極 baseline 此處用**同真實節點法**（2.2749/2.0363，與 1000-點內插 2.271/2.041 相符）。
