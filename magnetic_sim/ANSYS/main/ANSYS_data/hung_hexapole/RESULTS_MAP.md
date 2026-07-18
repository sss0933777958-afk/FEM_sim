# hung_hexapole — RESULTS_MAP（讀 FEM 結果前必查，凌駕 memory）

依 `.claude/rules/result-read-safety.md`：讀任何 hung_hexapole 結果前，先在此核路徑 + 指紋。

> **⚠ 2026-07-14 資料夾結構已翻轉：`data/coil<N>/<variant>/` → `data/<variant>/coil<N>/`**（使用者拍板，全包統一）。
> dir 現為 `data/<variant>/coil<N>/`（例 `data/gap_200um/coil1`、`data/no_gap/coil6`）。下表「dir」欄仍寫舊式，
> 實際請讀成 `<variant>/coil<N>`；**指紋（|B| / 節點數）不受改名影響**。
> 讀取碼一律 `import_ansys_data(fullfile(root,<variant>,'coil<N>'), <dataset>, 'coil<N>')`。

## 模型
hung_hexapole 6 極 electromagnetic microprobe，**native APDL primitive 重建**（非 CAD 匯入；匯入不能布林）。WP 在原點。詳見 memory `[[hung-hexapole-primitive-sim]]`。

## 網格（db/）
| 路徑 | 是什麼 | 指紋 |
|---|---|---|
| `db/mesh_graded/mesh_graded.db` | **canonical 網格**：SMRTSIZE5 + 二級 EREFINE（step1 全極 cone、step2 下極）+ 6 SOURC36 線圈（N70, r[4,7] 繞導柱） | 560,277 節點 / 3.35M 元素 / 7 鋼件 / 6 線圈 |
| `db/from_parasolid/hung_hexapole_full.db` | CAD 匯入真實幾何（97 vol，**非 primitive、非 mesh**，參考用） | 無 .dat，不是結果 |

## 結果（data/）＝ gap200µm μ_eff 逐極激發
**唯一 case**：`gap_200um`（支撐座等效導磁：**下極 μ=52.02 / 上極 μ=59.364**；1A 激發、N70）。μ_eff 只套 **6 整塊支撐座**（導柱/yoke/磁極排除；每極 EMODIF 下 10792 / 上 17349 元素）。**2026-07-07 重解**（db 清除 stray MAT99 導柱污染 + 加導柱排除 + μ_eff 改值）。

| dir（新式 `<variant>/coil<N>`） | 極（方位/層） | WP 中心 \|B\| | 檔案（每夾 4 個 .dat） |
|---|---|---|---|
| `data/gap_200um/coil1/` | P1（0°/下） | **5.88 mT** | coil1_{coord,bfield}_{all,wp}.dat |
| `data/gap_200um/coil2/` | P2（180°/上） | **7.53 mT** | coil2_… |
| `data/gap_200um/coil3/` | P3（120°/下） | **5.89 mT** | coil3_… |
| `data/gap_200um/coil4/` | P4（300°/上） | **7.50 mT** | coil4_… |
| `data/gap_200um/coil5/` | P5（60°/上） | **7.50 mT** | coil5_… |
| `data/gap_200um/coil6/` | P6（240°/下） | **5.95 mT** | coil6_… |

> **另有 `no_gap` 變體**（同 6 極、同網格、無氣隙對照）：`data/no_gap/coil<N>/`。

**指紋/對稱**：下極 P1/P3/P6 |B|≈5.9mT、Bz≈−3.4mT（下傾）；上極 P2/P4/P5 |B|≈7.5mT、Bz≈+3.45mT（上傾）；WP 徑向分量沿各極方位（3-fold 對稱）。（舊值 下7.3/上8.17mT 為 μ=88.41/101.36；新 μ 較低→|B| 較低。）
- `*_all.dat`：全 560k 節點 B（Tesla；PRNSOL 大值 >100mT 印製受 /FORMAT 上限，WP 區 <10mT 不受影響）。
- `*_wp.dat`：WP 球座標 R<2mm 區（~400k 節點，含匯聚細網格）。

## 座標/單位
公尺 MKS；B 輸出 Tesla；WP 在原點；下極 z<0（傾 −5.7°）、上極 z>0（傾 +34.7°）。

## 讀取注意
- 變體：`gap_200um`（μ_eff 氣隙）+ `no_gap`（無氣隙對照）；路徑 `data/<variant>/coil<N>/`（無 baseline）。
- 電流 = 1A（FEM 激發，見 `.claude/rules/fit-current-matches-sim.md`）；下游 fit 的 I 必須 = 1。
- source 慣例（電荷模型）：下極 P1/P3/P6 raw = sink，需翻號成 source（見 `.claude/rules/charge-model-source-convention.md`）。
