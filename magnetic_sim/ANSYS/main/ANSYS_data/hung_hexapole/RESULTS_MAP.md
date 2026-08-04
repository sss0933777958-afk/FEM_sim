# hung_hexapole — RESULTS_MAP（讀 FEM 結果前必查，凌駕 memory）

依 `.claude/rules/result-read-safety.md`：讀任何 hung_hexapole 結果前，先在此核路徑 + 指紋。

**路徑格式**：`data/<variant>/coil<N>/`（2026-07-14 起，`data/coil<N>/<variant>/` 舊式已淘汰）。
讀取碼一律 `import_ansys_data(fullfile(root,<variant>,'coil<N>'), <dataset>, 'coil<N>')`。

## 模型
hung_hexapole 6 極 electromagnetic microprobe，**native APDL primitive 重建**（非 CAD 匯入；匯入不能布林）。
WP 在原點（`SPH_OFST = 0`）。詳見 memory `[[hung-hexapole-primitive-sim]]`。

- **coil 編號 = identity**：mesh deck `fc_ra(1..6) = 0/180/120/300/60/240°` + `REAL,fc_pi` → real set k = 論文 P k。
- **號誌**：deck `R,N,1,TURNS*1`（正號）→ raw 六極**全 sink** → config `s_source = [-1×6]` 全域翻號成 source。
- 電流 = 1A（FEM 激發）；下游 fit 的 I 必須 = 1（`.claude/rules/fit-current-matches-sim.md`）。
- 單位：公尺 MKS，B 輸出 Tesla。

---

## ✅ 可用結果：R 家族（原始設計 + 兩個 R_sphere 變體）

三變體**同一支純版 deck**（`MT_Mesh_Graded.txt`）只改 `R_sphere`；`R500` = 原始設計值、
= `config/hung_hexapole/mt_constants.m` 的 `R_norm = 500e-6`，**現為 `default_variant`**。

| variant | R_sphere（tip→WP） | mesh 節點 / 元素 | `.dat` bfield_all 行數 | wp matched 節點 |
|---|---|---|---|---|
| `R300` | 0.3 mm | 589,264 / 3,520,338 | 589,248 | 439,015 |
| **`R500`**（default）| **0.5 mm** | **560,277 / 3,349,759** | **560,261** | **402,771** |
| `R700` | 0.7 mm | 791,452 / 4,745,571 | 791,436 | 622,074 |

⚠ **節點數對 R 非單調**（R300 > R500）：graded `EREFINE` 是 band-based，R 縮小時極尖端擠進內層細網格球
反而觸發更多細化。別拿「幾何變大→節點變多」當驗收判準。

### 指紋：WP 中心（最近節點，r_min ≈ 0.6–1.2 µm）1A 激發

| variant | 下極 P1 / P3 / P6 \|B\| | 下極 Bz | 上極 P2 / P4 / P5 \|B\| | 上極 Bz | \|B\|_max（wp 集） |
|---|---|---|---|---|---|
| `R300` | 16.58 / 16.70 / 16.61 mT | −9.1 mT | 15.75 / 15.56 / 15.59 mT | +8.4 mT | 1.31–1.41 T |
| **`R500`** | **9.72 / 9.73 / 9.79 mT** | **−5.33 mT** | **9.24 / 9.16 / 9.17 mT** | **+4.93 mT** | 1.19–1.27 T |
| `R700` | 6.96 / 6.95 / 6.96 mT | −3.81 mT | 6.58 / 6.52 / 6.53 mT | +3.53 mT | 1.13–1.22 T |

**對稱性判準**：同層三顆 \|B\| 散佈 < 1%；下極 Bz < 0、上極 Bz > 0；下極略強於上極（+5~6%）。
任一條不符 → 停下來查是否讀錯 variant。

### 校正結果（single / USE_BIAS=false、R_select=150 µm）
`R500`：ℓ̂ = 887.7 µm、ᴮĝ_I = 9.2948 mT/A、RMSPE 2.58%、𝒞 = 2481、κ = 0.8652，K̄_I 四項判準全過。
PDF：`…/Calibration_using_FEM_modeling/results/hung_hexapole/single/model_results_current_R500.pdf`。

### 網格 / 解 db
| 路徑 | 內容 |
|---|---|
| `db/mesh/mesh_graded_R{300,500,700}/mesh_graded.db` | canonical 網格（SMRTSIZE5 + 二級 EREFINE：step1 全極 cone、step2 下極）+ 6 SOURC36 線圈（N70, r[4,7] 繞導柱） |
| `db/sim/R{300,500,700}/sim_coil<N>.{db,rmg}` + `magsolv.out` | 逐極解（half-clean 後：主 `.db` + 主 `.rmg` + 主 log） |
| `db/geom/geom_full_R{300,700}/` | 幾何 db + IGES（**非結果**，無 `.dat`） |
| `db/geom/from_parasolid/hung_hexapole_full.db` | CAD 匯入真實幾何（97 vol，**非 primitive、非 mesh**，參考用） |

---

## ⛔ 不可用：`gap_200um` / `no_gap`（**兩層原因，別再誤用**）

**① 六顆 coil 網格不一致 → 管線直接 crash**
`extract_ansys_data.m` 以 coil1 的 N 開 `B = zeros(N,3,6)` 再塞其餘 coil：

| variant | 各 coil matched 節點 |
|---|---|
| `no_gap` | coil1/3/6 = 655,901（Jul 10, `mesh_graded.db`）；**coil2/4/5 = 853,500（Jul 8 `meshbg.db` 的 stale data，從未在現行網格重解）** |
| `gap_200um` | **三套網格混雜**：589,276 / 655,931 / 853,530 |

**② `no_gap` 根本不是原始設計幾何**
建在 `MT_Mesh_Graded_BaseGap.txt`（自述 `[ADDED 2026-07-09 EXPERIMENT]`、註解明寫「非 CAD baseline」），
與純版 `MT_Mesh_Graded.txt` 差 102 行：有 200 µm base slab、**下極鐵棒縮短 43→38 mm**（`ROD_LEN_LO`）、
`GLOW_R=40`、`GAP_HW_LO=5`。⚠ 舊表把 655,919 節點的 BaseGap 網格當 canonical，**已作廢**。

> 舊 `gap_200um` 表述（μ_eff 支撐座：下極 μ=52.02 / 上極 μ=59.364，WP \|B\| 下 5.9 / 上 7.5 mT）
> 僅供辨識舊資料，**不得用於新分析**；如何處置（保留 / 刪除）另案。

---

## 讀取注意
- **預設讀 `R500`**（config `default_variant`）；要 R300/R700 必在 `main.m` 明設 `VARIANT`。
- ⚠ **config 的 `R_norm` 固定 500 µm、不隨 variant 變**：`R_act`（單位向量）與 `Pc_base`（無因次）不受影響，
  但 `filter_iron_nodes` 的錐體頂點與 `build_V_matrix` 的 sensor 位置用的是絕對座標。
  R≤150 µm 取樣球碰不到鐵，故 current-base 的 R300/R700 結果不受影響；**要做 voltage 或加大 R_select 前必須先處理**。
- dataset：`*_all.dat`（全節點）/ `*_wp.dat`（WP 球區，含匯聚細網格）。
- source 慣例見 `.claude/rules/charge-model-source-convention.md`（hung = 全 sink → 全域翻號）。
