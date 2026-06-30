# apdl/long2016_hexapole_halfcut/sim/ — 求解腳本（baseline + 變體）

**用途**：半切 hexapole 的 ANSYS 求解 input。建幾何 + air domain + coil（SOURC36）+ MAG BC + /SOLU + 抽 WP 場。本層依「求解策略」分子資料夾。

> **分工慣例**：`sim/` 只放「mesh 完之後跑模擬」的 solve / RESUME 腳本；**mesh 產生統一用 `../mesh/` 的 mesh-only deck**（`MT_Mesh_Graded.txt` graded／`MT_Mesh_Baseline.txt` baseline smrt5／`MT_Mesh_LowerFilled*.txt`／`MT_Mesh_SinglePole.txt`）。舊散落的 mesh-only deck（`sim/mesh/`、`geom/mesh/`）已刪除。

**內容**：
- `baseline/` — `MT_Sim_P1..P6.txt`：Long Fei verbatim 基準（含 VADD 合併的正確 topology），6 顆極自激解。
- `gap100um_mueq/` — **舊單段式 μ_r 等效氣隙的唯一保留 deck**：6 顆 coil（MAT_PROT μ_r=56 模 7mm post + 上下各 50µm=100µm gap）。μ_eff=56 < 翻負臨界 ~72 → P2←P1 翻負（磁通匯進 P2）；coil2-6 的 gap100 場只此 deck 能產（`calib_gap100um.m` 在讀）。
  - 舊單段式 `gap200um_mueq/`(μ31)／`gap20um_mueq/`(μ156)／`gap2um_mueq/`(μ259) deck **已退場**（被 2 段式 `mueq_sweep` 取代；既有 `data/coilN/gap<X>um_mueq/` 仍在，其中 gap200 已是 2-seg 場）。
- `mueq_sweep/` — **2 段式 μ_eff gap sweep（2026-06-26 由 `gap_sweep_2seg` 併入，取代舊 μ∈{120,90,60} 單檔版）**：`MT_Sim_P1_mueq_sweep.txt`（P1 掃 gap∈{0,50,100,150,200}µm）+ `MT_Sim_coil2to6_gap200.txt`（coil2–6 各自激 @gap200）。RESUME `mesh_baseline.db`（baseline smrt5，同 baseline 網格）→ 加 `MAT_PROT_UP=5/MAT_PROT_LO=6` + EMODIF（上 3 protrusion μ_up／下 3 μ_lo，**上下分開**）→ coilN=+1 raw + BC → magsolv → 各匯 `data/coilN/gap<X>um_mueq/`。μ_eff 用 `doc/effective_permeability/effective_permeability.pdf` **兩段式串聯磁阻公式**（x=0 強制 280；gap200 μ_up59.3／μ_lo71.0，比舊單段 31 溫和）。\|B\|max 隨 gap 遞減 1.126→0.978 T（gap0=baseline 驗證 OK）。scratch CWD=`db/mueq_sweep`。
- `lower_filled/` — 下極填滿 full-cone 變體 solve：`MT_Sim_LowerFilled_coil1.txt`（RESUME graded `mesh_lowerfilled.db`）+ `MT_Sim_LowerFilled_smrt5_coil1.txt`（**2026-06-26**，RESUME 純 smrt5 `mesh_lowerfilled_smrt5.db`）。coil1=−1（P1 SOURCE，FEM 端翻號）→ `data/coil1/lower_filled/`。
- `graded/` — region-graded 密網格的 P1-as-SOURCE 求解：`MT_Sim_P1_graded.txt`（只匯 R<2mm wp）；`MT_Sim_P1_graded_p2reg.txt`（**[MODIFIED]** 加匯「P2 整根極區域」box x −55→5mm；Z 上界 `SPH_OFST+34e-3` 涵蓋 holder 頂，輸出 `*_p2reg_full`，給 P2 whole-pole 磁路圖補支撐座箭頭）。
- `resolve/` — RESUME P1 graded master `.db`、只改 coil 電流重解 P2..P6（跳過 meshing）。
- `singlepole/` — **單極模型** solve deck `MT_Sim_SinglePole.txt`：`RESUME` 由 `../../mesh/MT_Mesh_SinglePole.txt` 產的 `mesh_singlepole.db`（下極填回完整圓錐+支撐座+鐵柱、**均勻鐵件 0.3mm**）→ 設 coil 1A + `D,ALL,MAG,0` + `magsolv,3` → 匯 B 場到 `data/singlepole/`。mesh/solve 分離。
- `scripts/` — 產生上述 sim .txt 的 `.py` 輔助（`_build_gap100um_mu_eq.py`、`_build_from_long2016_verbatim.py`、`_generate_halfcut_sims.py`）。⚠ `_build_gap100um_mu_eq.py` 以「複製 `gap200um_mueq` decks」為模板，但該模板已退場 → 不可重跑，僅留作歷史（gap100 decks 已產）。

> **μ_eff 串聯磁阻公式**：`mu_eff = h / [ (h−g)/mu_steel + g ]`，h=7mm 頸、g=總氣隙、mu_steel=280。
> 對照：g=0→280（baseline）、2µm→259、20µm→156、**100µm→56**、200µm→31。氣隙 μ_r=1 主宰磁阻，故微小氣隙就大幅拉低 μ_eff。（此單段式公式現只剩 `gap100um_mueq/` deck 在用；gap2/20/200 deck 已退場，μ_eq sweep 改用 `mueq_sweep` 的**兩段式**公式。）
> **P2←P1（all-source）翻負臨界 ≈ μ_eff 72（總氣隙 ~72µm）**：μ_eff>72（小氣隙）讀正（回流）、<72（大氣隙）讀負（匯進 P2）。

**資料來源 / 流向**：幾何取自 `hexapole-long2016` Long2016 source verbatim；解出 `.dat`（場）/ `.db`（模型）存到 `ANSYS_data/long2016_hexapole_halfcut/coilN/{standard,graded,gap200um_mueq}/`。

**命名 / 慣例**：`MT_Sim_P<N>*`，6 顆極只差 `CURR_ARRAY`（excited coil=1 或依 source 慣例 ±1，其餘 0），其餘內容同步；`D,ALL,MAG,0` 邊界必存在於 `/SOLU` 前；source 慣例＝每顆激發極場朝 WP 射出（下極 P1/P3/P6=-1、上極 P2/P4/P5=+1）；改動標 `[ADDED]`/`[MODIFIED]`；改尺寸先對齊 CAD。

## Mesh 方式辨別（每個 deck 怎麼劃網格）

**兩個辨別管道**：
1. **每個 deck 頂部標頭** `! [MESH] <scheme> — <detail>`（緊接 `/CWD` 之後）。開檔即見。
2. **下表**（authoritative 索引）：

| deck 資料夾 | mesh 方式 | 細節 |
|---|---|---|
| `baseline/` | **baseline smrt5** | `smrt,5` free-tet SOLID96（Long2016 verbatim）；解出 `coilN/standard` |
| `gap100um_mueq/` | **baseline smrt5** | baseline mesh + `EMODIF` 6 protrusion → MAT_PROT μ_r=56（100µm gap，**不改 mesh**）；**6 顆 coil**（唯一保留的單段式 deck）|
| `mueq_sweep/` | **baseline smrt5** | RESUME `mesh_baseline.db`；EMODIF 上下 protrusion 各自 2 段式 μ_eff；P1 gap{0,50,100,150,200} + coil2–6 gap200 → `data/coilN/gap<X>um_mueq/` |
| `graded/` | **region-graded** | 極 0.3mm / WP sphere 0.3mm / yoke 1.5mm / outer air 4mm + 尖端 LESIZE 20µm；含 `_graded_p2reg`（加匯 P2 整根極區域）|
| `sensor_spheres/` | **baseline smrt5 + sensor 加密** | sim deck 自建 mesh：VOVLAP 6 顆 R0.16mm air 球 @ESIZE 0.04mm 進 V7，V7 維持 SmartSize |
| `resolve/` | **RESUME（不重劃）** | RESUME 既有 master `.db`、只改 coil 電流重解，跳過 meshing |
| `singlepole/` | **RESUME（不重劃）** | RESUME `mesh_singlepole.db`（由 `../../mesh/MT_Mesh_SinglePole.txt` 產；均勻鐵件 0.3mm + 空氣漸變）；解 1A → `data/singlepole/` |

> 慣例：新增 / 改 deck 時，務必在 `/CWD` 後補上 `! [MESH] <scheme>` 標頭，並同步本表。

**相關**：見 `../README.md`、`../../../CLAUDE.md`、`.claude/rules/{apdl-editing,fit-current-matches-sim,result-read-safety}.md`、`doc/workflows/apdl-fem-run.md`。
