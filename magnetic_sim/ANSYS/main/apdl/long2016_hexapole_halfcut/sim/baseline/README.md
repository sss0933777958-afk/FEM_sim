# apdl/long2016_hexapole_halfcut/sim/baseline/ — Long Fei verbatim 基準求解

**用途**：半切 hexapole 的「正確 topology」基準解。6 支腳本各激發一顆極（self-excitation），建完整幾何 + air + coil + MAG BC + /SOLU，抽 WP 區 B 場。這是 RESULTS_MAP 認定的正確 coil1–6 來源。

**內容**：`MT_Sim_P1.txt` … `MT_Sim_P6.txt`（6 支，內容除 `CURR_ARRAY` 與 `/CWD`、`/OUTPUT` 檔名外完全相同）。每支建幾何（含下極 VSBV 半切）、mesh、解一顆極。

**材料表（2026-08-06 [MODIFIED]：磁極拆成 6 個 body）**：

| MAT | 內容 | murx |
|---|---|---|
| 1 `MAT_AIR1` | 內氣球 R ≤ 7 mm | 1 |
| 2 `MAT_MT` | **yoke + protrusion（不再含磁極）** | 280 |
| 3 `MAT_AIR2` | 外氣柱 R = 80 mm | 1 |
| **4–9 `MAT_P1`–`MAT_P6`** | **6 顆磁極各一個 MAT**（4=P1、5=P2、6=P3、7=P4、8=P5、9=P6） | **280** |

六顆磁極與 yoke **同為 1018 steel、murx 均 280**，拆 MAT 只是為了能逐極選取 / 繪圖 / 後處理（`ESEL,S,MAT,,MAT_P1`），**物理與解完全不變**。元素 component 同步新增 `COMP_P1`…`COMP_P6`；`COMP_MT` 仍涵蓋全部鋼件（`MAT_MT` + `MAT_P1`–`MAT_P6`）。

**資料來源 / 流向**：幾何/求解抄自 `hexapole-long2016/apdl/MT_Modeling_..._Coil1.txt`（已含 VADD 合併 yoke+prot、cone+yoke，避免 full 版被 strip 掉造成的 30-70% 弱磁通）；output → `ANSYS_data/long2016_hexapole_halfcut/coilN/standard/`（`all`≈494873 行 / `wp`≈390579 matched，|B|max ~1.0–1.14 T）。

**命名 / 慣例**：`MT_Sim_P<N>.txt`；6 支**只差 `CURR_ARRAY`**（coil N=1 其餘 0），其餘必同步；`D,ALL,MAG,0` 必在 `/SOLU` 前；保留 `!****` 原註解；FEM 激發電流＝1A（fit 時模型 I 必對齊，見 `fit-current-matches-sim`）；改動標 `[ADDED]`/`[MODIFIED]`。

**相關**：見 `../README.md`、`.claude/rules/{apdl-editing,fit-current-matches-sim,result-read-safety}.md`、`doc/workflows/apdl-fem-run.md`。
