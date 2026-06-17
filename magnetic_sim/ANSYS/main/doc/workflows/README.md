# magnetic_sim/ANSYS/main/ Workflows

## 👉 怎麼用(給使用者)

跟 Claude 講自然句就會啟動對應流程,**不用記檔名、不用組長指令**。
缺參數 Claude 會主動問;中途遇到 ⏸ 檢查點,你批准才下一步。

### Input / Modeling(Round 1)

| 你要做什麼 | 對 Claude 講(任一句) | 啟動的流程 |
|---|---|---|
| 從 SolidWorks 出 STEP / IGES / IGES_converted | 「出 STEP」「重 export 模型」「更新 IGES」 | [cad-export](cad-export.md) |
| 從 STEP 抽幾何尺寸做成參數表 | 「解析 STEP」「STEP 變參數表」「讀 STEP」 | [step-to-apdl](step-to-apdl.md) |
| 寫 APDL 幾何腳本(`MT_Geom_*.txt`) | 「建 APDL 幾何」「寫 MT_Geom」「新增 variant」 | [apdl-geom-build](apdl-geom-build.md) |
| 單獨檢查模型 | 「檢查模型」「跑 model-check §N」 | [model-check](model-check.md) |

### Simulation / Post-processing(Round 2)

| 你要做什麼 | 對 Claude 講(任一句) | 啟動的流程 |
|---|---|---|
| 跑 FEM(ANSYS 解單 / 多 coil) | 「跑 FEM」「跑 sim」「跑 coil N」「ANSYS solve」 | [apdl-fem-run](apdl-fem-run.md) |
| 從 `.rst` 抽 B 場 / PATH / grid 成 `.dat` | 「抽 B 場」「跑 postproc」「抽 PATH」 | [apdl-postproc](apdl-postproc.md) |
| 跑 COMSOL `.mph`(DC / AC) | 「跑 COMSOL」「LiveLink」「COMSOL 頻掃」 | [comsol-livelink](comsol-livelink.md) |
| 跑 SEMulator 3D 製程 | 「跑 SEMulator」「跑製程」 | [semulator-process](semulator-process.md) |
| 同步單顆 IGES↔IGES_converted | 「同步 IGES」「重新轉 IGES」 | [iges-sync-quick](iges-sync-quick.md) |

### Analysis / Visualization(Round 2)

| 你要做什麼 | 對 Claude 講(任一句) | 啟動的流程 |
|---|---|---|
| H1 / H2 sensor 場比值 | 「跑 H1/H2」「算 H1H2」「dipole sensor 比較」 | [h1h2-analysis](h1h2-analysis.md) |
| 等效電荷模型擬合 | 「fit J」「擬合電荷」「跑 multipole」 | [charge-model-fit](charge-model-fit.md) |
| B_S 6×6 / 4×4 轉移矩陣 | 「算 B 矩陣」「算 B_S」「sensor 矩陣」 | [bs-matrix-derive](bs-matrix-derive.md) |
| 畫場圖(contour / quiver / streamline / placement) | 「畫場圖」「畫 quiver」「畫 sensor 位置」 | [field-plot](field-plot.md) |

> ⚠ 沒命中觸發片語時 Claude 會即興做、結果可能跟 SOP 不一致。
> 不確定講什麼時用上表那幾句最保險。

---

## Pipeline 全圖

```
SLDPRT                                           ┌─→ [comsol-livelink]   (替代 / 補強)
  │  → [cad-export]              ⏸ §1 + §2      │
  ▼                                              │
STEP + IGES + IGES_converted ───→ [semulator-process]
  │  → [step-to-apdl]             ⏸ §3
  ▼
<basename>_params.md                           ┌─→ [iges-sync-quick] (零碎更新)
  │  → [apdl-geom-build]          ⏸ §4         │
  ▼                                            │
MT_Geom_<variant>.txt ─── [apdl-fem-run] ──→ .rst/.rmg
                              │
                              ▼
                         [apdl-postproc] ──→ .dat
                              │
        ┌─────────────┬───────┴────────┬─────────────┐
        ▼             ▼                ▼             ▼
 [h1h2-analysis] [charge-model-fit] [bs-matrix-derive] [field-plot]
```

⏸ = 使用者必須明確批准才下一步。

## 13 份 SOP 一覽

| 階段 | 文件 | 做什麼 |
|---|---|---|
| input | [cad-export.md](cad-export.md) | SolidWorks 一次出 STEP + IGES + IGES_converted |
| input | [step-to-apdl.md](step-to-apdl.md) | 解析 STEP → 尺寸參數表(三維交叉檢查) |
| input | [apdl-geom-build.md](apdl-geom-build.md) | 參數表 + 樣板 → `MT_Geom_<variant>.txt` |
| check | [model-check.md](model-check.md) | 4 節獨立檢查 routine |
| sim | [apdl-fem-run.md](apdl-fem-run.md) | 跑 ANSYS solve |
| sim | [apdl-postproc.md](apdl-postproc.md) | 抽 B 場 / PATH / grid 成 `.dat` |
| sim | [comsol-livelink.md](comsol-livelink.md) | COMSOL DC/AC + LiveLink 自動化 |
| sim | [semulator-process.md](semulator-process.md) | SEMulator 3D 製程(GUI-only) |
| sync | [iges-sync-quick.md](iges-sync-quick.md) | IGES↔IGES_converted 快速同步 |
| analysis | [h1h2-analysis.md](h1h2-analysis.md) | H1/H2 場比值(雙極對稱性) |
| analysis | [charge-model-fit.md](charge-model-fit.md) | 等效電荷模型擬合 |
| analysis | [bs-matrix-derive.md](bs-matrix-derive.md) | B̄_S → B_S → V_out/V_in 矩陣鏈 |
| viz | [field-plot.md](field-plot.md) | 場視覺化(5 種圖型) |

## 觸發機制

本資料夾所有 SOP 都接到 `.claude/rules/main-workflows.md` 自然語 trigger。
要新增 SOP 時:寫 `<新名>.md` → 在本 README 加列 → 在 `main-workflows.md` 加觸發節。
