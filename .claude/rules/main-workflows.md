# magnetic_sim/ANSYS/main/doc/workflows/ 自然語觸發規則

當使用者用以下自然語要求對應動作時,**自動啟動該 SOP**。
工作涉及 `magnetic_sim/ANSYS/main/` 時搭配 `.claude/rules/main-workspace.md` 一起讀。

## 觸發 → SOP 對應

### cad-export.md(SolidWorks → STEP + IGES + IGES_converted)

觸發片語(任一即可):
- 「出 STEP」/「出 IGES」/「export 模型」
- 「從 SLDPRT 出檔」/「重新 export 幾何」
- 「更新 IGES」/「同步 IGES_converted」

啟動前主動問(缺什麼問什麼):
- `{topic}`(kuo_quadrupole / long_fei / zhang_quadrupole / 新 topic)
- `{sldprt}` 路徑(預設掃 `magnetic_sim/ANSYS/main/CAD/{topic}/SLDPRT/`)
- `{basename}`(預設用 sldprt 檔名去副檔名)

SOP 全文:`magnetic_sim/ANSYS/main/doc/workflows/cad-export.md`

---

### step-to-apdl.md(STEP → 參數表)

觸發片語:
- 「解析 STEP」/「讀 STEP」
- 「從 STEP 抽參數」/「STEP 變參數表」
- 「拆 STEP」/「分析 STEP 幾何」

啟動前主動問:
- `{topic}`
- `{step}` 路徑(預設掃 `magnetic_sim/ANSYS/main/CAD/{topic}/STEP/`)
- `{out}` 參數表輸出路徑(預設 `magnetic_sim/ANSYS/main/apdl/{topic}/geom/<basename>_params.md`)

SOP 全文:`magnetic_sim/ANSYS/main/doc/workflows/step-to-apdl.md`

---

### apdl-geom-build.md(參數表 → MT_Geom)

觸發片語:
- 「建 APDL 幾何」/「建 APDL 模型」
- 「寫 MT_Geom」/「做 geom 腳本」
- 「新增 variant」(在現有 topic 加新變體)

啟動前主動問:
- `{topic}`、`{variant}`(短 tag,例 `R0500_F20`)
- `{params}` 參數表路徑
- `{pole_config}`(quadrupole / hexapole / dipole / single)
- `{template}` 樣板腳本(從 apdl-geom-build.md 樣板表挑;不確定時讓使用者選)
- `{mu_r}`(常數,例 280)

SOP 全文:`magnetic_sim/ANSYS/main/doc/workflows/apdl-geom-build.md`

---

### model-check.md(檢查模型)

觸發片語:
- 「檢查模型」/「驗證模型」
- 「跑 model-check §<N>」(N = 1/2/3/4)
- 「看模型對不對」

通常**被其他 SOP 在 ⏸ 點自動呼叫**;單獨呼叫時問:
- 要跑哪一節?(§1 SolidWorks / §2 IGES round-trip / §3 參數表 / §4 APDL)
- 對應檢項的輸入(SLDPRT / IGES / 參數表 / .db 路徑)

SOP 全文:`magnetic_sim/ANSYS/main/doc/workflows/model-check.md`

---

### apdl-fem-run.md(跑 FEM)

觸發片語:
- 「跑 FEM」/「跑 sim」/「ANSYS solve」
- 「跑 coil <N>」/「batch sim」
- 「solve 模型」

啟動前主動問:
- `{topic}`、`{geom}`、`{sim_script}`(or 從樣板挑)
- `{coils}` 清單(預設 `[1]` 或全部 N 顆)
- `{case_tag}`(預設沿 `Lp<L_P>_T<thickness>_<excitation>` pattern)

SOP 全文:`magnetic_sim/ANSYS/main/doc/workflows/apdl-fem-run.md`

---

### apdl-postproc.md(後處理 / 抽 .dat)

觸發片語:
- 「抽 B 場」/「export bfield」
- 「post 處理」/「跑 postproc」
- 「抽 PATH」/「抽 grid」

啟動前主動問:
- `{topic}`、`{case_tag}`
- `{pattern}`(xyz_extract / FieldGrid / H1H2 / Diagnostics)
- 採樣範圍(box / grid 解析度 / PATH 點)

SOP 全文:`magnetic_sim/ANSYS/main/doc/workflows/apdl-postproc.md`

---

### h1h2-analysis.md(H1/H2 場比值)

觸發片語:
- 「跑 H1/H2」/「算 H1H2」/「H1/H2 比值」
- 「dipole sensor 比較」
- 「驗對稱性」(指雙極 sensor 對稱)

啟動前主動問:
- `{topic}`、`{case_tag}`、`{coils}`(forward / reverse)
- sensor 幾何(預設 Long Fei:cone surface + 0.41 mm normal)

SOP 全文:`magnetic_sim/ANSYS/main/doc/workflows/h1h2-analysis.md`

---

### charge-model-fit.md(等效電荷擬合)

觸發片語:
- 「fit J」/「跑 charge fit」/「擬合電荷」
- 「J fitting」/「跑 multipole」
- 「fit monopole/dipole/quadrupole」

啟動前主動問:
- `{topic}`、`{case_tag}`、`{sample_region}`(cube/cylinder/sphere)
- `{r_max}`、`{pole_shape}`(cylinder / cone-fillet — 決定 validity 安全範圍)
- `{model_order}`

**強制**:必先跑 validity sweep,err < 5% 才採信

SOP 全文:`magnetic_sim/ANSYS/main/doc/workflows/charge-model-fit.md`

---

### bs-matrix-derive.md(B_S 矩陣)

觸發片語:
- 「算 B 矩陣」/「算 B_S」/「B_bar 矩陣」
- 「轉移矩陣」/「sensor 矩陣」
- 「算 V_out/V_in」

啟動前主動問:
- `{topic}`、`{coil_results}` 路徑(N 個 coil 自激解齊)
- `{S_hall}`(預設 130 V/T)、`{I_in}`(預設 0.6 A)、`{k_A}`(預設 0.36)

**強制**:先跑 `diag_coil1_variants.m` 驗 coil1 資料夾陷阱

SOP 全文:`magnetic_sim/ANSYS/main/doc/workflows/bs-matrix-derive.md`

---

### field-plot.md(場視覺化)

觸發片語:
- 「畫場圖」/「畫 B field」/「畫 contour」
- 「畫 quiver」/「畫 streamline」
- 「畫 sensor 位置」/「畫 coil layout」

啟動前主動問:
- `{topic}`、`{case_tag}`、`{plot_type}`、`{plane}`、`{normalization}`

**強制順序**(`CLAUDE.md` Figure Production 規則):
1. 討論 content → 2. 討論 style → 3. MCP preview → 4. 使用者批准 → 5. 才存檔

SOP 全文:`magnetic_sim/ANSYS/main/doc/workflows/field-plot.md`

---

### comsol-livelink.md(COMSOL LiveLink)

觸發片語:
- 「跑 COMSOL」/「LiveLink」
- 「COMSOL 頻掃」/「COMSOL 抽 I/V」
- 「mphserver」/「跑 .mph」

啟動前主動問:
- `{topic}`、`{mph_path}`、`{task_script}`
- `{port}`(預設 2036)、`{KeepServer}` 旗標

**強制**:每個任務都要寫 `run_<task>.m` bootstrap wrapper,不能直接 `-batch <task>`

SOP 全文:`magnetic_sim/ANSYS/main/doc/workflows/comsol-livelink.md`

---

### semulator-process.md(SEMulator 製程)

觸發片語:
- 「跑 SEMulator」/「跑製程」
- 「SEMu build」/「製程模擬」

啟動前主動問:
- `{topic}`(目前只 `kuo_quadrupole`)
- `{model_file}`、`{process_file}`、`{layout}`

**現況**:paused at v2 / GUI-only;Claude 主要任務是引導 GUI + 更新 notes

SOP 全文:`magnetic_sim/ANSYS/main/doc/workflows/semulator-process.md`

---

### iges-sync-quick.md(IGES 快速同步)

觸發片語:
- 「同步 IGES」/「IGES_converted」
- 「重新轉 IGES」/「fix IGES 單位」

啟動前主動問:
- `{topic}`、`{iges_files}` 變動清單

**強制**:**kuo MKS 用 flag 2 / 重 export**,**不可抄 hung 的 sed 6→1**

SOP 全文:`magnetic_sim/ANSYS/main/doc/workflows/iges-sync-quick.md`

---

## 強制規則

1. 啟動前**先讀對應 SOP 全文**確認最新步驟,不靠記憶執行
2. ⏸ 檢查點**必須使用者明確批准**才下一步,不自行假設通過
3. 缺參數**主動問**,不假設預設值
4. 路徑放置遵守 `.claude/rules/main-workspace.md` 表格

## 不涉及 magnetic_sim/ANSYS/main/ 時忽略此規則
