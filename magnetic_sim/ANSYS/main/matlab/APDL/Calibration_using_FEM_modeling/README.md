# Calibration_using_FEM_modeling（三模型通用校正）

`long2016_hexapole_halfcut` / `hung_hexapole` / `NTU_hexapole` **共用**的一套校正碼。
**改一個參數（`MODEL`）切模型、改另一個（`VARIANT`）切資料類別、`BASE` 切 current/voltage**，跑同一支管線。

> ⚠ **結構凍結**：本夾的資料夾分層與檔案切分已定案，**不得擅自更動**；要改結構（新增/改名/移動/刪除資料夾或檔、
> 重組切分、加模型）**一律先問**。改既有檔內文不受限。見 `.claude/rules/calibration-shared-structure.md`。

## 管線（current + voltage）

```
main/main.m
  cfg = model_config(MODEL, GEOM)        載 config/<MODEL>/[<GEOM>/]mt_constants.m（幾何全在 config：R_act/Pc_base/sensor）
  raw = extract_ansys_data(cfg,DATASET,VARIANT)   純讀 .dat（measure, T）
  ── INTERP_TO='' ──  D=build_D(raw,cfg)（讀 cfg.R_act/Pc_base）→ [P,Bstack]=select_ball → Pc_base=D.Pc_base
  ── INTERP_TO≠'' ──  載參考 geom 點雲 → Bstack=interp_field_to_points(cfg,VARIANT,cfg.R_act,P,r_loc)；Pc_base=cfg.Pc_base
  [e,l_hat,J] = fitting(P,Bstack,Pc_base,l0,USE_BIAS)      優化器（e 開關；幾何無關）
  ── BASE=current ──  [KI_bar,gI_hat,G,rm] = solve_current(l_hat,e,Pc_base,P,Bstack,F)
  ── BASE=voltage ──  [V,~] = build_V_matrix(cfg,VARIANT,raw,...,V_METHOD)  （sensor 幾何來自 cfg；csv-tet|scattered）
                      [D_bar,gV_hat,G,rm] = solve_voltage(l_hat,e,Pc_base,P,Bstack,V)
  存 data/<model>/.mat/calib_<base>_<variant>_R<RRR>_<single|eighteen>.mat（自描述：結果+設定+GEOM/INTERP_TO/V_METHOD）
  emit_results(matfile)  → results/<model>/{eighteen|single}/model_results_<base>_<variant>.pdf
```
> **幾何通用化（2026-07-23）**：pipeline **不再假設魔術角**。幾何（`R_act`、`Pc_base`、per-pole `sensor_pos/n`）與
> 特例方法旗標（`v_method`、`interp_to`、`r_loc`、`sensor_r_loc`）**全部由 config 提供**，`build_D`/`build_V_matrix`
> 只「消費」。換一個 model / 角度 = 加一筆 `config/<model>/<geom>/mt_constants.m`，**零改 pipeline code**。

## 目錄結構（凍結）

```
config/<model>/[<geom>/]mt_constants.m   per-model[/幾何變體] 設定 + 幾何（R_act/Pc_base/sensor_pos/n）+ 方法旗標
                                         （long2016 分 tip40um/tip400um；hung/NTU 為 flat）
function/
  extract_ansys_data.m   純讀 .dat（raw, measure, T）
  build_D.m              raw → actuator frame（讀 cfg.R_act/Pc_base；不再假設魔術角）
  fitting.m              優化器（e 開關；build_S / make_Pc 為 local；幾何無關）
  solve_current.m        current：K̄_I / ĝ_I / 𝒞κ
  solve_voltage.m        voltage：D̄ / ĝ_V / 𝒞κ
  build_V_matrix.m       電壓提取（sensor 幾何 override>cfg>內建；V_METHOD csv-tet|scattered）
  model_config.m         dispatcher（model_config(model,geom)；variant 子夾+flat fallback；clear mt_constants 防快取）
  emit_tex.m             低階 LaTeX helper
  emit_results.m         讀自描述 .mat → LaTeX → PDF（current/voltage 兩分支）
  interp_field_to_points.m  把某變體 WP 場內插到參考 geom 的點雲（INTERP_TO；2026-07-23 由 utils 移入）
  import_ansys_data.m    純讀 coord/bfield .dat → struct（extract 的 helper；2026-07-23 由 backup 併入）
  filter_iron_nodes.m    幾何錐體濾鐵、回 air mask（build_D 的 helper；2026-07-23 由 backup 併入）
main/main.m              driver（頂部 per-run 調參 + BASE 開關）
common_path/ansys_path.m           共用路徑 resolver（model-agnostic；三模型載 FEM 場都經此，2026-07-23 由
                                   long2016/common 搬入，讓本夾自足、per-model 專案可退役）
data/<model>/{.mat,csv}/           結果存放（.mat 自描述；含保存的 single-pole .mat）
results/<model>/{eighteen,single}/ PDF（eighteen=USE_BIAS true、single=false）
utils/<model>/                     model-specific 一次性工具（如 tip 變體校正、hung R300/R700 疊圖 driver）
plot/<model>/<base>/               三模型 plot 腳本（<base>=current|voltage；不分 param。2026-07-23 由
                                   per-model 併入。⚠ 相依舊 per-model main_function/data、需改接才能重跑）
figures/<model>/<base>/{single,eighteen,common}/   三模型圖檔（single=fix、eighteen=bias、common=比較/診斷/場；
                                   由 legacy single_param/eighteen_param/shared 對應併入）
```

## 幾何變體：tip400（**走 main.m、無 bespoke driver**）

tip400（CNC 400µm 鈍尖、尖端後退）以前是 `utils/.../run_tip_calib.m`（已退役刪除）；現在 = **一筆 config + main.m 旗標**：
- `config/long2016_hexapole_halfcut/tip400um/mt_constants.m` 烘進：後退尖端、`R_act`=尖端到尖端（實測正交、=baseline
  magic-angle）、`Pc_base=R_act·dhat400`（非 canonical）、per-pole `sensor_pos/n`（tip+4.572·e2+(gap+0.41)·n̂、
  上極 CONE_ANG=35.49°、gap=rf(1−sinβ)）、`v_method='scattered'`（CSV≠solve mesh）、`interp_to='tip40um'`、`r_loc`/`sensor_r_loc`。
- 跑法：`main.m` 設 **`GEOM='tip400um', BASE='voltage', USE_BIAS=true`**（VARIANT/INTERP_TO/V_METHOD 自動由 cfg 填）。
- 結果 ℓ̂=3344µm/ĝ_V=5.23e−3（與舊 driver 逐位元一致）；tip40 baseline ℓ̂=857µm/ĝ_V=7.15e−3。
- `utils/` 現只留 `hung_hexapole/run_hung_current.m`（產 R300/R700 gain/iso 疊圖的正當 driver）。

## 現況（2026-07-22）

- **long2016**：current + voltage 兩路都完成、都與現有 `voltage_base`/`current_base` per-model 碼**逐位元等價**。
  - current：l̂=857.3µm、K̄_I 對角全正、K̄(1,1)=5/6、ĝ_I=9.66mT/A。
  - voltage：D̄(1,1)=5/6、ĝ_V=7.1466e-03mT/mV、`max|new.D̄−old.D̄|=0`。
- **NTU / hung**：config 已就緒；current 換 `MODEL` 即可跑（fit model-agnostic）；voltage 提取待 per-model 化。

## per-run 調參（在 main.m 頂部，不進 config）

`MODEL` / **`GEOM`**(config 幾何變體：long2016=tip40um|tip400um、hung/NTU='') / `VARIANT`(''→cfg.default_variant) /
`DATASET` / `BASE`(current|voltage) / `USE_BIAS` / `R_select`(150µm) / `l0`(0.5mm) / `I_actual`(1A)；
通用化旗標（''→由 cfg 填）：**`INTERP_TO`**(內插到參考 geom 點雲) / **`V_METHOD`**(csv-tet|scattered)；
voltage-only：`SOFF_upper`(4.572mm) / `n_uniform`(1e4) / `sensor_r`(0.15mm) / `axial_tol`(0.10mm)。

## 三模型差異的收法

- 純數值 / 幾何常數 / map / strategy → `config/<model>/mt_constants.m`。
- 檔名已正規化（NTU→coilN）。frame 物理（六極 magic-angle vs NTU 扁平板）→ `cfg.strategy`（`hex_magic`/`ntu_flat`，在 build_D/extract）。
- ✅ **本夾已完全自足、無跨專案相依**（2026-07-23）：pipeline 用到的 `import_ansys_data`/`filter_iron_nodes`
  已 copy 入 `function/`、`ansys_path` 入 `common_path/`，**不再 addpath `backup/…/analysis`**。
- ⚠ 三模型舊的 `<model>/Calibration_using_FEM_modeling/{current_base,voltage_base}/` 是**待刪除的 per-model 副本**。
  2026-07-23 已把其 **plot + figures** 併入本夾 `plot/`、`figures/`，**全部 `.mat`** 併入 `data/<model>/.mat/`
  （single-pole + 各 fit/calib，共 80 個；之後再整理），**`ansys_path.m`** 併入 `common_path/`；其餘
  （`code/main_function|function`、`main`、`results`、`reference`）隨 per-model 專案刪除。
- ⚠ 刪除範圍**限** `matlab/APDL/{long2016,hung,NTU}_hexapole*/`。`backup/` 樹**維持原狀不刪**（本夾已 copy 所需、
  不再依賴；backup 有自己的其他使用者）。
