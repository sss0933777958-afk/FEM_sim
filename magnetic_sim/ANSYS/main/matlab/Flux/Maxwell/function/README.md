# function/ — 校正管線函式庫

這夾是三模型（long2016 / hung / NTU）通用校正管線的**函式庫**，全部被 `../main/main.m` 呼叫。
本檔是**逐函式速查**（每支在幹嘛、吃什麼、吐什麼、有什麼坑）；整體管線總覽、目錄結構、跑法請看
**上一層的 `../README.md`**，兩者互補、不重覆。

> 結構凍結：檔案的切分與職責已定案，**不得擅自新增/改名/搬移/合併**（改既有檔內文不受限）。
> [MODIFIED 2026-08-23 使用者拍板] 取樣器**收斂成兩支 `conv_design_*`**：
> `conv_design.m`（原本一支兩用）拆成 `conv_design_ws.m`（球）與 `conv_design_sensor.m`
> （圓柱），並把 `sphere_grid_sample.m` 與 `Sensor_grid_sample.m` **整支併進對應那支**
> （兩檔已刪除）。兩支各自自足、**不共用外部引擎檔** —— 三線性內插那段兩邊各一份，
> 是這個決定的代價（改一邊要同步另一邊）。
> 見 `.claude/rules/calibration-shared-structure.md`。

---

## 管線一眼看懂（誰呼叫誰）

```
model_config
  └► extract_ansys_data ─(import_ansys_data)─► raw 場 [measure, T]
       └► build_actuator_data ─(filter_iron_nodes)─► ad [actuator, mT, all-source]
            └► select_ball(ad, R) ─► fitting ─► ℓ̂, e
                 ├─ current:  solve_current(F) ───────────────────► K̄_I, ĝ_I
                 └─ voltage:  build_V_matrix ─► solve_voltage(V) ──► D̄, ĝ_V
                      └► save data/<model>/.mat/ ─► emit_results ─(emit_tex)─► results/<model>/{single,eighteen}/*.pdf
```
- `build_actuator_data` 回傳的資料包這裡叫 `ad`（**不是** D̄ 矩陣 —— D̄ 是 `solve_voltage` 的產物）。
- ~~`interp_field_to_points`~~ 已於 2026-08-23 刪除（見下表）；`INTERP_TO≠''` 現在會直接報錯。
- 括號內 `import_ansys_data` / `filter_iron_nodes` 是內部 helper，不由 main 直接呼叫。

---

## 逐函式

### 設定 / 分派
| 檔 | 簽章 | 做什麼 |
|---|---|---|
| **model_config.m** | `cfg = model_config(model, geom)` | Dispatcher。載 `config/<model>/[<geom>/]mt_constants.m`（先 `clear mt_constants` 防同名快取），疊上共用預設（`N_I=6`），掛上 `cfg.select_ball` handle，檢查路由欄位（`strategy`/`apdl_to_paper_idx`/`default_variant`/`regions`）。**`select_ball(D,R)` 是它的 local**：取 `|p|<R` 的節點、把各 coil 場疊成 `Bstack`（3Np×6）。|

### 讀資料（純讀、不加工）
| 檔 | 簽章 | 做什麼 |
|---|---|---|
| **extract_ansys_data.m** | `raw = extract_ansys_data(cfg, dataset, variant)` | 迴圈 coil1..6 呼叫 `import_ansys_data`，回**measure/WP frame、特斯拉**的 raw 場（`.x/.y/.z`、`.B` N×3×6、`.node_id`）。**不濾鐵、不旋轉、不翻號、不換單位**——那些是 `build_actuator_data` 的事。印 `Matched N`（讀結果核指紋用）。|
| **import_ansys_data.m** | `data = import_ansys_data(results_dir, dataset, coil_name)` | 讀單一 coil 的 `*_coord_*.dat` + `*_bfield_*.dat`，處理兩個檔案怪癖（coord banner NaN 欄、bfield 連號負值黏在一起），依 node_id merge，回 `node_id,x,y,z[m],bx,by,bz,bsum[T]`。extract 的 helper。|

### 幾何 / 座標系
| 檔 | 簽章 | 做什麼 |
|---|---|---|
| **build_actuator_data.m** | `ad = build_actuator_data(raw, cfg)` | raw（measure、T）→ **actuator frame、mT、全 source** 的資料包 `ad`（**非 D̄ 矩陣**）。依 `cfg.strategy` 分兩策略：`hex_magic`（long2016/hung）用 `cfg.R_act` 旋轉 P、B + `filter_iron_nodes` 濾鐵 + ×1e3 + 套 all-source 號誌 `s_source(apdl_to_paper_idx)`；`ntu_flat`（NTU）取 WP 載荷球、不旋轉。回 `.Pa[m]/.r2/.Ba[mT]/.R_act/.Pc_base/.F`。斷言 `R_act` 為正交旋轉。|
| **filter_iron_nodes.m** | `[air_mask,dbg] = filter_iron_nodes(x,y,z,c,opts)` | 幾何錐體模型：每根極建錐包絡（tip/base 半徑 + 錐長）+ 尖端安全球，回 air mask（true=保留）。build_actuator_data 的 helper。|

### 擬合（幾何無關優化器）
| 檔 | 簽章 | 做什麼 |
|---|---|---|
| **fitting.m** | `[e, l_hat, J] = fitting(P, Bstack, Pc_base, l0, USE_BIAS)` | Variable-projection 解有效長 `l_hat`（+ `USE_BIAS` 時 17 個偏移 `e`）。`lsqnonlin`：殘差裡由 `(l,e)` 建 kernel `S`、對每個激發欄 LS-profile 掉電荷（只迭代幾何 `l(+e)`）。`USE_BIAS=false`=single/fix、`true`=18-param/eighteen。回 `e(17×1)`、`l_hat[m]`、`J`（SSE[mT²]）。|

### 解 current 分支
| 檔 | 簽章 | 做什麼 |
|---|---|---|
| **solve_current.m** | `[KI_bar, gI_hat, G, rm] = solve_current(l_hat, e, Pc_base, P, Bstack, F)` | 重建 kernel、profile 電荷 `G`，`H_I=G·Fᵀ(FFᵀ)⁻¹`；`gI_hat=(6/5)H_I(1,1)[mT/A]`、gauge `KI_bar=(5/(6·G(1,1)))·H_I`（`K̄(1,1)=5/6`）。另算控制指標 `rm`（每點 `svd(S·Ĥ)`→𝒞/κ）+ 擬合 `RMSPE[%]`。|

### 解 voltage 分支
| 檔 | 簽章 | 做什麼 |
|---|---|---|
| **build_V_matrix.m** | `[V, exc_sign] = build_V_matrix(cfg, variant, raw, S_hall, SOFF_upper, n_uniform, sensor_r, axial_tol, sensor_override, V_METHOD)` | 電壓提取。sensor 幾何優先序：`sensor_override` > `cfg.sensor_pos/n` > 內建錐面公式。每 sensor 圓柱撒 `n_uniform` 點內插 raw 場，`V=S_hall·⟨B·n̂⟩[mV]`。兩內插法：`csv-tet`（預設，重建 tet mesh 重心內插，需 CSV mesh==solve mesh）、`scattered`（`scatteredInterpolant`，CSV≠solve mesh 時如 tip400）。最後套 all-source 欄號誌。|
| **solve_voltage.m** | `[D_bar, gV_hat, G, rm] = solve_voltage(l_hat, e, Pc_base, P, Bstack, V)` | solve_current 的電壓版，改由 6×6 sensor 電壓 `V[mV]` 驅動：`H_V=(G·Vᵀ)/(V·Vᵀ)`、`gV_hat=(6/5)H_V(1,1)[mT/mV]`、gauge `D_bar=(5/(6·H_V(1,1)))·H_V`（`D̄(1,1)=5/6`）+ 同樣 𝒞/κ。|

### 取樣器（產點：各維度數量 + 內插位置 + 配比增量）
兩支同一套方法（等分累積測度、取格心），差在幾何：球 vs 圓柱。都**只產查詢點的位置**，
配比決定「各維度要幾個點」，階梯決定「加點時該加哪一軸」。

| 檔 | 簽章 | 做什麼 |
|---|---|---|
| **conv_design_ws.m** | `[P,Bstack,tri,info] = conv_design_ws(N_r,N_φ,N_θ, R, opt)` | **球**：①設計階梯決定三整數 ②產點 ③濾鐵 ④規則格三線性內插取場 ⑤轉 actuator frame。配比固定 `N_r : N_φ : N_θ = 1 : 3 : 3π`；`r = R·u^(1/3)`、`φ = acos(1−2w)`、`θ = 2πv`（皆取格心 `(index−0.5)/N`）。回**已是 `fitting`/`solve_*` 要的形狀**（`P` Np×3、`Bstack` 3Np×6）；原始形狀場在 `info.B`。五種給法：`.ladder` 只出階梯表 / `.step` 取階梯第 q 級 / 明給三整數 / `.N_target`(或 `.c`) 由配比反解 / `.query` 任意點取場。**不做校正、不判收斂** —— 那兩件事在 `main.m`。|
| **conv_design_sensor.m** | `[x,y,z,B,tri,info] = conv_design_sensor(N_r,N_θ,N_z, opt)` | **圓柱**（Hall 感測器體積平均用）：同樣①階梯 ②產點 ③三線性內插取場，只差**不濾鐵、不轉 frame**（sensor 法線是全域框的量）。配比 `N_r : N_θ : N_z = 1 : 2π : √2·(H/R)`（**多一個形狀參數 H/R**，不是常數）；`r = R·√u`、`θ = 2πv`、`z = H·w`。⚠ 內建退化地板 `N_θ≥3`、`N_z≥2`（見檔頭實測表）。`.center` `.axis` `.span` 決定擺位；正交基底與 `build_V_matrix` 的亂數版逐字相同。`nargout<4` 時只產點、不讀場。**組 V 交給 `build_V_matrix`**。|
⚠ **配比階梯**（兩支各有一份逐字相同的 local `ladder_`）：
`s = t./w; [~,j] = min(s); t(j) = t(j)+1;`。因為格子邊長 `edge_i ∝ w_i/t_i = 1/s_i`，
`min(s)` 就是**最長的那條邊** → 切它一刀。平手時 `min` 取第一個索引 → 優先序固定。

⚠ 感測器圓柱 `R=150µm、H=100µm`（`main.m` 的 `sensor_r`/`axial_tol`）→ `H/R=0.667`、
配比 `1 : 6.283 : 0.943`。**N_θ 權重最大**，階梯絕大多數步都在加 N_θ；若場主要沿軸向變化，
誤差會由 N_z 主導，要用 `.NRTZ` 明給較大的 N_z。

---

### 內插（公平比較路徑）
| 檔 | 簽章 | 做什麼 |
|---|---|---|
| ~~**interp_field_to_points.m**~~ | — | **[REMOVED 2026-08-23]** 它呼叫只存在於 APDL 分支的 `extract_ansys_data`，在純 Maxwell path 上必定 `Undefined function`。`main.m` 的 `INTERP_TO` 分支改成明確報錯。要復活的話：來源是規則格，應改用 `conv_design_ws` 的 `.query` 模式做三線性，而不是把 APDL 的 `scatteredInterpolant` 版本搬回來。|

### 輸出（LaTeX→PDF）
| 檔 | 簽章 | 做什麼 |
|---|---|---|
| **emit_tex.m** | `T = emit_tex()` | 低階 LaTeX helper handle：`T.mat`（6×6 bmatrix）、`T.e`（3×6 偏移矩陣）、`T.scalar_unit`。數字格式慣例：`|指數|≥2` 才抽 `×10ⁿ`、正值不加 `+`。|
| **emit_results.m** | `emit_results(matfile)` | 讀自描述 `.mat` → 寫 LaTeX → inline `xelatex` → PDF 到 `results/<model>/<single|eighteen>/model_results_<base>_<variant>.pdf`。current 印 `K̄_I / ᴮĤ_I[mT/A] / ℓ̂ / G / F / ĝ_I / RMSPE`；voltage 印 `D̄ / ᴮĤ_V[mT/mV] / ℓ̂ / G / V / ĝ_V`（欄重排回 paper 序）；`USE_BIAS` 另印 `e[µm]`；都印 𝒞/κ。清中間檔只留 pdf。|

---

## 讀者必知的慣例（跨函式）

- **actuator frame**：整包載入→擬合→呈現一律在磁極軸座標系。六極用 `cfg.R_act` 旋轉、NTU 扁平板不旋轉（`.claude/rules/actuator-frame.md`）。
- **all-source 號誌**：每根極激發時 B 一律從尖端**射出**。對 raw 套 `cfg.s_source(apdl_to_paper_idx)`（只翻下極 sink）。判準＝K̄_I / D̄ 對角全正（`.claude/rules/charge-model-source-convention.md`）。
- **coil→paper map 是 per-model**：`cfg.apdl_to_paper_idx`（long2016/NTU `[1,3,6,5,2,4]`、hung identity）。**別互抄**（`.claude/rules/pole-coil-numbering.md`）。
- **單位**：ℓ̂ µm、b/場 mT、ĝ_I mT/A、ĝ_V mT/mV、V mV、`S_hall=130 mV/mT`（EQ-730L）。幾何 deck 內部仍 mm-magnitude（`.claude/rules/unit-reference.md`）。
- **取樣半徑 R**：`select_ball` 的 R **不在 config**，是 `main.m` 的 per-run `R_select`（預設 150µm）。
- **gauge 前物理矩陣**：結果 PDF 除 gauge 後 K̄_I / D̄，一律加印 `ᴮĤ_I = ĝ_I·K̄_I` / `ᴮĤ_V = ĝ_V·D̄`（`.claude/rules/calibration-transfer-matrix-output.md`）。
- **常數來源**：所有幾何/物理常數 + 方法旗標在 `config/<model>/[<geom>/]mt_constants.m`；pipeline 只消費、不假設魔術角。
- ⚠ **換機提醒**：`emit_results.m` 的 xelatex 路徑寫死（`C:\Users\Kuo\...\MiKTeX\...\xelatex.exe`），換機器要改。
