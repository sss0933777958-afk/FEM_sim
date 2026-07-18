# 電荷模型符號慣例：全 source（強制讀取）

**使用者拍板（2026-06-12）**：做電荷模型擬合 / 場圖時，**每顆磁極激發時，磁場一律從該極尖端「射出去」（field out，全部當 source）**。

**🔑 拍板更新（2026-06-29，所有專案通用）**：canonical 做法 = **對 raw FEM 場直接套 `s_source=[-1,+1,-1,+1,+1,-1]`——只翻下極（sink）P1/P3/P6 成 source；上極（source）P2/P4/P5 原封不動、不翻**。**不要再用「全域 negate（六極全翻 −B_FEM）再額外補翻上極」那條 legacy 路**（它等價但易踩坑）。範本實作 = `Hall_sensor_base_{fix,no_fix}_dir/code/main/main_Dmatrix.m` 的 `as_sign`（`apdl_to_paper_idx(j)∈[2 4 5] → −1`，對 Bstack 翻上極激發欄；因 `load_coils_actuator` 已回全域 −B_FEM，「翻上極欄」= 把上極補回 +B_FEM = 只翻下極 sink）。**判準：真·全 source ⇒ 電荷自激發對角（D^v / K̂）全正；上極對角若跑出負 = 用了全域 negate，要翻回。**

對應 memory：`feedback_charge_model_source_convention.md`
相關文件：`magnetic_sim/ANSYS/backup/hexapole-long2016/docs/coil-winding-sign-convention.md`（記載 raw FEM 的 sink/source 物理）
相關規則：**`pole-coil-numbering.md`（編號 per-model；號誌錯常常是 map 抄錯造成的，兩者要一起查）**
相關 memory：[[sensor-sign-convention-toward-wp]]、[[long2016-halfcut-KI-fit]]

---

## 🔴 各 model 的 raw FEM 現況（**per-model，不可互抄**；2026-07-17 實測建表）

**raw FEM 是 sink 還 source 由該 deck 的電流方向決定，每個 model 不同。**

| model | deck 激發 | **raw FEM 尖端** | 正確處理 | 實測佐證 |
|---|---|---|---|---|
| **long2016** | `R,N,1,TURNS*1`（+） | **下極 sink / 上極 source（混合）** | 對 Bstack 套 `s_source=[-1,+1,-1,+1,+1,-1]`（只翻下極） | doc coil-winding §2/§5 |
| **hung** | `R,N,1,TURNS*1`（+） | **六極全 sink** | `load_coils_actuator.m:46` 的**全域 negate 一次到位**；**不得再 per-pole 翻號** | R700 六顆 coil 在 WP 的 raw B̂ 對各自 P k tip 方向內積**全 +0.999**；`diag(G)` 全正 |
| **NTU** | `R,N,1,-TURNS*...`（−，deck 註「電流反向」） | **六極全 source** ✅ 已符合新 model 標準 | 不需翻號（loader 也不該全域 negate） | full_assembly coil1 在 WP 的 raw B 與 P1 tip→WP 方向內積 **+0.957** |

⚠ **hung 的教訓**：`emit_model_results.m` 抄了 long2016 的 `coil_sign=[1 -1 1 -1 -1 1]` 去翻 hung 的上極 → 把已經正確的翻成錯的。**號誌處理必須依「該 model raw 的實測事實」，不是抄慣例。**

## ✅ 驗證判準（使用者拍板 2026-07-17：**要嘛就是看 K̄_I 的正負號**）

不要用眼睛猜、不要靠記憶。全 source + 正確 map ⇒ K̄_I 必須同時滿足：
- **對角占優**（`argmax|row|` = identity）← 驗 **map**（見 `pole-coil-numbering.md`）
- **對角全正**、**off-diag 全負**、**每列和 ≈ 0**（電荷中性，emergent） ← 驗 **號誌**
- `K̄(1,1) = 5/6`（gauge）

任何一條不過 ⇒ 停下來查 map 或號誌，**不要硬翻號把它壓成好看的樣子**。
K̂ 是自由 LS 解（`solve_KI_bar_gain.m` 唯一 gauge 是純量 `5/(6·g11)`）→ 此檢驗非循環論證。

**判準的獨立驗證法**（不靠擬合）：拿 raw `.dat`，在 WP 小球內取真實節點平均 B，與該極 tip 方向做內積 —— `>0` = B 指向尖端 = **sink**、`<0` = 從尖端射出 = **source**。
⚠ 但**對極共軸**（`dhat(P2) = −dhat(P1)`）：光看 WP 場方向**無法區分「P1 sink」與「P2 source」**，要靠 `G` 的對角占優（電荷落在哪顆尖端）才能定出是哪顆被激發。

---

## 🔒 核心慣例

**所有 6 顆極激發時，B 場方向都從尖端向外（source / 發散）。**

這是一個**呈現/建模符號慣例**（不是宣稱 raw 物理），對每顆極的 raw FEM B 套一個 per-pole sign `s_j` 使其朝外。

### raw FEM 的事實（不可改，來自 coil-winding doc §2/§5）

6 顆 coil 的 SOURC36 繞法相同、鐵芯內磁通都朝 +z；但因尖端相對 coil 的幾何位置不同：

| 層 | 紙上極名 | raw FEM 尖端磁通 | 要變 source 的 `s_j` |
|---|---|---|---|
| 下層 | **P1, P3, P6** | **SINK**（B 朝尖端進） | **−1**（翻號使朝外） |
| 上層 | **P2, P4, P5** | **SOURCE**（B 離尖端出） | **+1**（保持） |

→ **套在 raw FEM B 上的 source 慣例符號向量**（依 P1…P6 順序）：

```
s_source = [ -1, +1, -1, +1, +1, -1 ]
```

即「**翻下極 P1/P3/P6**」。這是論文「全 sink」慣例（`s_sink = [+1,-1,+1,-1,-1,+1]`）的整體 ×(−1)。

---

## canonical vs legacy（全域 negate）

**canonical（新 code 一律這樣）**：對 Bstack 直接套 `s_source`（只翻下極 sink）。實作 = `main_Dmatrix.m` 的 `as_sign`。
⇒ 電荷對角全正、off-diag 全負、每列和 ≈ 0（電荷中性 `K̂ = I − ones/6` 結構），不需任何「補翻」。

**legacy（等價但別再用）**：`load_coils_actuator` 回**全域 −B_FEM**（`s = [-1,-1,-1,-1,-1,-1]`）。在那 gauge 下下極對角正 ✓、**上極對角負 ✗**；要補回全 source 得**額外翻上極 P2/P4/P5 三欄**（或對 K̂/D 顯示時套 `coil_sign=[1 -1 1 -1 -1 1]`）。已驗證：legacy 補翻後 = canonical，**數值逐位相同**（field error 不變，純翻號）。

> ⚠ 上面這段的 raw sink/source 描述**只適用 long2016**。hung = 全 sink、NTU = 全 source（見檔頭現況表）。
> **新 model 一律在 deck 端就做成「raw 全 source」**，不靠後處理翻號（見 `pole-coil-numbering.md` 第 3 層）。

### 現況符合度（2026-06-30 更新）
| 狀態 | 腳本 | 做法 |
|---|---|---|
| ✅ canonical | `Hall_sensor_base_{fix,no_fix}_dir/main_Dmatrix.m` | flip-sink on Bstack（`as_sign`）→ D^v 對角全正 |
| ✅ canonical | **`fix_dir`（`load_coils.m` 2026-06-30 轉 flip-sink、`write_KI_tex.m` 移除顯示翻號）** | load_coils 只翻下極 sink → K̂_I 對角全正；存檔 Khat 直接 flip-sink。**已驗 G=gB·K̂·F == Hall D^v（max diff 7.4e-11，無需 coil_sign）**，PDF=`fix_dir/results/G_vs_Dv_gap200um_mueq.pdf` |
| ⚠ legacy 等價 | `no_fix_dir/main.m`（電荷 fit）、`sweep_KI_trend*.m`、`gen_KI_trend_per_radius.m` | 全域 −B_FEM + display 翻號（`coil_sign`）；數值等價、**暫不強制改**，但**新增/重寫時改用 canonical** |

> 注意：K̂/D^v 對角「正/負」本身含 gauge 自由度；**判斷對錯一律以「B 是否朝尖端外」物理方向為準**（全 source ⇒ 對角全正），不要只看對角正負硬壓。

---

## 強制規則

1. **所有專案**電荷模型擬合 / 場圖 / D 矩陣**預設一律全 source**（B 從每顆激發極尖端射出），除非使用者明確要論文 sink 慣例
2. **canonical = 只翻 sink（下極 P1/P3/P6）on Bstack（`as_sign`）；上極 source 不翻；不走全域 negate**（新 code 強制）
3. 下極 = P1/P3/P6、上極 = P2/P4/P5（依 coil-winding doc，不可記反）
4. legacy（全域 negate）需**額外翻上極 P2/P4/P5** 才等價；數值與 canonical 逐位相同（純翻號）
5. 不可硬壓對角全正當「物理」——必須回到「B 朝尖端外」的物理判據；off-diag 本來有結構（≈ −1/6）
6. 驗證法（doc §5）：看每顆 coil 在 WP center 的 B 方向，**6 顆都應背離激發尖端（朝外）**；或看電荷對角是否全正

---

## 觸發片語

- 「電荷模型符號」/「K̂ 對角為什麼負」/「磁極是 source 還 sink」
- 「跑 charge fit」/「畫激發場圖」時自動套用本慣例
- 「為什麼擬合收斂了還不物理」

## 何時不適用

- 使用者明確要求用論文 sink 慣例做對照
- 純幾何 / mesh / 後處理抽 .dat（不涉及場方向呈現）
