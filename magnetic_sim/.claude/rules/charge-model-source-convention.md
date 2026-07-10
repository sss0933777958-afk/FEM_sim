# 電荷模型符號慣例：全 source（強制讀取）

**使用者拍板（2026-06-12）**：做電荷模型擬合 / 場圖時，**每顆磁極激發時，磁場一律從該極尖端「射出去」（field out，全部當 source）**。

**🔑 拍板更新（2026-06-29，所有專案通用）**：canonical 做法 = **對 raw FEM 場直接套 `s_source=[-1,+1,-1,+1,+1,-1]`——只翻下極（sink）P1/P3/P6 成 source；上極（source）P2/P4/P5 原封不動、不翻**。**不要再用「全域 negate（六極全翻 −B_FEM）再額外補翻上極」那條 legacy 路**（它等價但易踩坑）。範本實作 = `Hall_sensor_base_{fix,no_fix}_dir/code/main/main_Dmatrix.m` 的 `as_sign`（`apdl_to_paper_idx(j)∈[2 4 5] → −1`，對 Bstack 翻上極激發欄；因 `load_coils_actuator` 已回全域 −B_FEM，「翻上極欄」= 把上極補回 +B_FEM = 只翻下極 sink）。**判準：真·全 source ⇒ 電荷自激發對角（D^v / K̂）全正；上極對角若跑出負 = 用了全域 negate，要翻回。**

對應 memory：`feedback_charge_model_source_convention.md`
相關文件：`magnetic_sim/ANSYS/backup/hexapole-long2016/docs/coil-winding-sign-convention.md`（記載 raw FEM 的 sink/source 物理）
相關 memory：[[sensor-sign-convention-toward-wp]]、[[long2016-halfcut-KI-fit]]

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
