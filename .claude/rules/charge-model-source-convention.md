# 電荷模型符號慣例：全 source（強制讀取）

**使用者拍板（2026-06-12）**：做電荷模型擬合 / 場圖時，**每顆磁極激發時，磁場一律從該極尖端「射出去」（field out，全部當 source）**。

對應 memory：`feedback_charge_model_source_convention.md`
相關文件：`magnetic_sim/hexapole-long2016/docs/coil-winding-sign-convention.md`（記載 raw FEM 的 sink/source 物理）
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

## 套用到既有 trend 擬合（重要）

`magnetic_sim/ANSYS/main/analysis/long2016_hexapole_halfcut/fit/sweep_KI_trend*.m` 與
`magnetic_sim/ANSYS/main/doc/fitting_trend/scripts/gen_KI_trend_per_radius.m` **目前對 B 統一 negate**（`Bn = -[bx,by,bz]`，等於 `s = [-1,-1,-1,-1,-1,-1]`），並固定 `k̂11 = +5/6`。

在那個 gauge 下：
- 下極 P1/P3/P6 → 已成 source → 對角正 ✓
- 上極 P2/P4/P5 → 變 sink → 對角負 ✗

**要符合本慣例（全 source）**：在該 fit 結果上**再翻上極 P2/P4/P5 三欄**，
使 K̂ **對角全正、off-diag 全負、每列和 ≈ 0**（電荷中性 `K̂ = I − ones/6` 結構）。
field error 一個數字都不變，純翻號。

> 注意：K̂ 對角「正/負」本身含 gB gauge 自由度；**判斷對錯一律以「B 是否朝尖端外」這個物理方向為準**，不要只看對角正負。

---

## 強制規則

1. kuo 電荷模型擬合 / 場圖**預設一律全 source**（B 從每顆激發極尖端射出），除非使用者明確要論文 sink 慣例
2. 下極 = P1/P3/P6、上極 = P2/P4/P5（依 coil-winding doc，不可記反）
3. 既有 trend fit 因統一 negate，需**額外翻上極 P2/P4/P5** 才符合本慣例
4. 不可硬壓對角全正當「物理」——必須回到「B 朝尖端外」的物理判據；off-diag 本來有結構（≈ −1/6）
5. 驗證法（doc §5）：看每顆 coil 在 WP center 的 B 方向，**6 顆都應背離激發尖端（朝外）**

---

## 觸發片語

- 「電荷模型符號」/「K̂ 對角為什麼負」/「磁極是 source 還 sink」
- 「跑 charge fit」/「畫激發場圖」時自動套用本慣例
- 「為什麼擬合收斂了還不物理」

## 何時不適用

- 使用者明確要求用論文 sink 慣例做對照
- 純幾何 / mesh / 後處理抽 .dat（不涉及場方向呈現）
