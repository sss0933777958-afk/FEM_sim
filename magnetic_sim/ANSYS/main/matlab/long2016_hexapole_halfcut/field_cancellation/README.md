# …/long2016_hexapole_halfcut/field_cancellation/ — 工作空間磁場抵銷掃描

**用途**：測「哪個 source/sink 電流組合會讓工作空間磁場抵銷」。每極只有 source/sink 兩態 ⇒ 對 6 極 raw FEM 場各乘 ±1，共 **2⁶=64 種**符號組合（等幅、純變號、**不套 all-source**、用真實物理極性）。用**真實 FEM 節點**評估每組合在 **WP 中心 R=50µm 球**（實際操作體積）的 **mean|B|**，找最會抵銷與最強的組合，在 500µm xy 切面畫磁路圖。

**內容**：`code/main/sweep_field_cancellation.m`（driver）、`code/plot/plot_field_xy.m`、`figures/`（3 張 xy：mincancel / max / allsource）。**無 `results/`**。

**判準**：**R=50µm 球 mean|B|**（使用者關注的操作體積；非 500µm、非中間切面、非 min|B|）。評分區(50µm 球) 與繪圖區(500µm 切面) 刻意不同。

**資料來源 / 流向**：讀 `../../ANSYS_data/long2016_hexapole_halfcut/coil1..6/standard`（`'all'`，6 個單線圈 1A 解，用 `../common/ansys_path`）→ 線性疊加 `B=Σ s_j·B_FEM,j` → 圖寫 `figures/`。表印 console，不存 `.mat`。

**結果（2026-06-17，over R=50µm sphere, 348 節點）**：
- **min cancel = [+ + + − − −]**：mean|B|=**0.196 mT**（球內幾乎均勻近零，比次佳 2.36mT 好 ~12×）。
- **max = [− − − − − −]（全同向）**：mean|B|=32.2 mT。
- 全 source `[−+−++−]`（每極尖端皆射出）≈ rank 3：mean|B|=2.38mT（中心零、50µm 邊緣升到 4.9mT）→ 中心抵銷+梯度，**非最佳**。

**命名 / 慣例**：raw FEM、真實節點不內插、符號=物理電流方向（不翻 all-source）；線性疊加合法（μ_r 常數、無飽和）；s 與 −s 在 |B| 上簡併（64→32）。

**相關**：見上層 `../README.md`、`../../../CLAUDE.md`（繪圖腳本規則）。
