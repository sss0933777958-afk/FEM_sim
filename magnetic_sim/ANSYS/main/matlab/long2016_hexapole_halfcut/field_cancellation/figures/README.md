# …/field_cancellation/figures/ — 圖檔

**用途**：source/sink 組合掃描的磁路圖輸出（直接放此層）。判準＝**WP 中心 R=50µm 球的 mean|B|**（實際操作體積）；圖為 **500µm xy 切面**（視覺脈絡，中央 50µm 才是評分區）。

**內容**（WP 中心 xy 切面 quiver；真實 FEM 節點不內插、格點抽樣、箭頭=面內單位方向、色=|B|(3D)；標 WP 中心 `+` 與 6 極方位 ±）：
- `field_xy_mincancel.png` — 最會抵銷 **[+ + + − − −]**（mean|B|_R50=0.196 mT，50µm 球內幾乎均勻近零 0.16–0.24mT）。
- `field_xy_max.png` — 最強 **[− − − − − −]（全同向）**（mean|B|_R50=32.2 mT）。
- `field_xy_allsource.png` — 全 source 參考 **[− + − + + −]**（每極尖端皆射出）：mean|B|_R50=2.38 mT（中心 0.37、50µm 邊緣升到 4.9mT）→ 中心零但有梯度，**非最佳**。

**判準結果（2026-06-17，over R=50µm sphere, 348 節點）**：`[+++−−−]` mean|B|=**0.196 mT**（唯一近零且均勻，比次佳 2.36mT 好 ~12×）；全 source `[−+−++−]` ~rank 3（2.38mT，中心零+梯度）；全同向 32.2mT。

**產生**：`../code/main/sweep_field_cancellation.m`（跑完自動畫三張）。raw FEM 疊加、符號=物理電流方向、不套 all-source。

**相關**：見上層 `../README.md`、`../code/plot/README.md`。
