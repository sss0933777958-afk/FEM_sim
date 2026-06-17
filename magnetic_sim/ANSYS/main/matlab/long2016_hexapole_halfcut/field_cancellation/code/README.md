# …/field_cancellation/code/ — 程式碼

**用途**：source/sink 組合掃描的所有 MATLAB 碼，依角色分組。

**內容**：
- `main/` — driver `sweep_field_cancellation.m`（載 6 coil → 取 WP 切面圓盤真實節點 → 掃 64 組合、依 min|B| 排序、印表 → 呼叫 plot 畫 min/max）。
- `plot/` — `plot_field_xy.m`（給一個 sign 組合，畫 WP 切面 xy quiver；真實節點不內插、格點抽樣、色=|B|、標 WP 中心 + 6 極方位與 ±）。

**資料來源 / 流向**：`main` 讀 `../../ANSYS_data/.../coil1..6/standard`（`ansys_path`）→ 疊加算場 → `plot` 出圖到 `../figures/`。

**命名 / 慣例**：`code/{main,plot}/`；driver 在 `main/`、繪圖在 `plot/`。

**相關**：見上層 `../README.md`。
