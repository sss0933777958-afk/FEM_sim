# apdl/hung_hexapole/gui/ — GUI 檢視 / render 用腳本（非求解）

**用途**：檢視 hung_hexapole **網格**（graded mesh 變體）的鋼件結構 —— 互動 GUI 或 batch 出 PNG。
不求解（無 `/SOLU`）；RESUME **拋棄式 `view.db` 副本**（不覆蓋 canonical mesh db）。**通用**：吃啟動 `-dir`，
可看任一 mesh 變體（如工作空間 R700 / R300）。

**內容**：
- `MT_View_SteelMesh.txt` — **互動 GUI 檢視**：RESUME `view.db` → `ESEL,S,MAT,,2`（鋼件 MAT_MT=2）→ `EPLOT`；
  header 含近-WP 6-tip 加密區 / 全域 / 內層細空氣的 hint。
  ⚠ **一定在獨立 scratch 夾開、不可用 mesh 夾當 -dir**（互動 GUI「Save」會寫回 -dir → 用 mesh 夾會覆蓋
  掉 `mesh_graded.db`，踩過）。啟動：`mkdir <scratch>`；`cp <mesh 夾>\mesh_graded.db <scratch>\view.db`；
  `MAPDL -g -dir <scratch> -j view -i MT_View_SteelMesh.txt`。
- `MT_Render_SteelMesh.txt` — **batch 多視角 PNG**（iso / front x-z / top x-y / 近-WP 6-tip 加密盒）→ 落 `-dir`。
  啟動：目標圖夾放 `view.db`（mesh_graded.db 副本），`MAPDL -b -dir <圖夾> -j vimg -i MT_Render_SteelMesh.txt -o vimg.out`。
  用於 `figures/hung_hexapole/mesh_R<tag>/`（`vimg000.png…`）。

**鋼件 = MAT_MT=2**（6 磁極 + 6 導柱/上核 + yoke）；空氣 = MAT_AIR1(內)/MAT_AIR_MID(中)/MAT_AIR2(外)。
hung 幾何有傾角（非平面板），故近-WP 檢視用「小盒 ±1.5mm」涵蓋 6 tip 匯聚區（R=0.3~0.7mm tip 皆含）。

**命名 / 慣例**：`MT_*`；不求解故無 `/SOLU`；互動 session 產的 `.lock/.log/.err` 為暫存（關閉後可清）。

**相關**：`../README.md`、`../mesh/MT_Mesh_Graded.txt`（產 mesh_graded.db 的 graded mesh deck）、
`apdl/NTU_hexapole/gui/`（同款範本）。
