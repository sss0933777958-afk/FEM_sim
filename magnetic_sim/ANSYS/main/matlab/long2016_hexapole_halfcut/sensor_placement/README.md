# long2016_hexapole_halfcut/sensor_placement — Hall sensor 幾何 / 位置 / 靈敏度繪圖（純繪圖）

**用途**：畫 Long Fei 下極半切六極上 **Hall (B_surface) sensor 的幾何、擺放位置與靈敏度**示意圖 —— 每顆極（P1..P6）sensor 放哪、正法向 n+ 方向、disc 取樣、上/下極不同擺法、不同 sensor 距離變體。為 B_S 矩陣推導做準備。**純繪圖功能組**。

**內容**：
- `code/plot/` — 例 `plot_pole_sensor_placement.m`（6 極 sensor 擺位側視）、`plot_upper/lower_pole_sensor_placement*.m`（上/下極 + 1p572 / 4p472 距離變體）、`plot_P1_geom_vectors_nhat_ahat.m`（n̂/â 幾何向量）、`plot_disc_sampling.m`（disc 取樣）、`plot_sensor_sensitivity.m`（靈敏度）、`plot_btip_placement.m`。
- `figures/` — 輸出 PNG（`Bsurf_placement_P1..P6.png`、`upper/lower_pole_sensor_placement*.png`、`P1_geom_vectors_*.png`）。

**資料來源 / 流向**：依幾何常數（`mt_constants`，極尖端座標 / 傾角 / WP frame）+ 視需要讀 ANSYS_data 場資料 → 畫真實節點 / 幾何 → PNG 到本組 `figures/`。

**命名 / 慣例**：純繪圖組 schema = `code/plot/` + `figures/`（無 `scripts/` / `results/`）；一張定案圖一支腳本。sensor 正法向 n+ 朝向見 `.claude/rules`（sensor 符號慣例：n+ 朝向 / 背離 WP，下極 n+=+z 出半切平面、上極 n+=cone 外法向）。場圖真實節點不內插。

**相關**：見 `../README.md`、`../../CLAUDE.md`。
