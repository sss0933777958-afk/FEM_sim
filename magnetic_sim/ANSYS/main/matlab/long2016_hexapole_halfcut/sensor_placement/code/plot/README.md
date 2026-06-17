# long2016_hexapole_halfcut/sensor_placement/code/plot — sensor 擺放繪圖腳本

**用途**：放 sensor_placement 的繪圖腳本（一張定案圖一支）：sensor 幾何 / 擺位 / 正法向 n+ / disc 取樣 / 靈敏度。

**內容**：`plot_pole_sensor_placement.m`（6 極側視擺位 + n+）、`plot_upper/lower_pole_sensor_placement*.m`（上/下極 + 1p572 / 4p472 距離變體）、`plot_P1/P2_geom_vectors_nhat_ahat.m`（n̂/â 向量）、`plot_disc_sampling.m`、`plot_sensor_sensitivity.m`、`plot_btip_placement.m`。

**資料來源 / 流向**：`mt_constants` 幾何 + 視需要 ANSYS_data 場 → 真實節點 / 幾何 → PNG 到 `../../figures/`。

**命名 / 慣例**：`plot_<極/部位>_<內容>.m`；PNG 放 `../../figures/`。**新增繪圖腳本前須依 `../../../../CLAUDE.md` 繪圖腳本規則**：先確認功能組、一任務一腳本（原地改到定案、定案前不另開）、定案後才存圖、真實節點不內插。

**相關**：見 `../../README.md`、`../../../../CLAUDE.md`。
