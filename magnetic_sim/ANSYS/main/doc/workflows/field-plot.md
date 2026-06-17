# field-plot

出 B 場圖:contour / quiver / streamline / sensor placement / coil layout。
**強制先 preview 才存檔**(`CLAUDE.md` Figure Production 規則)。

## 強制：資料來源 = 真實節點

場圖一律畫 FEM **節點原值**（節點實際位置 + 節點 Bx/By/Bz），**不可用 scatteredInterpolant / 格點內插，除非使用者明確要求**;要求內插時須在回覆/圖說標明「內插」,不可把內插當 raw。可讀性減量用**節點抽樣**（每格挑最近 y=0 平面的節點）,非內插。對應 memory `plot-real-nodes`;範例 `magnetic_sim/ANSYS/main/analysis/long2016_p1_only/plot_Bvector_p1_only_nodal.m`。

## 何時用

- 跑完 FEM / postproc 要視覺化驗證
- 出論文 / 報告圖
- 驗 sensor / coil 擺位

## 輸入

| 參數 | 範例 |
|---|---|
| `{topic}` | `long2016_hexapole_halfcut` |
| `{case_tag}` | `coil1_pre_fine_mesh` |
| `{plot_type}` | `contour` / `quiver` / `streamline` / `placement` / `layout` |
| `{plane}` | `xy` / `xz` / `yz`(2D);3D 為 `quiver_3d` |
| `{normalization}` | 絕對值 T;報告圖常用 mT 或 Gauss |

## 既有風格腳本

| 類型 | 樣板 |
|---|---|
| Contour | `magnetic_sim/ANSYS/main/analysis/kuo_quadrupole/plot/plot_Bcontour_xy_*.m` |
| Quiver | `plot_Bvector_topview_*.m`(top)、`plot_Bvector_sideview*.m`(side) |
| Streamline | `magnetic_sim/ANSYS/main/analysis/long2016_charge_fit/plot/plot_streamlines_yzero_with_charges.m` |
| Sensor placement | `magnetic_sim/ANSYS/main/analysis/long2016_hexapole_halfcut/plot/plot_pole_sensor_placement.m`<br>`plot_disc_sampling.m` |
| Coil layout | `magnetic_sim/ANSYS/main/analysis/kuo_quadrupole/plot/plot_ScaleDown_coil_layout.m`<br>`preview_coil_placement.m` |
| Flux check | `magnetic_sim/ANSYS/main/analysis/long2016_dipole_lower/check_flux_conservation.m` |

## 步驟(**強制順序**)

1. **討論 content** — 哪個 case、哪平面、軸範圍、normalize 用什麼?
2. **討論 style** — 字體大小、title 文字、legend、colormap、line thickness
3. **MATLAB MCP preview** — render 草圖,**先看不存**
4. ⏸ **使用者批准** content + style
5. **存到 `magnetic_sim/ANSYS/main/figures/{topic}/{case_tag}/<file>.png`**(`main-workspace.md` 規則)

## 產物

- [ ] preview 互動確認過
- [ ] `magnetic_sim/ANSYS/main/figures/{topic}/{case_tag}/<file>.png`(300 dpi)
- [ ] (報告圖)同步 mention path 進 `magnetic_sim/ANSYS/main/pdf/{topic}/scripts/` LaTeX

## 常見坑

- **跳過 preview 直接存** → 違反 `CLAUDE.md` Figure Production 規則
- 圖存到 `G:\my_workspace\report\` / git root → 違反 `main-workspace.md`(報告需要圖時
  改 reference `magnetic_sim/ANSYS/main/figures/...`)
- 標題用 LaTeX `\int` 等 → MATLAB OK 但本 chat 不渲染;用文字版「∫ B dz」
  (memory `feedback_no_latex`)
- 圖名帶方法後綴(`_jfit` / `_post` / `_fit`)→ 違反 `main-workspace.md`(case_tag 不可帶方法)
- preview / final 用不同腳本 → 後續找不到那張圖是誰生的

## 適用 topic

所有 topic。風格全 project 統一(`CLAUDE.md` 規則)。
