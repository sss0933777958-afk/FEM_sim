# …/no_fix_l/code/main/ — no-fix-ℓ 主程式

**用途**：18-param bias 校正的唯一驅動程式。
**內容**：`main.m` — 頂部 config（MODE `single`/`sweep`、R 範圍 50:5:500 µm 或單一 150、I_actual=1 A、SHAPE=ball、ell0=0.5e-3、dataset='all'）；流程 load_coils_actuator → select_ball → fit_bias → make_Pc → gauge_KI → region_field_err → write_KbarI_tex。

**資料來源 / 流向**：讀 `ANSYS_data/long2016_hexapole_halfcut/`（經 `../../../common/ansys_path`）6-coil FEM（1 A）→ 每個 R 出一支 `.tex` 到 `../../results/no_fix_l/`（`fit_ball_R<R>um_<I>A.tex`）。

**命名 / 慣例**：單一主程式組 → 只有 `main.m`；數學在 `../function/`。

**相關**：見上層 `../README.md`。
