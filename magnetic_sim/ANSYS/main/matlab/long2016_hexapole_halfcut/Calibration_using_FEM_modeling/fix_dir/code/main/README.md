# …/fix_l/code/main/ — fix-ℓ 主程式

**用途**：fix-ℓ 校正的驅動程式（六極 + 單極）。
**內容**：
- `main.m` — 六極 fix-ℓ：頂部 config（MODE `single`/`sweep`、R 範圍 50:5:500 µm 或單一 150、I_actual=1 A、SHAPE=ball）；流程 load_coils → select_ball → fit_KI_fixl → region_field_err → write_KI_tex。
- `main_singlepole.m` — 單極 **無-bias** 點電荷 fit `{ℓ̂, ĝ_I}`（`load_singlepole`+`fit_singlepole`，R≤150µm、+x 框、SGN=+1）。頂部 `GROUP` 開關：`'shapes'`（filled/halfcut/tipcut 尖端形狀，預設）→ `data/singlepole_fit.mat`+`results/singlepole_fit.pdf`；`'gap_pos'`（支撐座氣隙同截面積 100mm²、位置 x=42.5/47.5）→ `data/singlepole_gap_pos.mat`+`results/singlepole_gap_pos.pdf`。（bias 版單極在 `no_fix_dir`。）

**資料來源 / 流向**：讀 `ANSYS_data/long2016_hexapole_halfcut/` 6-coil FEM（1 A）→ 每個 R 出一支 `.tex` 到 `../../results/fix_l/`（`fit_ball_R<R>um_<I>A.tex`）。

**命名 / 慣例**：單一主程式組 → 只有 `main.m`；數學在 `../function/`。

**相關**：見上層 `../README.md`。
