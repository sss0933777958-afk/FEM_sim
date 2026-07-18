# …/no_fix_l/code/main/ — no-fix-ℓ 主程式

**用途**：bias 校正的驅動程式（六極 18-param + 單極）。
**內容**：
- `main.m` — 六極（USE_BIAS 統一 fix / 18-param bias）：頂部 config（MODE `single`/`sweep`、R 50:5:500 µm 或單一 150、I_actual=1 A、TURNS、SHAPE=ball、ell0=0.5e-3、dataset='all'、VARIANT）；流程 `load_coils_actuator` → `select_ball` → `build_S_matrix` → `fitting` → `solve_KI_bar_gain` → `calc_ellipsoid` → 存 `.mat` → **`emit_model_results` 產 PDF**（末端一步）。
- `main_singlepole_bias.m` — 單極 **有-bias** 點電荷 fit `{ℓ̂, e_y, e_z, ĝ_I}`（電荷離軸橫向位移 `pc=ℓ·[1,e_y,e_z]`，ℓ 抓軸向、(e_y,e_z) 抓形狀不對稱）；3 尖端形狀 filled/halfcut/tipcut、R≤150µm、reuse fix_dir 的 `load_singlepole`+本層 `fit_singlepole_bias`。→ `data/singlepole_bias_fit.mat`+`results/singlepole_bias_fit.pdf`。對照 fix_dir 無-bias：halfcut 誤差 6.24%→1.03%。

**資料來源 / 流向**：讀 `ANSYS_data/long2016_hexapole_halfcut/`（經 `../../../common/ansys_path`）6-coil FEM（1 A）→ 每個 R 出一支 `.tex` 到 `../../results/no_fix_l/`（`fit_ball_R<R>um_<I>A.tex`）。

**命名 / 慣例**：driver 在 `main/`；模型函式一律在 `../main_function/`（唯一函式夾）。

**相關**：見上層 `../README.md`。
