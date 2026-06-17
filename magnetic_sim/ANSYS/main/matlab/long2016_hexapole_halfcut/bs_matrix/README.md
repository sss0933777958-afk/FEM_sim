# long2016_hexapole_halfcut/bs_matrix/ — B̄_S / B_S / V_out_V_in 6×6 轉移矩陣功能組

**用途**：Long Fei 半切六極的 **6×6 sensor 轉移矩陣推導**。從 6 顆 coil 自激 FEM 解，在每個 Hall sensor 位置（0.3 mm 直徑感測盤面積平均）取沿 sensor 法向的 B，組出 `B̄_S`（每安培）、`B_S`（操作電流）、以及 `V_out/V_in`（含 Hall 增益 `S_hall=130 V/T`、放大器 `k_A=0.3614 A/V`、操作電流 `I_in=0.6 A`），並產出對應 LaTeX。

**內容**（功能組 root，下分 `code/scripts` + `results`）：
- `code/scripts/` — 矩陣計算 + LaTeX 產生器（`compute_B_*`、`gen_Vout_Vin_*`、`gen_Bs_latex` 等）。
- `results/` — auto-gen `.tex`（`B_bar_matrix_0p6A`、`B_s_*`、`Vout_Vin_*`，含多 sensor 距離變體與 gap200um_mueq 變體）。

**資料來源 / 流向**：腳本讀 `ANSYS_data/long2016_hexapole_halfcut/coilN`（FEM `.dat`，FEM 激發 0.6A → 直接對齊操作電流，免線性縮放）→ 感測盤面積平均算 B → `.mat` 寫 `MATLAB_data/<model>/bs_matrix/`（如 `B_bar.mat`）→ `.tex` 進 `results/` 與 `doc/Solve B_matrix/.../scripts/`。
- ⚠ sensor 資料夾陷阱（baseline vs gap200um_mueq vs pre_fine_mesh 等）：讀結果前依 `result-read-safety` 規則核指紋。

**命名 / 慣例**：功能組 schema = `code/scripts` + `results`（此組無 `plot/`）；`gen_Vout_Vin_<dist>` 的 `<dist>`（如 `4p572`）= sensor-to-tip 距離變體；`0p6A` = 操作電流；`results/` 只放 auto-gen `.tex`。

**相關**：見上層 `../README.md`、`../common/README.md`（resolver）、`../../../CLAUDE.md`、`.claude/rules/result-read-safety.md`。
