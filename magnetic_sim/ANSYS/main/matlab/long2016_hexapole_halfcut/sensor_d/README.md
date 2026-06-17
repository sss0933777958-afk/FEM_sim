# long2016_hexapole_halfcut/sensor_d/ — Hall-sensor 模型 per-pole 常數 d 功能組

**用途**：Long Fei 半切六極 **Hall-sensor 模型 per-pole 常數 `d` 的求解 / 驗證 + LaTeX**。從 calibration 結果（`calib_sensor_d.mat`）求每極感測常數 `d`，含兩版：`d_v2`（無增益，`b_ij=S_i V_j d`，d~1e-3）與 `d_final`（含增益 `g_H=1/(4πℓ̂²)`，d~1e-8），兩者精確關係 `d_v2 = d_final × g_H`；另產 `K̂_H` 並做殘差 / NRMSE 驗證。

**內容**（功能組 root，下分 `code/scripts`，無 `results/`/`figures/`——LaTeX 輸出到 doc 樹）：
- `code/scripts/` — `gen_d_latex`（出 `d_v2.tex` + `d_final.tex`）、`gen_KH_latex`（出 `K̂_H` LaTeX）、`verify_d_*`（驗證 d 與 QS / Q-form 殘差）。

**資料來源 / 流向**：腳本讀 `MATLAB_data/long2016_hexapole_halfcut/charge_fit/calib_sensor_d.mat`（calib_fem.m page-2 產出，含 `d/gH/Vmat/ℓ̂/S_hall/...`）→ 算/排欄 → LaTeX 寫 `doc/charge_model_fitting/long2016_hexapole_halfcut/scripts/`。

**命名 / 慣例**：功能組 schema = `code/scripts`；`gen_*`=出 LaTeX、`verify_*`=驗證；欄序由「激發 coil」重排成「激發 paper pole P1..P6」（`apdl_to_paper_idx`）。

**相關**：見上層 `../README.md`、`../common/README.md`（resolver）、`../../../CLAUDE.md`。
