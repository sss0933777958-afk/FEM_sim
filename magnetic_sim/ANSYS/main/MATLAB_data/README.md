# magnetic_sim/ANSYS/main/MATLAB_data — MATLAB 分析結果資料根

main 設計所有 **MATLAB 分析輸出**（`.mat` / `.csv` / `.npz`）集中放在這裡，與 **FEM 原始資料**
（`magnetic_sim/ANSYS/main/ANSYS_data/`，存 `.dat` / `.db` 等）**分開**。

```
資料流：ANSYS 求解 → magnetic_sim/ANSYS/main/ANSYS_data/<model>/...（.dat 場/.db 模型）
                       ↓ MATLAB 讀 .dat 做擬合/算矩陣/校正
                     magnetic_sim/ANSYS/main/MATLAB_data/<model>/<功能>/...（.mat 成果）
```

## 結構：`<model>/<功能>/`

第一層 = **物理模型**，第二層 = **分析功能**：

| model | 狀態 |
|---|---|
| `long2016_hexapole_halfcut/` | Long Fei 2016 下極半切六極（主力，有完整資料） |
| `kuo_quadrupole/` | Kuo 4-pole MEMS 四極（**空骨架**，FEM 資料已刪） |
| `zhang_quadrupole/` | Zhang 4-pole 四極（**空骨架**，FEM 資料已刪） |

## 功能子資料夾（以 `long2016_hexapole_halfcut/` 為例）

### `charge_fit/` — 點電荷 K_I 模型擬合 / 校正 / 驗證
- `fit_KI_full.mat`：document 版電荷模型擬合（K̂_I、ℓ̂、ĝ_B、R_a）。
- `calib_sensor_d.mat`：Hall-sensor 模型 per-pole 常數 d（含 g_H、Vmat、ℓ̂=0.856mm、殘差）。
- `calibration_final.mat`：最終校正（R*=150µm，ℓ̂=0.856、gB=8.43e-3、‖K̂‖=2.436）。
- `joint_6coil_40um_fit.mat`：R=40µm 球的 6-coil 聯合擬合。
- `validate_combos.mat` / `validate_combos_R150.mat`：電流組合閉環驗證（R150 為選定半徑）。
- `fit_KI_ball/`：各取樣半徑的 K_I 擬合 + 收斂（`fit_KI_R{040,050,100,150,500}.mat`、
  `fit_KI_ball_sweep.mat`、`KI_convergence_gB50.mat`、`KI_perpoint_gB50_R040.mat`）。
  （個別半徑 R200–R450 已刪——可由 `fit_KI_ball_sweep` / `fitting_trend/` 重生。）
- `fitting_trend/`：R 掃描趨勢與選取（`sweep_alln_vs_R.mat` 主檔、`KI_trend_sweep*`、
  `nrmse_window_maxmin`、`objective_R_select`、`Rlo/Rhi_decision`、`combo_nrmse` 等）。

### `bs_matrix/` — B̄_S / V_out_V_in 轉移矩陣
- `B_bar.mat`、`Bbar_S.mat`、`Bbar_S_4p572.mat`(+`_gap200um`)：6×6 B̄_S 感測矩陣
  （sensor 在錐尖後 4.572mm；gap 變體為氣隙對照）。
- `B_s_matrix.mat`、`B_s_final.mat`、`B_s_kA_final.mat`、`B_bar_s_corrected.mat`、`B_bar_S_toWP.mat`、
  `Bs.mat`：B_S 系列中間/修正版。
- `Vout_Vin_*.mat`：V_out/V_in 矩陣，含 **sensor 距離掃描**（1p572…4p572 mm）+ `_gap200um`/`_mueq` 變體
  + `V_out_V_in_kA0p3` / `Vout_Vin_kA0p36`（放大器增益 k_A 版）。

### `flux_profile/` — 沿軸磁通剖面 Φ(s)
- `P1_flux_profile.mat` / `P2_flux_profile.mat`（+ `_smrt4` 舊批次）：P1/P2 沿極軸磁通。

### `freq_response/` — Fig 4.4 頻率響應 / cross-act
- `fig44_P1.mat`、`fig44_P1_30freq.mat`、`fig44_P1_v2.mat`：P1 self-act 頻響。
- `P1field_P5act.mat`(+`_1kHz`)、`fig44_P5act_P1sensor.mat`、`fig44_P5act_freqsweep.mat`、
  `wp_P5act.mat`：P5 激發、P1 感測的 cross-act。
- `I_in_30freq.mat`：30 頻點輸入電流。`comsol_wp_plane_P1.mat`：COMSOL WP 平面對照。

### `bh_saturation/`
- `bh_saturation_sweep.mat`：鐵材 B-H 非線性飽和掃描結果。

## 路徑慣例
- MATLAB 結果：resolver `matlab_path('<model>','<功能>')`
  （例 `matlab_path('long2016_hexapole_halfcut','charge_fit')`）。
- 對應 FEM 資料：resolver `ansys_path('<model>','coil1','standard')`（→ `magnetic_sim/ANSYS/main/ANSYS_data/`）。
- 兩個 resolver 都在 `magnetic_sim/ANSYS/main/matlab/<model>/common/`（相對解析、資料夾名各集中一處，改名只改一行）。
- `.mat`/`.csv`/`.npz` 在 git 追蹤（小、是成果）；FEM 重產物（.dat/.db…）在 `magnetic_sim/ANSYS/main/ANSYS_data/` 被 gitignore。
