# MATLAB_data/long2016_hexapole_halfcut — Long Fei 2016 下極半切六極（halfcut hexapole）的 MATLAB 分析成果

**用途**：本 model（主力設計，FEM 求解器 = ANSYS MAPDL）的所有 **MATLAB 分析輸出**（`.mat`/`.csv`/`.npz`，git 追蹤），與 **FEM 原始資料**（`ANSYS_data/long2016_hexapole_halfcut/` 的 `.dat`/`.db`，gitignore）分開存放。

**內容**：第二層 = 分析功能（activity），每個功能一個子夾，各自有 `README.md`：

| 功能子夾 | 是什麼 |
|---|---|
| `charge_fit/` | 點電荷模型 K_I 擬合 / Hall-sensor 校正 / 閉環驗證（含 `fit_KI_ball/`、`fitting_trend/` 兩個子層）|
| `bs_matrix/` | 6×6 B̄_S 感測矩陣與 V_out/V_in 轉移矩陣（含 sensor 距離掃描、gap/kA 變體）|
| `flux_profile/` | 沿極軸磁通剖面 Φ(s)（P1 / P2）|
| `freq_response/` | Fig 4.4 頻率響應 / cross-act（P1 self、P5→P1）|
| `bh_saturation/` | 鐵材 B-H 非線性飽和掃描 |

**資料來源 / 流向**：由 `matlab/long2016_hexapole_halfcut/<功能組>/` 的腳本讀 `ANSYS_data/long2016_hexapole_halfcut/coilN/*.dat`（1A FEM 場）做擬合 / 算矩陣 / 校正後寫入這裡。寫讀路徑用 resolver `matlab_path('long2016_hexapole_halfcut','<功能>')`；對應 FEM 用 `ansys_path('long2016_hexapole_halfcut','coil1','standard')`（兩支 resolver 在 `matlab/long2016_hexapole_halfcut/common/`）。

**相關**：見 [../README.md](../README.md)（MATLAB_data 總覽）與 [../../CLAUDE.md](../../CLAUDE.md)（main 工作目錄規則）。
