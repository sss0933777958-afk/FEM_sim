# h1h2-analysis

在 sensor placement 處(下極 cone surface + 0.41 mm normal offset)量 B,
算 H1 / H2 ratio,驗證雙極磁路對稱性 / 漏磁 / 反向耦合。

## 何時用

- 驗證新 dipole 設計兩端 sensor 對稱(理論 ratio = ±1)
- 比較同 / 反向激勵 H1 vs H2 漏磁(`*_opp.txt`)
- hexapole 推廣:對任兩相鄰 / 對極比較

## 輸入

| 參數 | 範例 |
|---|---|
| `{topic}` | `long2016_dipole_lower` / `long2016_p1p2_yoke` |
| `{case_tag}` | `coil1` / `coil1f`(forward)/ `coil1r`(reverse) |
| `{sensor_geom}` | 預設 Long Fei:cone surface + 0.41 mm normal(memory `project_long_fei_sensor_placement`) |
| `{coils}` | 要激勵的 coil 列表 |

## 步驟

1. **跑 FEM**(若沒跑過):走 [apdl-fem-run.md](apdl-fem-run.md),sim 腳本用
   `magnetic_sim/ANSYS/main/apdl/{topic}/sim/MT_Sim_H1H2.txt` 樣板
2. **抽 PATH**:走 [apdl-postproc.md](apdl-postproc.md) `H1H2` pattern
   → `magnetic_sim/ANSYS/main/ANSYS_data/{topic}/{case_tag}/coil*_{H1,H2}_*.dat`
3. **MATLAB 算 ratio**:
   ```matlab
   cd 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\analysis\{topic}'
   run('ratio_H1_H2.m')           % 既有腳本,讀 dat 算 H1/H2 → 印表
   % 或反向版本
   run('ratio_H1_H2_opp.m')
   ```
4. ⏸ **檢查點 → [`model-check.md` §3](model-check.md#3-參數表-sanity-check)** 變體:
   ratio 合不合理?(對稱設計理論 ≈ 1;不對稱要對得上電路圖)
5. **出表 + 圖**:走 [field-plot.md](field-plot.md) 出 `ratio_H1_H2_*.png`

## 既有腳本

| topic | ratio 腳本 |
|---|---|
| long2016_dipole_lower | `ratio_H1_H2.m` + `ratio_H1_H2_opp.m` |
| long2016_dipole_tilted | `ratio_H1_H2.m` |
| long2016_p1_only / p1p2 / p1p2_yoke / p1p2p3_yoke | 各有 `ratio_H1_H2.m` |
| long2016_h1h2 | `ratio_H1_H2_validation.m` |

新 topic 要複製 `analysis/long2016_dipole_lower/ratio_H1_H2.m` 當樣板。

## 產物

- [ ] `magnetic_sim/ANSYS/main/ANSYS_data/{topic}/{case_tag}/coil*_{H1,H2}_*.dat`
- [ ] MATLAB console 印的 H1 / H2 / ratio 數值
- [ ] (可選)`magnetic_sim/ANSYS/main/figures/{topic}/{case_tag}/ratio_H1_H2.png`

## 常見坑

- sensor normal 方向約定錯 → ratio 號型錯(memory `feedback_sensor_sign_convention_toward_wp`:
  上極 n+ ⊥ 傾斜錐面,不可抄下極)
- PATH 點 hardcode long_fei 4.572 / 0.41 mm → 新幾何要同步改
- 「同向激勵」vs「反向激勵」沒分清 → forward 用 `MT_Sim_H1H2.txt`,reverse 用 `*_opp.txt`
- raw 數值保留物理號;對稱化是 cosmetic 層的事(memory `project_long_fei_B_bar`)

## 適用 topic

dipole(主)、hexapole(per-pair 推廣)、quadrupole(對角 vs 鄰角)。
