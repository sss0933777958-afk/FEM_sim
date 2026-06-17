# comsol-livelink

用 MATLAB LiveLink 跑 `.mph`(DC/AC 皆可):場探針、頻率掃描、係數抽取、模型檢查。
一鍵啟動見 `run_matlab_with_comsol.ps1`(memory `feedback_comsol_automation`)。

## 何時用

- ANSYS 不夠細(渦電流、多 turn coil 分布電流)→ COMSOL 補
- 頻率掃描(`Frequency Domain` study)
- 跨 `.mph` 任務(diag / probe / extract / fix)

## 輸入

| 參數 | 範例 |
|---|---|
| `{topic}` | `long_fei` / `kuo_quadrupole` |
| `{mph_path}` | `magnetic_sim/ANSYS/main/mph/{topic}/long_fei_model_DC.mph`(or `_AC.mph`) |
| `{task_script}` | `magnetic_sim/ANSYS/main/comsol/{topic}/<task>.m`(實際任務) |
| `{wrapper}` | `magnetic_sim/ANSYS/main/comsol/{topic}/run_<task>.m`(bootstrap,addpath + mphstart + 呼叫真任務) |
| `{port}` | 預設 2036 |
| `{KeepServer}` | 旗標 — 跨任務重用 server 加這個 |

## 既有任務(`magnetic_sim/ANSYS/main/comsol/long_fei/`)

| 類型 | 腳本 |
|---|---|
| Smoke / 幾何檢查 | `diag_geom.m` |
| 材料檢查 | `check_materials.m`、`inspect_sigma.m` |
| Coil 激勵驗證 | `check_coil_excitation.m` |
| AC 修正 | `fix_ac_*.m` |
| 頻掃 | `fig44_P5act_freqsweep.m`、`extract_I_in_30freq.m` |
| 場探針 | `fig44_P5act_P1sensor.m`、`probe_wp_center.m`、`P1field_P5act.m` |

## .mph 命名 pattern

`magnetic_sim/ANSYS/main/mph/{topic}/<name>_{DC,AC}.mph` + `.bak.mph` 備份。
DC = magnetostatic, AC = Frequency Domain。

## 啟動方式(2 種,擇一)

### A. 一鍵(推薦)

```powershell
cd kuo\comsol
.\run_matlab_with_comsol.ps1 -ScriptName run_<task> -KeepServer
```

啟動 mphserver(port 2036)→ launch MATLAB `-batch run_<task>` → 跑完(KeepServer = 不關 server)。

### B. 手動分開(偵錯用)

```powershell
& "G:\my_workspace\software\COMSOL62\Multiphysics\bin\win64\comsolmphserver.exe" `
   -login auto -user kuo
# 另開 MATLAB shell:
matlab -batch "run('magnetic_sim/ANSYS/main/comsol/{topic}/run_<task>.m')"
```

## Bootstrap wrapper pattern(每個任務都要寫一個 `run_<task>.m`)

```matlab
addpath('G:\my_workspace\software\COMSOL62\Multiphysics\mli');
mphstart(2036);
cd('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\comsol\{topic}');
<task_script_name>;   % 不加 .m
```

## 7 個坑(memory `feedback_comsol_automation`)

1. `-mlroot` 引號被 PS `-ArgumentList` 吃掉 → 分開啟 server + `mphstart(2036)`
2. R2025b 沒 `bin\win32\matlab.exe`(legacy 路徑)→ 同 1
3. `-batch` stdout 在 PS 啟動下空白 → 腳本內 `diary('out.txt')`
4. PS here-string newline 不穩 → 一律用 `.m` wrapper 當單 token
5. **Multi-Turn Coil 自動建 `ccc1`** → `CoilType='Numeric'` 不要手動 `coil.create('ccc1', ...)`;
   軸對稱簡單 loop 改 `CoilType='Circular'`
6. `mphselectcoords` 偶爾回空 → 用 `mphselectbox(model, [xmin xmax ymin ymax zmin zmax]', 'domain')`
7. server 等 credentials 卡住 → 啟動加 `-login auto -user <name>`

## 產物

- [ ] task 輸出 `.mat` / `.csv` / 圖 → `magnetic_sim/ANSYS/main/comsol/{topic}/` 或 `magnetic_sim/ANSYS/main/figures/{topic}/`
- [ ] `diary('out.txt')` log

## 常見額外坑(本 session 新加)

- 看 kHz 頻率相依說「artifact」前**先 `mphinterp` 查 σ**(memory
  `feedback_comsol_multiturn_ac_artifact`:Long Fei AC mph 鋼 σ=7e6 是真導電,
  kHz 飄移是真實渦電流,不是 artifact)
- WP probe 取錯座標 → 假象「100× 偏差/斷路」;WP 不一定在 z=0(memory
  `project_long_fei_DC_comsol`)

## 適用 topic

任何 `.mph`(`long_fei` / `kuo_quadrupole` 已實戰;新 mph 沿這套即可)。
