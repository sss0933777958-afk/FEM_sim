# apdl-fem-run

用既有 `MT_Geom_*.txt` + sim 腳本跑 FEM,產出 `.rst` / `.rmg` / `.db` 到
`magnetic_sim/ANSYS/main/ANSYS_data/{topic}/{case_tag}/coilN/`。

## 何時用

- 新建好 geom 要試解
- 改 coil 激勵 / 改材料參數重跑
- 多 coil 批次掃(用 `apdl/sweep/<topic>/run_*.sh`)

## 輸入

| 參數 | 範例 |
|---|---|
| `{topic}` | `kuo_quadrupole` / `long2016_hexapole_halfcut` / `zhang_quadrupole` |
| `{geom}` | `magnetic_sim/ANSYS/main/apdl/{topic}/geom/MT_Geom_<variant>.txt` |
| `{sim_script}` | `magnetic_sim/ANSYS/main/apdl/{topic}/sim/MT_Sim_*.txt`(從樣板挑) |
| `{coils}` | `[1]` / `[1..6]` / `[1,3,6]` |
| `{case_tag}` | 例 `Lp0p46_T55_TURNS6`(沿 `main-workspace.md` pattern) |

## 樣板

| pole_config | 可抄 |
|---|---|
| quadrupole | `magnetic_sim/ANSYS/main/apdl/kuo_quadrupole/sim/MT_Sim_TURNS6_Final.txt` |
| hexapole | `magnetic_sim/ANSYS/main/apdl/long2016_hexapole_halfcut/sim/MT_Sim_P1.txt`(6 顆 P1-P6 各一檔) |
| dipole | `magnetic_sim/ANSYS/main/apdl/long2016_dipole_lower/sim/MT_Sim_dipole.txt` |
| H1/H2 變體 | `magnetic_sim/ANSYS/main/apdl/long2016_p1_only/sim/MT_Sim_H1H2.txt`(配 `apdl/postproc/.../MT_Post_H1H2.txt`) |

## 前置

- 共用 Pre-flight 5 項(尤其 #3 **geom/sim 參數一致性**;漂移就停手問)
- `apdl-geom-build.md` 跑完 + model-check §4 通過

## 步驟

1. **確認 geom 與 sim 參數一致**:`grep "POLE_R\|POLE_L\|TURNS\|R_sphere"` 兩檔比對
2. **準備結果目錄**:
   ```powershell
   $RES = "G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\{topic}\{case_tag}\coil$N"
   New-Item -ItemType Directory -Force $RES | Out-Null
   ```
3. **跑 ANSYS batch**(**`-dir` 一定要絕對路徑** — memory `feedback_keep_topdirs_clean`):
   ```powershell
   $ANSYS = "G:\ANSYS Inc\v252\ansys\bin\winx64\MAPDL.exe"
   & $ANSYS -b -np 4 -m 24000 -dir $RES -j "coil$N" `
     -i "<absolute path to {sim_script}>" `
     -o "$RES\solve.out"
   ```
4. **多 coil 批次**:複製 `magnetic_sim/ANSYS/main/apdl/{topic}/sweep/run_sweep.sh` 改 coil 列表
5. **檢查 solve.out** 結尾 `* END *` + 無 `*** ERROR ***`
6. **登記 case** 到 `magnetic_sim/ANSYS/main/ANSYS_data/{topic}/_CASES.md`(沒有就建)

## 產物

- [ ] `magnetic_sim/ANSYS/main/ANSYS_data/{topic}/{case_tag}/coilN/{coil$N}.{rst,rmg,db,out}`
- [ ] `solve.out` 正常結束
- [ ] `_CASES.md` 多一行

## 常見坑

- 忘記 `-dir <絕對>` → 副產品落 cwd / `magnetic_sim/ANSYS/main/` 根
- geom 跟 sim 參數不一致 → 解得出但物理錯(Pre-flight #3 攔截)
- 缺 `D,ALL,MAG,0` BC → DSP solver 不唯一解
- 上下極 SOURC36 winding 方向錯 → 對角符號反(memory `project_long_fei_B_bar`)

## 適用 topic

quadrupole / hexapole / dipole / single 全吃。
