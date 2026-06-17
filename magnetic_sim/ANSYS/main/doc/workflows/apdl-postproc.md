# apdl-postproc

從 `.rst` / `.rmg` 抽 B 場、座標、PATH profile 等 `.dat`,給 MATLAB 端用。

## 何時用

- 跑完 FEM 要把資料拉到 MATLAB 做 fit / plot / matrix
- 想抽特定 PATH 或 probe grid 而不是整個 model

## 輸入

| 參數 | 範例 |
|---|---|
| `{topic}` | `long2016_hexapole_halfcut` |
| `{case_tag}` | `coilN` 或 `Lp0p46_T55_*` |
| `{pattern}` | `xyz_extract` / `FieldGrid` / `H1H2` / `Diagnostics`(見下表) |

## 4 個樣板

| pattern | 既有腳本 | 輸出 | 用途 |
|---|---|---|---|
| `xyz_extract` | `magnetic_sim/ANSYS/main/apdl/kuo_quadrupole/postproc/MT_Post_R0500_F20_xyz_extract.txt`<br>`magnetic_sim/ANSYS/main/apdl/zhang_quadrupole/postproc/MT_Post_Zhang_xyz_extract.txt` | `coil1_coord_all.dat` + `coil1_bfield_all.dat`(NLIST + PRNSOL) | 抽 ±400 µm cube 給 `scatteredInterpolant` / fit |
| `FieldGrid` | `magnetic_sim/ANSYS/main/apdl/long2016_p1_only/postproc/MT_Post_FieldGrid.txt`<br>`magnetic_sim/ANSYS/main/apdl/kuo_quadrupole/postproc/MT_Post_R0500_grid21.txt` | `<tag>_coord_all.dat` + `<tag>_bfield_all.dat`(21×21 grid) | probe grid |
| `H1H2` | `magnetic_sim/ANSYS/main/apdl/long2016_dipole_lower/postproc/MT_Post_H1H2.txt` | `coil1f_H1_*.dat`(PRPATH;cone surface + 0.41 mm offset) | 雙極 H1/H2 sensor profile,配 [h1h2-analysis.md](h1h2-analysis.md) |
| `Diagnostics` | `magnetic_sim/ANSYS/main/apdl/kuo_quadrupole/postproc/MT_Diagnose_Energy_Shells.txt`<br>`MT_Diagnose_Material_Split.txt` | 能量分殼 / 材料分區分布 | 驗 BC / 驗材料指派 |

## 前置

- `apdl-fem-run.md` 跑完(`.rst`/`.rmg`/`.db` 在 `magnetic_sim/ANSYS/main/ANSYS_data/{topic}/{case_tag}/coilN/`)
- 已知要抽什麼(grid / path / 整域)

## 步驟

1. **挑樣板**(上表)→ 複製到 `magnetic_sim/ANSYS/main/apdl/{topic}/postproc/MT_Post_<變體>.txt`
2. **改 *SET 區塊**:`case_tag`、`out_dir`、`jobname` 改成參數化
   (既有腳本 hardcode,改的時候順便修一下,SOP 要求 parametrise)
3. **跑 ANSYS batch**(skip-solve 模式,只讀 `.rmg`):
   ```powershell
   & $ANSYS -b -np 4 -m 24000 -dir $RES -j "coil$N" `
     -i "<absolute path to MT_Post_*.txt>" -o "$RES\post.out"
   ```
4. **檢查輸出 `.dat`** 行數合理(NLIST + PRNSOL 通常數萬行)

## 產物

- [ ] `magnetic_sim/ANSYS/main/ANSYS_data/{topic}/{case_tag}/coilN/<tag>_coord_all.dat`(座標)
- [ ] `magnetic_sim/ANSYS/main/ANSYS_data/{topic}/{case_tag}/coilN/<tag>_bfield_all.dat`(B 場)
- [ ] (H1H2)`coil*_{H1,H2}_*.dat` PATH 輸出
- [ ] `post.out` 無 `*** ERROR ***`

## 常見坑

- 腳本 hardcode `/CWD` → 改 case 要記得改路徑(SOP 要求參數化)
- 沒先 `RESUME, jobname, db` → PRNSOL 抓不到 result
- `xyz_extract` 取樣 box 太大 → 檔案動輒幾百 MB
- `H1H2` PATH 座標寫死成 long_fei 數值(4.572 / 0.41 mm)→ 改新設計要同步改

## 適用 topic

所有 topic。`H1H2` pattern 主用於 dipole / 上下極對比;`xyz_extract` 主用於 fit。
