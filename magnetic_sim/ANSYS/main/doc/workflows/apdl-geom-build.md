# apdl-geom-build

把 `step-to-apdl.md` 產出的參數表 + 既有腳本當樣板,寫成可跑的
`MT_Geom_<variant>.txt`(`/PREP7` → 幾何 → 材料 → 元素 → 網格 → SAVE)。

## 何時用

- 新 topic / 新 variant 要建幾何
- 既有 variant 改尺寸(small delta)→ copy + 改 *SET 區塊

## 輸入

| 參數 | 範例 |
|---|---|
| `{topic}` | `zhang_quadrupole` |
| `{variant}` | `R0500_F20` |
| `{params}` | `magnetic_sim/ANSYS/main/apdl/{topic}/geom/<basename>_params.md` |
| `{pole_config}` | `quadrupole` / `hexapole` / `dipole` / `single` |
| `{n_poles}` | 4 / 6 / 2 / 1 |
| `{template}` | 從下表挑既有腳本當起手樣板 |
| `{mu_r}` | 常數(例 280);**不可**用 B-H curve 若做線性疊加 |

## 樣板挑選

| pole_config | 既有可抄的腳本 |
|---|---|
| quadrupole | `magnetic_sim/ANSYS/main/apdl/kuo_quadrupole/geom/MT_Geom_R0500.txt`<br>`magnetic_sim/ANSYS/main/apdl/zhang_quadrupole/geom/MT_Geom.txt` |
| hexapole | `magnetic_sim/ANSYS/main/apdl/long2016_hexapole_halfcut/geom/MT_Geom_Export.txt`<br>`magnetic_sim/ANSYS/main/apdl/long2016_hexapole_full/geom/MT_Geom_Export.txt` |
| dipole | `magnetic_sim/ANSYS/main/apdl/long2016_dipole_lower/geom/MT_Geom_Export.txt`<br>`magnetic_sim/ANSYS/main/apdl/long2016_dipole_tilted/geom/MT_Geom_Export.txt` |
| single | `magnetic_sim/ANSYS/main/apdl/kuo_quadrupole/geom/MT_Geom_ScaleDown.txt`(改激勵單極) |

## 前置

- 共用 Pre-flight 5 項(尤其第 3:geom / sim 參數一致性 — 漂移就停手問)
- `step-to-apdl.md` 跑完 + model-check §3 通過
- 讀過 `.claude/rules/apdl-editing.md` + `simulation-constraints.md`

## 步驟

1. **挑樣板**(上表)→ 複製到 `magnetic_sim/ANSYS/main/apdl/{topic}/geom/MT_Geom_{variant}.txt`
2. **填 *SET 區塊** — 直接把 `{params}` 表的 `name=value`(MKS,m)灌入
3. **改幾何建構**:依 `{pole_config}` 用 CYL / SPH / BLOCK / VGEN / VROTAT
   建 pole + 氣域(inner sphere + outer cylinder)
4. **材料 / 元素** — SOLID96(volumes) + SOURC36(coils),`MP,MURX,1,{mu_r}`
5. **網格** — inner-sphere 細網控制 plot 平滑度
   (參考 memory `project_long2016_hexapole_halfcut` 的 ESIZE 10 µm 經驗)
6. **SAVE** + 可選 IGES export 自查
7. ⏸ **檢查點 → [`model-check.md` §4](model-check.md#4-apdl-模型驗證)**
8. **登記 variant** — 在 `magnetic_sim/ANSYS/main/apdl/{topic}/geom/_VARIANTS.md`(沒有就建)加一行

## 產物

- [ ] `magnetic_sim/ANSYS/main/apdl/{topic}/geom/MT_Geom_{variant}.txt`
- [ ] 跑一次出 `.db`(可選),確認 element count 在預期區間
- [ ] `_VARIANTS.md` 多一行

## 常見坑(完整 11 項見 memory `feedback_ansys_apdl_pitfalls`)

- `*` 前後加空格 → 被當註解
- mu_r 用 B-H curve + 做線性疊加 → 結果無物理意義(`simulation-constraints.md` §5)
- 缺 `D,ALL,MAG,0` 邊界 → 解不出 / 不唯一(`simulation-constraints.md` §7)
- SOURC36 上 / 下極 COIL_H 用了 ± 相反 → 上極自激極性反
  (前車之鑑:memory `project_long_fei_B_bar` coil1 陷阱)
- ANSYS 副產品(`.rst`/`.db`/`.out`)落在 `magnetic_sim/ANSYS/main/` 根 → 違反 `feedback_keep_topdirs_clean`
  → batch 一定加 `-dir <絕對路徑>`

## 適用 pole 配置

dipole(`long2016_dipole_*`)/ quadrupole(`kuo_quadrupole`、`zhang_quadrupole`)/
hexapole(`long2016_hexapole_halfcut`)/ single-pole(`long2016_p1_only`)。
