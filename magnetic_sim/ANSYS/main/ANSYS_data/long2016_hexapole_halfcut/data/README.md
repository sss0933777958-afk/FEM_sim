# ANSYS_data/long2016_hexapole_halfcut/data/ — FEM 場結果（.dat）

**用途**：半切 hexapole 各 coil 各變體的 FEM 解出來的**場 / 座標 `.dat`**（交付檔）。
MATLAB 分析（Vmat / d / D_H / charge fit / 場圖）一律從這裡讀。

**結構**：`coil<N>/<variant>/coil<N>_<kind>_<dataset>.dat`
- `<N>` = 1..6（APDL coil；對照紙上極名 {P1,P3,P6,P5,P2,P4}）。
- `<kind>` = `coord`（節點座標 NLIST）/ `bfield`（B 場 PRNSOL）。
- `<dataset>` = `all`（全域 ~494,873 節點）/ `wp`（WP 區 ~390,579）/ `p2reg`（P2 整根極區域，graded 用）/ `circuit`（graded 的 P1 極區域）。

**變體（`<variant>`）**：
| variant | 是什麼 | 哪些 coil |
|---|---|---|
| `standard` | Long2016 verbatim **baseline**（smrt5、μ_r=280 全鋼） | coil1–6 |
| `sensor_spheres` | baseline + 6 顆 sensor 加密球（每 sensor 524~586 節點） | coil1 |
| `gap200um_mueq` | μ_r 等效 200µm 氣隙（protrusion μ_r=31） | coil1–6 |
| `gap20um_mueq` | μ_r 等效 20µm 氣隙（μ_r=156） | coil1 |
| `gap100um_mueq` | μ_r 等效 100µm 氣隙（μ_r=56；P2←P1 翻負區） | coil1–6 |
| `mueq_s1/s2/s3` | μ_eff 掃描 120/90/60（找 P2←P1 翻負臨界 ~μ72） | coil1 |
| `graded` | region-graded 密網格、P1-as-SOURCE（dataset `circuit`/`wp`） | coil1 |
| `graded_p2` | graded 密網格、**P2 整根極區域** dump（dataset `p2reg`，畫 P2 whole-pole 磁路用） | coil1 |

**讀取防呆 / 指紋**：載入前**先查 `../RESULTS_MAP.md`**（single source of truth：哪個 dir 可信、
載入後該看到的 matched 節點數 + |B|max）。⚠ `coilN` vs `coilN_*_mueq` 節點數相同，**靠 |B|max 區分**
（gap 版較低）。規則：`.claude/rules/result-read-safety.md`。

**慣例**：只放交付 `.dat`（歸檔白名單，不放 ANSYS 中間檔）；sim 解出後由 `apdl/.../sim/` 的 deck
`/CWD` 寫 scratch、再複製 `.dat` 進此處。

**相關**：見 `../RESULTS_MAP.md`、`../README.md`、`apdl/long2016_hexapole_halfcut/sim/README.md`。
