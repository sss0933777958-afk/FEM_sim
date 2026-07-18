# …/Hall_sensor_base_fix_dir/figures/ — 圖檔

**用途**：本包（求 d / sign 診斷 / P2 磁路）的 PNG 圖。皆由 `../code/plot/` 腳本產生（真實 FEM 節點、
不內插，除非檔名標 interp/schematic）。**改圖一律原地改該腳本重出，不手改 PNG。**

**圖檔家族 → 來源腳本**：
| 圖（檔名樣式） | 出自 | 看什麼 |
|---|---|---|
| `P2sensor_Braw_P1exc_<mode>[_standard][_raw].png` | `plot_P2sensor_Braw_P1exc.m` | P2 sensor 局部磁路箭頭（mode=circuit/zoom/pole/why/nodes） |
| `P2sensor_air_circuit_3d.png` / `P1P2_air_circuit_3d_*.png` | `plot_P1P2_air_circuit_3d.m`（FOCUS sensor/full） | 3D 空氣磁路（近 WP） |
| `P2pole_circuit_3d_graded_p2.png` | `plot_P1P2_air_circuit_3d.m`（FOCUS p2pole, graded） | P2 整根極 3D 磁路 |
| `P2pole_circuit_2d.png` | `plot_P2pole_circuit_2d.m` | **P2 整根磁路 2D 剖面（WP→支撐座）+ 鋼件輪廓**（讀 graded `p2reg_full`：Z 擴到 holder 頂，含支撐座上半段箭頭） |
| `P2_Bn_sign_map_standard.png` | `diag_P2_Bn_map.m` | P2 cone 上 B·n+ 匯進/回流 sign 分布 |
| `P2sensor_tets_3d.png` | `plot_P2sensor_tets_3d.m` | sensor 圓柱被真實 FEM tet 包住 |
| `interp_tet_schematic.png` | `plot_interp_tet_schematic.m` | 重心法內插示意（非真實資料） |

**相關**：見 `../code/plot/README.md`、`../code/README.md`。
