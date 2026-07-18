# apdl/NTU_hexapole/sim/ — 求解腳本

NTU_hexapole 的 ANSYS solve deck。mesh/solve 分離：mesh 由 `../mesh/` 產 `.db`，本層 deck `RESUME` 它再加線圈/BC/求解。兩支：單極 `singlepole/`、上層完整 `upper_assembly/`。

## `singlepole/MT_Sim_SinglePole.txt`（單磁極薄板 + 短導柱）
比照 `long2016_hexapole_halfcut/sim/singlepole/MT_Sim_SinglePole.txt` 的 RESUME 型 solve+extract。

- **RESUME** `ANSYS_data/NTU_hexapole/db/mesh_graded/mesh_graded.db`（幾何+網格+材料+2 顆 sensor 加密球）。
- **SOURC36 線圈**（mesh 未含 coil，於此加）：與短導柱同軸 (x=20.5,y=0,+z)、內半徑 3.5 / 外半徑 7.5 / 徑厚 4.0 / 高 4.2mm、z 中心 9.15mm（z 7.05–11.25）、**1 A × 50 匝、電流反向**（`R2=-TURNS*CURR`；使用者：磁路要往磁極跑→板內磁通 −x 朝極尖）。⚠ 半徑值（非直徑）；內緣 R3.5 > 導柱 R2.5 clear。
- **遠場 BC** `D,ALL,MAG,0`：外圓柱側面 r≈80mm + 上下端蓋 z≈±35mm。
- `magsolv,3`（DSP，278,887 eq，收斂、pivot 全正）。
- **抽 B 場**（`data/singlepole/coil1/`，PRNSOL,B,COMP）：全節點 `*_all` + 2 顆 sensor 球 `*_sensor_S1/S2`（coord + bfield 各一）。⚠ Fortran 定寬輸出：負 Bz 會吃掉與 By 間的空白 → 下游解析要用 regex/定寬、勿單純 split。

## 跑法
```bash
cd ANSYS_data/NTU_hexapole/db/singlepole
"G:/ANSYS Inc/v252/ansys/bin/winx64/MAPDL.exe" -b -np 1 -m 24000 \
  -dir <該夾> -j sim_singlepole -i apdl/NTU_hexapole/sim/singlepole/MT_Sim_SinglePole.txt -o solve.out
```
產物：`db/singlepole/sim_singlepole.db`(+主 .rmg, solve.out) ; `data/singlepole/coil1/*.dat`。

## 結果（1 A × 50 匝）
| sensor | mean \|B\| | mean B_z（翻電流後） | n+ |
|---|---|---|---|
| S1 (板下 z=−0.41) | 2.271 mT | −2.150 mT | −z |
| S2 (板上 z=+0.66) | 2.041 mT | +1.854 mT | +z |

翻電流後 B_z 相對舊版反號（板內 in-plane 磁通朝 −x 極尖＝往磁極跑）。1000 點面積平均 + B·n+ PDF = `matlab/NTU_hexapole/field_viz/`（⟨B·n+⟩ S1 +2.146 / S2 +1.853 mT）。

## `upper_assembly/MT_Sim_UpperAssembly.txt`（上層完整：3 極 + 6 導柱 + yoke）
**沿用單極單組激發**（使用者拍板）：只 post0(0°) 一顆線圈 + pole0 兩顆 sensor，設定與單極完全相同（同座標、同 R3.5/R7.5/4.2mm/z7.05–11.25、反向 1A×50 匝）。量 **yoke + 鄰極 + 6 導柱對 pole0 sensor 的磁路耦合影響**。

- **RESUME** `db/mesh_graded_upper/mesh_graded_upper.db`（見 `../mesh/`）。
- **遠場 BC** `D,ALL,MAG,0`：外圓柱側面 r≈100mm + 端蓋 z≈±55mm（enlarged 含 yoke）。
- `magsolv,3`（DSP，710,193 eq，收斂、0 error）。
- **抽 B 場** → `data/upper_assembly/coil1/upper_assembly_{coord,bfield}_{all,sensor_S1,sensor_S2}.dat`（S1 1421 / S2 1419 節點）。
- **跑法**：`-j sim_upper_assembly -m 40000`，產物 `db/upper_assembly/`（half-clean：.db + 主 .rmg + .out + magsolv.out）。

### 結果（1 A × 50 匝，post0 激發）＋ 對單極 baseline
| sensor | ⟨\|B\|⟩ upper | ⟨B·n+⟩ upper | ⟨\|B\|⟩ 單極 baseline | 倍率 |
|---|---|---|---|---|
| S1 (板下 z=−0.41, n+ −z) | **5.386 mT** | +5.109 mT | 2.271 mT | **2.37×** |
| S2 (板上 z=+0.66, n+ +z) | **4.316 mT** | +3.866 mT | 2.041 mT | **2.11×** |

yoke + 6 導柱 + 鄰極形成低磁阻 **return path** → pole0 板內磁通大增 → sensor 場 ~2.2–2.4× 於單極。B·n+ 皆正（沿 n+ 離板）、Bz 兩側反號（磁通穿板），與翻電流單極慣例一致。

## 後續（未做）
sensor 1000 點面積平均 postproc（比照單極 field_viz）；multi-coil / 每極矩陣（如需）。
