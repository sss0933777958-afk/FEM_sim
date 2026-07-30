# RESULTS_MAP — NTU_hexapole

**Single source of truth：這個 topic 底下每個 result dir 是什麼、能不能信、載入後該看到什麼指紋。**

讀任何結果前**先查這張表**（規則：`.claude/rules/result-read-safety.md`）。本表凌駕 memory。
建立：2026-07-22（NTU `.dat` 前綴正規化成 `coilN` 之後）。

> **⚠ 2026-07-22 檔名正規化**：NTU 的 `.dat` 前綴由「類別名」改成 **`coilN`**，與 long2016/hung 一致。
> 例：`data/full_assembly/coil1/full_assembly_coord_all.dat` → `.../coil1_coord_all.dat`。
> 讀取碼一律 `import_ansys_data(fullfile(root,<variant>,'coil<N>'), <dataset>, 'coil<N>')`
> （第 3 引數＝`coilN`，**不再是類別名**）。所有 NTU 分析腳本已同步更新（13 處）。

---

## 結構

`data/<variant>/coil<N>/coil<N>_{coord,bfield}_<dataset>.dat`
- **variant**（= 幾何/組裝案例）：`full_assembly`（6 極全組裝，主力）、`full_assembly_sleeve`（加 6 襯套）、
  `singlepole`（單極板，只 coil1）、`upper_assembly`（上層組裝，只 coil1）。
- **dataset**：NTU 只有 `all`（無 `wp`）；`singlepole`/`upper_assembly` 另有 `sensor_S1`/`sensor_S2`（板下/板上 sensor 加密球）。
- NTU 為**扁平板幾何、WP 在原點附近**（WP 形心 ≈ [0,0,9.3] mm）；WP 附近取樣靠 MATLAB 端 `R_load`（2mm 球）預篩，
  不靠 `_wp` 檔。載入策略 = `ntu_flat`（不旋轉、不濾鐵），見 `matlab/APDL/Calibration_using_FEM_modeling/`。

---

## CANONICAL

| variant | coils | dataset | matched 節點（coil1） | WP 球 \|B\| 指紋（coil1，2mm 球內 max） |
|---|---|---|---|---|
| `full_assembly` | coil1–6（**六 coil 節點一致** 1787195） | all | **1787195** | **~839 mT** |
| `full_assembly_sleeve` | coil1–6 | all | ~1806621（coord 1806639） | （μr=1 襯套 ≡ 空氣 ≡ full_assembly，物理等價） |
| `singlepole` | coil1 | all + sensor_S1/S2 | ~279238（coord 279256） | sensor 球：S1 1433 / S2 1431 節點 |
| `upper_assembly` | coil1 | all + sensor_S1/S2 | ~710335（coord 710353） | sensor 球：S1 1433 / S2 1431 節點 |

- **map**：NTU `apdl_to_paper_idx = [1,3,6,5,2,4]`（同 long2016，deck 下極先/上極後）。
- **號誌**：NTU deck 電流反向 → raw FEM 全 source（不需翻號）；載入 `-1e3`（T→mT）。
- 驗證（2026-07-22）：共用 extractor（`ntu_flat`）對 `full_assembly` 載出 N=34543（2mm 球）、
  與 NTU per-model `load_coils_actuator` 逐值相同（diff=0）。

---

## 怎麼用（層②核指紋）

`import_ansys_data(<dir>, <dataset>, 'coil<N>')` 讀 `<dir>/coil<N>_coord_<dataset>.dat` +
`coil<N>_bfield_<dataset>.dat`，印 `Matched %d nodes`。載入後對照上表 matched 節點數；
`full_assembly` vs `full_assembly_sleeve` 節點數不同可直接分辨。
