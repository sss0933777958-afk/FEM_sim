# main/ — 工作目錄全域規則（給 Claude）

`magnetic_sim/ANSYS/main/` 是目前**唯一活躍的設計**：4-pole MEMS Quadrupole（原 `kuo/`，FEM 求解器 = ANSYS MAPDL）。
本檔是 Claude 在此目錄工作的導覽 + 規則。動手前先讀：要找什麼資料、放哪、用哪支 resolver、繪圖怎麼做。

> 補充規範（不在此重複）：輸出位置細表見 `../../../.claude/rules/main-workspace.md`；
> 操作 SOP（出 STEP / 跑 FEM / 抽場 / fit / 畫圖…）見 `doc/workflows/README.md`；
> 讀 ANSYS 結果防呆見 `../../../.claude/rules/result-read-safety.md`。

---

## 🔒 鐵則：不擅自更動檔案架構

**未經使用者明確指示，不得更動檔案／資料夾架構** —— 包含**移動、改名、刪除、新建資料夾、重組目錄、搬移檔案**。
- 需要這類動作時 → **先問使用者**，得到明確同意才做。
- **例外**：編輯既有檔案的**內文**（改程式/文件內容）不受此限；但**新建/刪除/移動檔案或資料夾**一律先問。
- 繪圖時「新增功能組資料夾」也屬架構變動 → 依下方繪圖規則，**開新組前先問**。

---

## 📄 鐵則：改動同步 README

**改了東西，就把「受影響的」README 一起更新**（不必每次重掃全部）。
- 改了某資料夾的**內容**（新增/刪除/修改裡面的檔）→ 更新**該夾的 `README.md`**（內容清單、用途若變）。
- 若**新增/改名/移動資料夾**（依上方鐵則須先問使用者）→ 一併更新**上層的索引 README** 與本檔 **`## 資料夾架構地圖`**。
- 範圍 = **受影響的那幾份**；整個 repo 每層都有 README（根到葉），別讓它們跟實況脫節。

---

## 資料夾架構地圖（9 個頂層夾）

每夾都有自己的 `README.md`，要動該夾前先讀。第二層幾乎都是 **model topic**：
`long2016_hexapole_halfcut`（主力）/ `kuo_quadrupole` / `zhang_quadrupole`（CAD 用 `long_fei`）。

| 資料夾 | 是什麼 / 放什麼 | 要找資料去這裡 |
|---|---|---|
| `CAD_model/` | SolidWorks 原檔 + STEP（幾何 **source of truth**） | 量尺寸、出圖前對齊 CAD |
| `IGES/` | ANSYS 匯出的 IGES（公尺 + mm 兩版） | 幾何匯出原件 |
| `IGES_converted/` | 單位轉換後給 ANSYS `IGESIN` 的 IGES | 建 mesh 前匯入用 |
| `apdl/` | APDL 腳本：`<model>/{geom,sim,postproc}/`（+ sweep） | 改幾何/參數/重跑 sim 的 input |
| `ANSYS_data/` | FEM 輸出 `<model>/<case>/`（`.dat` 場 / `.db` 模型 / `.cdb`） | **讀 FEM 結果**（.dat） |
| `matlab/` | MATLAB 分析碼 `<model>/<功能組>/code/...` + `figures/` + `results/` | 跑分析、畫圖、resolver |
| `MATLAB_data/` | 分析成果 `<model>/<功能>/`（`.mat`/`.csv`/`.npz`） | **讀/寫分析結果**（.mat） |
| `doc/` | LaTeX 原稿 + 編譯 PDF + `workflows/`(SOP) | 推導、報告、流程 SOP |
| `.claude/` | Claude Code 本地設定（`settings.local.json`） | 非工作產物，通常不動 |

## 資料流 pipeline（一條龍）

```
CAD_model (SLDPRT/STEP)
   → IGES (ANSYS 匯出)  →  IGES_converted (單位轉換)
   → apdl/<model>/geom + sim (APDL input)
   → [ANSYS MAPDL 求解]
   → ANSYS_data/<model>/<case>/*.dat (場) + *.db (模型)
   → matlab/<model>/<功能組>/code (讀 .dat 做 fit/矩陣/校正/畫圖)
   → MATLAB_data/<model>/<功能>/*.mat (成果)  +  matlab/.../<功能組>/figures/*.png (圖)
   → doc/<主題>/ (LaTeX/PDF 報告)
```

## Resolver（路徑解析，不要硬寫絕對路徑）

都在 `matlab/<model>/common/`，相對自身定位（資料夾改名/搬移自動沿用）：
- `ansys_path('<model>'[, 'coilN', ...])` → `ANSYS_data/<model>/...`（讀 FEM `.dat/.db`）
- `matlab_path('<model>', '<功能>'[, file])` → `MATLAB_data/<model>/<功能>/...`（讀寫 `.mat`）

## matlab/ 功能組 schema

`matlab/<model>/` 第二層是**功能組**（activity），每組視性質含：
- 單一主程式組 → `code/main/main.m`（如 `Calibration using FEM modeling/{fix_l,no_fix_l}`）
- 多腳本組 → `code/scripts/`（如 `fixl_fit, bias_fit, bs_matrix, sensor_d, validation`）
- 純繪圖組 → `code/plot/`（如 `field_viz, sensor_placement`）
- 每組另有 `figures/`（圖）與 `results/`（auto-gen `.tex`）。

---

## 🎨 繪圖腳本規則（強制，畫任何圖前先讀）

1. **動手畫圖前，先跟使用者確認這支繪圖腳本屬於哪個功能組**（`field_viz` / `sensor_placement` / `flux_profile` / …），不要自己猜。
2. **不屬於任何現有功能組** → 要**新增一個功能組資料夾**（`matlab/<model>/<新組>/code/plot/` + `figures/`）；**開新組前也要先問使用者**。
3. **每個圖 / 每個繪圖任務只維護「一支」腳本**：在同一支腳本上**原地反覆修改**到使用者定案；**定案前不可另開新腳本**。
4. **使用者沒明講「新增（另一支）」就不開第二支腳本**（一個功能組底下可有多支，但各自對應一張已定案的圖）。
5. **圖檔要等使用者把腳本定案後，才存最終版到該組 `figures/`**；定案前一律用 MATLAB MCP **preview** 討論，不落地最終檔。
6. 沿用 repo 既有 Figure Production 慣例：**場圖一律畫真實 FEM 節點原值，不用 scatteredInterpolant / 格點內插**（除非使用者明確要求，且須在圖說標示為內插）。

> 一句話：**先問功能組 → 一任務一腳本、原地改 → 定案後才存圖**。
