# matlab/ — MATLAB 分析程式碼

讀 FEM 場（APDL `.dat` / Maxwell `.fld`）做點電荷模型擬合、解轉移矩陣、出結果 PDF。

## 第一層 = **求解器分支**（不是 model）

```
matlab/
├── APDL/       ← MAPDL（純量位）求解的場
│   └── Calibration_using_FEM_modeling/   ← 三模型共用校正管線（有自己的 README，權威）
└── Maxwell/    ← Ansys Maxwell 求解的場（結構鏡像上者，但少一層 Calibration 夾）
```

兩分支的 `function/` 是**鏡像**（同檔名、同職責、只差讀檔那一支：`extract_ansys_data` vs
`extract_maxwell_data`）。**改動要兩邊同步**——尤其 `emit_results.m`、`build_V_matrix.m`、
`utils/pole_sensor_geometry.m`。

> 兩分支對同一顆 long2016 模型有 **~14% 的 ĝ_I 落差**（9.66 vs 10.75 mT/A），已排除材料 / 幾何 /
> 安匝 / 網格收斂，尚未結案。比較兩邊數字前先確認基準一致（同 variant、同 R、同 dataset）。

## 共用結構（兩分支皆同）

```
config/<model>/[<geom>/]mt_constants.m   幾何 + 物理常數 + 路由旗標（唯一常數來源）
function/                                管線函式庫（見各分支 function/README.md）
main/main.m                              driver：只改頂部 per-run 旗標
utils/pole_sensor_geometry.m             ★ sensor 貼附幾何**唯一來源**（CAD 實測錐體 + 0.41mm 氣隙）
utils/<model>/                           後處理 compute 腳本（讀 calib .mat → 算 → 存 computed .mat）
plot/<model>/<base>/                     繪圖腳本（讀 .mat → 出圖，不重運算）
data/<model>/.mat/                       所有 .mat（自描述：結果 + 產生它的設定）
results/<model>/{single,eighteen}/        結果 PDF（single=USE_BIAS false、eighteen=true）
figures/<model>/<base>/{single,eighteen,common}/   圖檔
common_path/ansys_path.m                 路徑 resolver（model-agnostic）
```

三段式分工（不可越界）：
**① 校正** `main/main.m` + `function/` → `data/.../calib_*.mat`
→ **② compute** `utils/<model>/` 讀 calib `.mat` → 存 computed `.mat`（**不重跑校正**）
→ **③ 畫圖** `plot/<model>/<base>/` 讀 `.mat` → `figures/`（**不做重運算**）

## 跑法

```matlab
% 改 main.m 頂部旗標即可切換，其餘不用動
MODEL='long2016_hexapole_halfcut'; GEOM='tip40um';
VARIANT='maxwell';        % APDL 分支用 'graded'
BASE='voltage';           % 'current' | 'voltage'
USE_BIAS=true;            % false=single(fix) / true=eighteen(18-param bias)
R_select=150e-6;          % 定案取樣半徑
```
跑完自動存 `.mat` 並呼叫 `emit_results` 出 PDF。

## 產物落點（規則）

- `.mat` → `data/<model>/.mat/`（**不要**再寫 `MATLAB_data/`，該夾 2026-06-26 已移除；
  `matlab_path()` resolver 已 deprecated）。
- PDF → `results/<model>/{single,eighteen}/`，**只放 `.pdf`**（`.tex/.aux/.log` 由 emit 端編完即刪）。
- 圖 → `figures/<model>/<base>/{single,eighteen,common}/`。
- ⚠ **論文圖不在這裡**：`figures/paper_fig_plot/{plot,data}` → `figures/paper_fig/Section*/`。

## 規則（動手前先讀全文）

- **禁止使用 `backup/` 的資料與程式**（`no-backup-data.md`）——常數一律 `model_config(model, geom)`。
- 全程 **actuator frame**（`actuator-frame.md`）、**all-source 號誌**（`charge-model-source-convention.md`）。
- **coil→paper map 是 per-model 不可互抄**（`pole-coil-numbering.md`）；驗證看 K̄_I 對角占優且全正。
- 擬合電流必須 = FEM 激發電流 1 A（`fit-current-matches-sim.md`）。
- 結構凍結：新增 / 改名 / 移動 / 刪除檔或夾**一律先問**（`calibration-shared-structure.md`）。
- 畫圖先問風格選項、輸出實檔覆蓋迭代（`figure-style.md` / `figure-output.md`）。
- 單位：ℓ̂ µm、場 mT、ĝ_I mT/A、ĝ_V mT/mV、V mV（`unit-reference.md`）。

## 相關

`APDL/Calibration_using_FEM_modeling/README.md`（管線總覽）、同夾 `function/README.md`（逐函式速查）、
`../CLAUDE.md`（資料夾架構地圖）、`../ANSYS_data/<model>/RESULTS_MAP.md`（讀結果前必查）。
