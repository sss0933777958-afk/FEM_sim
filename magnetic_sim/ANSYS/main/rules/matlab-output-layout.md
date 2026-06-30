# 規則 #2：MATLAB `.mat` 放在「產生它的程式旁的 `data/`」

**使用者拍板（2026-06-26）**：MATLAB 分析產出的 `.mat`，放在**產生它的程式所在功能組旁邊的 `data/` 子資料夾**：

```
matlab/<model>/<activity>[/<subfolder>]/data/*.mat
```

不再丟中央 `MATLAB_data/<model>/...`。`.mat` 跟產它的 `code/` 同層、就近管理。

## 🔒 規則
1. **`Calibration_using_FEM_modeling` 的子資料夾一律建 `data/`**（`fix_dir`、`no_fix_dir`、`Hall_sensor_base_fix_dir`、`Hall_sensor_base_no_fix_dir`，即使暫無 `.mat` 也先建）。
2. **`matlab/<model>/` 其他活動夾**（`bs_matrix`、`bias_fit`、`flux_profile`…）：**有產 `.mat` 才建 `data/`**；純繪圖組（`field_viz`、`sensor_placement`、`field_cancellation`）不產 `.mat`、不必建。
3. **寫法**：腳本用自己已定義的 root 變數（`TREE`/`CAL`/相對自身）算出 `fullfile(<本功能組夾>,'data')`，**不要硬寫 `MATLAB_data` 絕對路徑**。
4. **跨夾讀取**：A 夾的腳本要讀 B 夾產的 `.mat` → 指 B 的 `data/`（相對或絕對），讀取方跟著產出方走。
5. **同步 README**：建/搬 `data/` 後，更新該功能組 `README.md` 的 pipeline／輸出路徑。

## 與既有慣例的關係（重要）
- **覆寫**：本規則**刻意覆寫**全域 `…/FEM_sim/.claude/rules/main-workspace.md`「`.mat → MATLAB_data/<model>/<功能>`」那條——限 `matlab/<model>` 產出物。
- **`MATLAB_data/` = legacy / 過渡**：尚未遷移的 shared namespace（如 `charge_fit/fitting_trend`、`fit_KI_ball`、`bs_matrix`、`flux_profile`、`freq_response`、`bh_saturation`、`calib_bias`）**暫留** `MATLAB_data/`；其 reader 暫不改。分波遷移。
- **`matlab_path()` resolver 不改**：仍解析到 `MATLAB_data/`，服務上述未遷移 / 跨 model 讀取。遷移後的腳本**改用 local `data/`**、不再經 `matlab_path` 寫自己的 `.mat`。
- **與 ANSYS_data 的 `data/` 區分**：
  - `ANSYS_data/<model>/data/` = FEM **`.dat` 場**（+ `.db/.cdb` 等交付）白名單（見 sim-cleanup「歸檔資料夾保留原則」）。
  - `matlab/<model>/<activity>/data/` = **MATLAB `.mat` 分析成果**（本規則）。兩者不同樹、不同內容物。

## 觸發片語（任一即套用）
- 「寫新 fit / 分析腳本要存 `.mat`」/「`.mat` 放哪」
- 「建 data 資料夾」/「`.mat` 放程式旁」
- 改既有 fit/矩陣/校正腳本的 `save()` 目的地時

## 何時不適用
- 純繪圖（只出 `.png`）/ 純讀 FEM `.dat`、不產 `.mat` 的腳本。
- 還沒分波遷移的 shared / 跨 model `.mat`（暫仍 `MATLAB_data/`）。

相關：`db-folder-retention.md`、`results-pdf-only.md`、全域 `main-workspace.md`、memory `feedback_matlab_local_data_layout`。
