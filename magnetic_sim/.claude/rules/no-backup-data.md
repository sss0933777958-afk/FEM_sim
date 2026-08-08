# 禁止使用 `backup/` 的資料與程式（強制讀取）

**使用者拍板（2026-08-08）**：**以後不要再使用 `backup/` 的資料。**
活躍工作一律只用 `magnetic_sim/ANSYS/main/` 底下的 config / function / 資料。

當工作涉及：
- 寫 / 改任何 `.m` 腳本，而它需要 `mt_constants`、幾何常數、或任何模型設定
- 看到腳本裡有 `addpath(... backup ...)`、讀 `backup/**/*.dat`、`backup/**/*.mat`
- 新增 `figures/paper_fig_plot/` 的繪圖腳本

**動手前先讀完此規則。**

對應 memory：`feedback_no_backup_data.md`
相關規則：`calibration-shared-structure.md`（共用校正夾已自足）、`modify-existing-files.md`、`result-read-safety.md`

---

## 🔒 核心規則

1. **禁止 `addpath` 到 `backup/`**，禁止從 `backup/` 讀 `.m`、`.mat`、`.dat`、`.db`。
2. **模型設定一律走 live config**：
   ```matlab
   % ✅ 正確
   cfg = model_config('long2016_hexapole_halfcut', 'tip40um');
   % ❌ 禁止
   addpath('...\ANSYS\backup\hexapole-long2016\analysis');  cnst = mt_constants();
   ```
   兩個分支各有一份 dispatcher：
   - `matlab/APDL/Calibration_using_FEM_modeling/function/model_config.m`
   - `matlab/Maxwell/function/model_config.m`
3. **FEM 結果**一律讀 `ANSYS_data/<topic>/...`（照 `result-read-safety.md` 三層防呆），
   不從 `backup/**/results/` 讀。
4. `backup/` 只保留**歷史歸檔**用途（追溯舊結果、對照舊做法），**不可成為活躍程式的相依**。

## 為什麼

`backup/hexapole-long2016/analysis/mt_constants.m` 是第三份平行的常數檔，
**不會跟著 live config 更新**。2026-08-08 發現 `figures/paper_fig_plot/` 有 12 支腳本吃它，
導致 CAD 實測的真實錐體幾何（`pole_cone_slope` / `pole_tip_axial`）只有校正管線看得到、
論文圖看不到 —— 同一顆模型在兩處用不同幾何。這正是要禁止的失效模式。

## ✅ 已全數清除（2026-08-08 一次做完）

**86 支腳本**的 backup 相依已全部改掉，`matlab/` 與 `figures/` 底下**零殘留**
（唯一保留：`config/long2016.../tip40um/mt_constants.m` 檔頭一行**註解**，記錄數值出處，非相依）。

改法（範本）：把 `addpath('...backup\hexapole-long2016\analysis')` 換成

```matlab
CALROOT = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling';
addpath(fullfile(CALROOT,'function'), fullfile(CALROOT,'utils'), fullfile(CALROOT,'common_path'));
cnst = model_config('<model>','<geom>');     % long2016 用 ('long2016_hexapole_halfcut','tip40um')
```

🔴 **踩過的坑**：
1. **變數名 `CAL` 會撞名**（多支腳本自己有 `CAL`，覆蓋後 `outdir` 會指錯）→ 一律用 **`CALROOT`**。
2. **只換 `addpath` 不換 `mt_constants()` 會炸**：live `function/` 沒有 `mt_constants`
   （它在 `config/<model>/<geom>/`），必須同時改成 `model_config(...)`。
3. **per-model 不可互抄**：hung/NTU 的腳本要用自己的 `model_config('hung_hexapole')` /
   `('NTU_hexapole')`，不是 long2016。（它們原本 addpath 到 long2016 的 backup 常數，本身就是錯的。）
4. python 批次替換時 Windows 路徑的 `\m` 等會被當跳脫序列 → **替換靜默失敗**。
   改用**逐行比對**、且事後一定 `grep` 驗證，不要信 replace 的回傳值。

⚠ 順帶發現（**與 backup 無關、尚未修**）：`plot/hung_hexapole/` 底下多支腳本仍指向
**已刪除的舊 per-model 樹**（`matlab/APDL/hung_hexapole/.../code/main_function`、`load_coils_actuator`）
→ 那些本來就跑不起來，要復活需另行重寫。

## 觸發片語

- 「addpath backup」/「backup 的 mt_constants」/「從 backup 讀」
- 新寫 / 改寫任何需要模型常數的腳本時

## 何時不適用

- **刻意的歷史追溯**（要重現某個舊結果、對照舊做法）—— 此時必須在回覆與註解中
  明確標示「這是 backup 歷史資料、非 live」。
- `backup/` 樹內部自己的檔案互相引用（那棵樹已凍結，不動它）。
