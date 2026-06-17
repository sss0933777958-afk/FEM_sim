# long2016_hexapole_halfcut/common/ — ★ 路徑 resolver（唯一來源）

**用途**：所有分析腳本的**路徑解析單一來源（single source of truth）**。把「FEM 結果在哪、MATLAB 成果在哪」集中在兩支 function，腳本不再硬寫絕對路徑（drive letter）；資料夾改名／整個 repo 搬移時，只改這兩支裡各一行即可。

**內容**（兩支 resolver，皆相對自身 `.m` 位置定位，不寫死磁碟路徑）：
- `ansys_path.m` — 解析 **FEM 輸入/結果根**（`.dat` 場 / `.db` 模型 / `.cdb`）。
- `matlab_path.m` — 解析 **MATLAB 分析成果根**（`.mat` / `.csv` / `.npz`）。

**怎麼用**（兩支簽章相同：`(model, varargin)`，`varargin` 一路 `fullfile` 串到底）：
```matlab
addpath('.../main/matlab/long2016_hexapole_halfcut/common');   % 先把 common 加入路徑
model = 'long2016_hexapole_halfcut';

% 讀 FEM 結果（ANSYS_data/<model>/...）
ansys_path()                              % .../main/ANSYS_data
ansys_path(model)                         % .../main/ANSYS_data/<model>
ansys_path(model,'coil1','standard')      % .../main/ANSYS_data/<model>/coil1/standard
d = import_ansys_data(ansys_path(model,'coil1','standard'), 'all', 'coil1');

% 讀/寫 MATLAB 成果（MATLAB_data/<model>/<功能>/...）
matlab_path(model)                        % .../main/MATLAB_data/<model>
matlab_path(model,'charge_fit')           % .../<model>/charge_fit
matlab_path(model,'bs_matrix','B_bar.mat')% .../<model>/bs_matrix/B_bar.mat
```

**自身定位機制**：兩支都用 `here = fileparts(mfilename('fullpath'))` 取得自己所在的 `.../common`，再 `fileparts ×3` 回到 `.../main`，最後接 `ANSYS_data` / `MATLAB_data`。資料夾名字各寫在檔內**唯一一行常數**（`DATA_DIRNAME` / `MAT_DIRNAME`）——要改名只動那一行。

**資料來源 / 流向**：本層不算任何東西，只回傳路徑字串。`ansys_path` → 指向 `ANSYS_data`（FEM `.dat/.db`，讀）；`matlab_path` → 指向 `MATLAB_data`（分析 `.mat`，讀寫）。功能組腳本（`bias_fit` / `bs_matrix` / `sensor_d`）一律透過這兩支取路徑。

**命名 / 慣例**：
- ⚠ **`ansys_path.m` / `matlab_path.m` 絕不可隨意改動**——全功能組腳本都依賴其簽章與回傳結構；唯一允許的修改是「資料夾改名時改那唯一一行常數」。
- 不要在腳本裡硬寫 `G:\...\ANSYS_data\<model>` 絕對路徑，一律走 resolver。

**相關**：見上層 `../README.md`（matlab/ 功能組 schema）、`../../../CLAUDE.md`（main/ 工作目錄全域規則、Resolver 段）。
