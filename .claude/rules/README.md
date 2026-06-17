# .claude/rules/ — path-scoped 編輯規則

**用途**：每個 `.md` 是一條 Claude Code 工作規則，依工作情境（觸發片語 / 涉及的目錄）自動載入；強制規則須「動手前先讀全文」。對應觸發片語總表見 `../../CLAUDE.md` 的 Quick Triggers。

**內容**（逐條一行用途）：
- `main-workspace.md` — `magnetic_sim/ANSYS/main/`（Kuo Quadrupole）工作目錄：各類產物的正規輸出位置表 + 命名規則。
- `main-workflows.md` — `main/doc/workflows/` 各 SOP 的自然語觸發對應（出 STEP / 跑 FEM / 抽場 / fit / 畫圖…）。
- `ansys-cad-alignment.md` — ANSYS 幾何尺寸必對齊 CAD（STEP/IGES = source of truth）；改幾何前必量 CAD、不一致必通報。
- `result-read-safety.md` — 讀 ANSYS 結果三層防呆（讀前回報消歧、讀後核指紋、以 RESULTS_MAP 為準），避免讀錯 case。
- `sim-cleanup.md` — 清 sim 副產物 SOP：6 項不可影響 + 2 項不可失去能力、強制 dry-run、預設 half-clean、`--full` 須明確同意。
- `apdl-editing.md` — APDL 腳本編輯規則（`[ADDED]`/`[MODIFIED]` 標記、6 coil 只差 CURR_ARRAY、保留 `D,ALL,MAG,0`）。
- `comsol-livelink.md` — COMSOL LiveLink 連線通用法（拆兩 process：獨立 server + `matlab -batch` 內 `mphstart(2036)`）。
- `fit-current-matches-sim.md` — 擬合電流必須等於 FEM 激發電流（目前 1A），不可把操作電流 0.6A 灌進 fit。
- `charge-model-source-convention.md` — 電荷模型符號慣例：每顆極激發時 B 一律從尖端射出（全 source）。
- `hung-docs.md` — 涉及 `magnetic_sim/ANSYS/backup/hung/` 時的必讀文件清單（troubleshooting、build workflow…）。
- `hexapole-build.md` — 「建 hexapole」觸發流程：依序收集 4 個參數、嚴格照 build-workflow 執行。
- `iges-model-id.md` — 貼 `IGES_converted/<topic>/*.iges` 路徑時，從路徑識別物理模型，不問「這是哪個模型」。
- `simulation-constraints.md` — APDL 幾何/材料/求解器硬約束（alpha=54.74° FIXED、tip 公式鎖定、元素型別、BC）。

**相關**：規則上層設定見 `../README.md`；觸發片語總表與 repo 總覽見 `../../CLAUDE.md`。
