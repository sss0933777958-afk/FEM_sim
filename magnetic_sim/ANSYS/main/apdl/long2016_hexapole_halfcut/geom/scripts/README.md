# apdl/long2016_hexapole_halfcut/geom/scripts/ — 產生 geom 匯出 APDL 的 .py 輔助

**用途**：以 Python 自動「寫出」`geom/export/` 的 APDL 幾何匯出腳本（mm + metre 兩版），避免手刻同步出錯。`.py` 只負責產生 .txt，不跑 ANSYS。

**內容**：
- `_generate_geom_export_gap.py` — 從 halfcut baseline 產生 gap200um 的兩支 export 腳本（mm → `IGES_converted` flag 6→2、metre → `IGES`）；套用 protrusion CYL4 拆兩段 split，與 `sim/scripts/_generate_halfcut_sims_gap.py` 一致。

**資料來源 / 流向**：讀 baseline geom 樣板 → 寫出 `geom/export/MT_Geom_Export*_gap200um.txt`；實際 ANSYS run + units-flag patch 由對應的 `run_geom_export_gap200um.ps1` wrapper 完成（不在此層）。

**命名 / 慣例**：產生器以底線 `_` 開頭（`_generate_*.py`）標示為輔助腳本；產生的 .txt 才是 FEM input；改 .py 後須重跑產生對應 .txt，勿手改產物與產生器脫鉤。

**相關**：見 `../README.md`、`../export/README.md`、`../../sim/scripts/README.md`、`doc/workflows/cad-export.md`。
