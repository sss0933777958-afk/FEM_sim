# apdl/long2016_hexapole_halfcut/sim/scripts/ — 產生 sim APDL 的 .py 輔助

**用途**：以 Python 自動「寫出」`sim/` 各變體的 6 支 APDL 求解腳本，確保 6 顆極只差 `CURR_ARRAY`、其餘同步、CWD 正確。`.py` 只產生 .txt，不跑 ANSYS。

**內容**：
- `_build_from_long2016_verbatim.py` — 直接從 Long2016 source（`hexapole-long2016/apdl/MT_..._Coil1.txt`，已含 VADD 合併的正確 topology）verbatim 寫出 6 支 baseline，每支 3 處最小改動（CWD、CURR_ARRAY、…）。修掉 full 版被 strip VADD 造成的 30-70% 弱磁通。
- `_build_gap200um_mu_eq.py` — 幾何不動，加 `MAT_PROT` μ_r=31、mesh 後把 protrusion 元素 reassign，產生 gap200um_mueq 6 支。
- `_generate_halfcut_sims.py` — 從 full-hexapole 樣板產生 halfcut（每下極 VROTAT 後加 VSBV 切半錐）。

**資料來源 / 流向**：讀 `hexapole-long2016` / `long2016_hexapole_full` 樣板 → 寫出 `sim/{baseline,gap200um_mueq}/MT_Sim_P*.txt`；CWD 指向 `ANSYS_data/long2016_hexapole_halfcut/coilN/`。

**命名 / 慣例**：產生器以底線 `_` 開頭標示輔助；產物 .txt 才是 FEM input；6 支只差 `CURR_ARRAY`；改 .py 後重跑產生 .txt，勿手改產物與產生器脫鉤。

**相關**：見 `../README.md`、`../baseline/README.md`、`../../geom/scripts/README.md`、`.claude/rules/apdl-editing.md`。
