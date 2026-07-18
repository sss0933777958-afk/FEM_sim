# apdl/hung_hexapole/geom/scripts/ — CAD 匯入用的 .py / .bat 輔助腳本

**用途**：hung_hexapole 從 CAD（`CAD_model/hung_hexapole/STEP/Full_Assembly.STEP`）匯入真實幾何到 ANSYS 的非-APDL 輔助腳本。依 apdl 規則，`geom/import/` 只放 `.txt`（APDL decks），故 `.py`/`.bat` 一律歸此 `scripts/`。

**內容**：
- `build_ws_variants.sh` — **工作空間半徑 (R_sphere) 變體產生器**：從 canonical `../export/MT_Geom_FullAssembly.txt`（R=0.5mm）sed 出指定 R 的暫存 deck → 跑 MAPDL 建幾何 → 存 geom `.db`（`db/geom/geom_full_R<tag>/`）+ IGES 中繼(flag 6→1) → 呼叫 `iges_to_step.py` 轉 **mm STEP → `model_check/`**（交付檢查一律 STEP，per deliver-step-for-check）。用法：`bash build_ws_variants.sh 0.7 0.3`。canonical deck 不改；tip 依 magic-angle 線性隨 R 縮放、外結構固定、tip-to-tip=2R。
- `iges_to_step.py` — OCC(OCP)：flag=1 的 ANSYS IGES → 明確 mm STEP（`write.step.unit=MM`），讀回驗 bbox + 6 tip 距原點=R。用法：`python iges_to_step.py <in.iges> <out.step> [R_mm]`。
- `step_to_iges.py` — python/**OCP**：STEP → IGES（mm, unit flag 2），供目視檢查（輸出 `IGES/` + `model_check/hung_hexapole/`）。
- `sc_step_to_parasolid.py` — **SpaceClaim headless** IronPython：STEP → Parasolid `.x_t`（給 ac4para 轉 ANF）。
- `ac4para252.bat` — MAPDL `~PARAIN` 用不上（此環境不穩），改由此 wrapper 直接呼叫 `ac4para.exe`：`.x_t` → `.anf`（APDL 幾何檔）。**bat 純 ASCII、內含設 PATH（含 `ansys\bin\winx64` 核心 DLL）**。

**流向（canonical，全文見 `.claude/rules/cad-import-ansys.md`）**：
STEP →(`step_to_iges.py`)→ IGES(檢查)；STEP →(`sc_step_to_parasolid.py`)→ `.x_t` →(`ac4para252.bat`/`ac4para.exe`)→ `.anf` →(MAPDL `/INPUT` = `../import/MT_Input_ANF.txt`)→ `.db`。

**慣例**：中間檔（`.step`/`.x_t`/`.anf`）不留 apdl；`.db` 進 `ANSYS_data/hung_hexapole/db/`。

**相關**：`../README.md`、`../import/README.md`、`.claude/rules/cad-import-ansys.md`、memory `cad-import-step-to-ansys`。
