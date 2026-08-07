# CAD(STEP)→ANSYS 真實幾何匯入（跨設計通用規則）

**這份規則回答一個問題：「怎麼把一個複雜 CAD（STEP）的真實幾何可靠地匯進 ANSYS？」**
任何設計（hung / long_fei / kuo / zhang / 未來新模型）要把外部 CAD 幾何進 ANSYS 時，照這份做，不要重新踩坑。

**使用者拍板（2026-07-05）**：**複雜 CAD 一律走「匯入真實幾何」，不要再用 primitive（CYL4/CONE/BLOCK）硬拼**
（primitive 逼近複雜 B-spline 錐/T 字/半切件不收斂、來回「錯很多」；hung_hexapole 驗證通過此匯入法）。

對應 memory：`reference_cad_import_step_to_ansys.md`
相關規則：`ansys-cad-alignment.md`（CAD=source of truth、原點對齊）、`sim-cleanup.md`
相關 memory：[[iges-nurbs-measure-trap]]、[[ansys-cad-alignment]]

---

## 🔑 兩條 canonical 管線（原點不平移＝CAD 全域原點）

### 管線 A：STEP → IGES（mm，目視檢查用）
python + **OCP**（`cadquery` 內含 OpenCascade）headless：
```python
from OCP.STEPControl import STEPControl_Reader
from OCP.IGESControl import IGESControl_Writer
from OCP.Interface import Interface_Static
r = STEPControl_Reader(); r.ReadFile(step); r.TransferRoots(); shape = r.OneShape()
# 量 bbox；若被誤放大 25.4（inch 標籤）才 scale x1/25.4（實測 OCC/SpaceClaim 多半直接讀 raw=mm，不必）
Interface_Static.SetCVal_s("write.iges.unit","MM"); Interface_Static.SetIVal_s("write.iges.brep.mode",1)
w = IGESControl_Writer("MM",1); w.AddShape(shape); w.Write(out)   # unit flag=2 (MM)
```
輸出 `IGES/<model>/` + `model_check/<model>/<name>.iges`。**這份給使用者疊 CAD 目視檢查**。
範本：`magnetic_sim/ANSYS/main/apdl/hung_hexapole/geom/scripts/step_to_iges.py`。

### 管線 B：STEP → ANSYS `.db`（公尺/MKS，模型）＝三步
**1) STEP → Parasolid `.x_t`（SpaceClaim headless）**
```
"G:\ANSYS Inc\v252\scdm\SpaceClaim.exe" /Headless=True /Splash=False /ExitAfterScript=True /RunScript="<sc.py>"
```
`sc.py`（IronPython）：`DocumentOpen.Execute(step)` → `DocumentSave.Execute(out_xt)`。
⚠ **SpaceClaim headless 常不自動退出** → 跑完 `Stop-Process -Name SpaceClaim -Force`。
範本：`.../geom/scripts/sc_step_to_parasolid.py`。

**2) Parasolid `.x_t` → ANF（`ac4para.exe`，ANF 其實是 APDL 幾何檔）**
```
ac4para.exe  <in.x_t>  <out.anf>  SOLIDS  <dbname>
```
- exe：`G:\ANSYS Inc\v252\ansys\ac4\bin\para\winx64\ac4para.exe`。
- ⚠ **PATH 必須含 `G:\ANSYS Inc\v252\ansys\bin\winx64`**（ANSYS 核心 DLL），否則 exit `-1073741515`＝`0xC0000135` DLL-not-found。
- `P_SCHEMA` 環境變數要在（ANSYS 安裝已設）。

**3) MAPDL `/INPUT` 那個 ANF → `.db`**
```
/CLEAR,NOSTART
/INPUT,<name>,anf        ! ANF 是 APDL 幾何指令，直接 /INPUT 即建幾何
/PREP7
NUMCMP,ALL
SAVE
```
hung 出 **97 vol / 540 area / 680 kp**、bbox ±62.5mm、原點置中。
範本：`.../geom/import/MT_Input_ANF.txt`。輸出 db 進 `ANSYS_data/<model>/db/from_parasolid/<name>.db`。

---

## ❌ 走過但不通（此環境，別再試）
1. **MAPDL `IGESIN`**：2025 R2 只剩 `SMOOTH`（`IOPTN,IGES,ALTERNATE` 已移除，log 說「Label ALTERNAT is not recognized」）。
   對複雜 trimmed NURBS 面報「Cannot project lines to surface / Poorly defined area」→ **匯入 0 volume**。別用 IGESIN 匯複雜 CAD。
2. **MAPDL `~PARAIN`**：靠呼叫 translator wrapper bat（`ac4para252.bat`），此環境 bat 找得到但 ANF 產不出/未真執行 → 不穩。
   **改用「ac4para 直接產 ANF + `/INPUT`」繞過**（＝管線 B 步驟 2-3，最可靠）。
3. wrapper/任何 `.bat` **不可含中文**：cmd（Big5 codepage）讀 UTF-8 中文 → 亂碼、還會被當指令執行報錯。

---

## 單位（重要）
- **IGES = mm**（SolidWorks / 目視檢查用；unit flag 2）。
- **ANSYS 模型 = 公尺（MKS）**：`0.0625 m = 62.5 mm`，符合本專案「模擬用 MKS」慣例（見 backup/hung docs）。
- ANF 座標**內嵌在 `kpt` / `lcurv` / `asurf`（NURBS 控制點）**，**不可文字硬 scale**；要 mm 版模型須另議（非文字取代）。
- SpaceClaim 與 OCC 都會把「inch-標籤但實為 mm」的 STEP 讀成 raw=mm（不誤放大 25.4）——但仍**量 bbox 確認**。

---

## 路徑表（本機，換機先驗證）
| 元件 | 路徑 |
|---|---|
| SpaceClaim | `G:\ANSYS Inc\v252\scdm\SpaceClaim.exe` |
| Parasolid 轉譯器 | `G:\ANSYS Inc\v252\ansys\ac4\bin\para\winx64\ac4para.exe` |
| ANSYS 核心 DLL（加 PATH） | `G:\ANSYS Inc\v252\ansys\bin\winx64` |
| MAPDL | `G:\ANSYS Inc\v252\ansys\bin\winx64\MAPDL.exe` |
| OCP | `cadquery`（python，`from OCP...`） |

---

## 產物落點 + 清理（照 db-folder-retention / sim-cleanup + apdl 規則）
- **`.py` / `.bat` 腳本**（`step_to_iges.py` / `sc_step_to_parasolid.py` / `ac4para252.bat`）→ `apdl/<model>/geom/scripts/`。
  ⚠ **不可放 `geom/import/`**：apdl 的 `geom/import/` 依規則**只放 `.txt`**（+ 標準 README）。
- **`.txt` APDL deck**（`MT_Input_ANF.txt`）→ `apdl/<model>/geom/import/`。
- IGES：`IGES/<model>/`（原始）+ `model_check/<model>/`（mm 檢查）。
- 模型：`ANSYS_data/<model>/db/from_parasolid/<name>.db`（**db/ 只留 .db**）。
- 中間檔（`.step`/`.x_t`/`.anf`）**不留 apdl**（可重生）；別留 db/。

---

## 觸發片語（任一即啟動此規則）
- 「匯入 STEP」/「CAD 匯入 ANSYS」/「STEP 進 ANSYS」/「把這個 CAD 建進 ANSYS」
- 「IGESIN 匯不進」/「~PARAIN」/「Parasolid 匯入」/「ac4para」
- 「建真實幾何模型」/「別 primitive、直接匯入」

## 何時不適用
- **純參數化/idealized 幾何**（用 primitive 建的乾淨可控幾何、要做參數掃描）——那用 primitive 建模是對的。
- 純後處理 / 純 mesh / 純 MATLAB。
