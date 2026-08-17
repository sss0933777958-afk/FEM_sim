# hung_hexapole — ANSYS_data

Hung hexapole 的 ANSYS 幾何/模型資料。來源 CAD：`../../CAD_model/hung_hexapole/STEP/Full_Assembly.STEP`
（mm 尺寸；檔案 inch 標籤視為誤標）。

## data/（FEM 場 `.dat`，2026-07-14 結構＝`<variant>/coil<N>`）
`data/<variant>/coil<N>/coil<N>_{coord,bfield}_{all,wp}.dat`。變體：`gap_200um`（μ_eff 氣隙）+ `no_gap`（無氣隙對照），
各 coil1–6 自激 1A 解。**指紋/物理意義查 `RESULTS_MAP.md`**（讀取前必查，per `../../../.claude/rules/result-read-safety.md`）。
> 2026-07-14 由 `coil<N>/<variant>` 翻轉為 `<variant>/coil<N>`（`ANSYS_data/` 全包統一）。

> 交付檔白名單（照 `.claude/rules/ansys-db-cleanup.md`）：`db/<case>/` 只留 `.db`（+主 `.rmg`）；
> log/`.out`/`.err`/`.iges` 一律不留（IGES 交付檔在 `../../model_check/hung_hexapole/`）。

## 目前幾何來源：真實 CAD 匯入（✅ 已成功進 ANSYS）
- **ANSYS 模型（真實幾何）**：`db/from_parasolid/hung_hexapole_full.db`
  — **97 volumes / 540 areas / 680 keypoints**（與 STEP 完全一致）、bbox ±62.5mm、原點=CAD 原點。
  **單位=公尺（MKS）**（0.0625m=62.5mm；符合本專案「模擬用 MKS」慣例）。
- **IGES（真實幾何、供目視檢查、mm）**：`../../model_check/hung_hexapole/hung_hexapole_full.iges`
  （OCP 從 STEP 轉出，mm、97 solid、540 face、原點=CAD 原點；疊到 Full_Assembly 應完全吻合）。

### 匯入管線（可重現）
STEP →(OCP `step_to_iges.py`)→ IGES(mm, 檢查用)；
STEP →(SpaceClaim `sc_step_to_parasolid.py`)→ `.x_t`(Parasolid, 公尺) →(`ac4para.exe`)→ `.anf`(APDL 幾何) →(MAPDL `/INPUT` = `MT_Input_ANF.txt`)→ `.db`。
- ⚠ 走過但**不通**：MAPDL `IGESIN`（2025 R2 只 SMOOTH，對 9 個複雜 NURBS 面失敗）、`~PARAIN`（呼叫 translator bat 在此環境不穩）。
  → canonical = **ac4para 直接產 ANF + `/INPUT`**（見 `../../apdl/hung_hexapole/geom/import/README.md`）。

## db/
- `db/from_parasolid/hung_hexapole_full.db` — ✅ **真實幾何 ANSYS 模型**（97 vol，公尺/MKS）。
- `db/geom_schematic/hunghex.db` — 舊 primitive 示意重建（CYL4/CONE/BLOCK），**已被真實匯入取代、僅留參考**。

## 相關路徑
- APDL：`../../apdl/hung_hexapole/geom/{export(primitive schematic), import(STEP→IGES + IGESIN)}/`
- IGES：`../../IGES/hung_hexapole/`（原始）、`../../model_check/hung_hexapole/`（交付/檢查）
