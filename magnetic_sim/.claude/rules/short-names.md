# 檔名一律要短（強制）

**使用者拍板（2026-08-24）**：**這個專案所有檔名都不要太長，簡單就好。**
既有的短檔名慣例（`results-pdf-only.md` 把結果 PDF 由 38 字元縮到 17、
memory `short-filenames` 的 IGES/APDL 命名）從此**升格為全 repo 通則**，
不再只綁某幾個資料夾。

當工作涉及：
- **新建**任何檔案（`.m` / `.txt` deck / `.py` / `.mat` / `.pdf` / `.png` / `.step` / `.md`）
- **改**輸出檔名的組法（`sprintf` 拼 stem、`save()` 目的地、`figdir` 檔名）
- 命名新的 variant / case / 資料夾

**動手前先讀完此規則。**

對應 memory：`short-filenames`
相關規則：`results-pdf-only.md`（結果 PDF 檔名格式）、`gap-nogap-folder-convention.md`
（變體夾只准兩個名字）、`modify-existing-files.md`（別複製成 `_v2` 變體檔）。

---

## 🔒 核心規則

1. **路徑已經講過的事，檔名不要再講一次。**
   `results/zhi_peng/single/` 之下就叫 `current_R150.pdf`，
   **不要** `model_results_current_zhi_peng_single_maxwell.pdf`。
   model、base 分支、single/eighteen 由**資料夾**表達。
2. **變體後綴 3–6 字元**：`_gap`、`_split`、`_mm`、`_hp`、`_slab`。
   **不要** `_HollowProt_Plain`、`_Section_Partition` 這種完整詞組。
   單字太長就用縮寫，但縮寫要在該檔的註解 / config 裡寫明是什麼。
3. **只留「會用來區分」的欄位**。同一夾裡不會撞名的資訊就不要進檔名
   （點數 `_N528`、半徑 `_R150` 會撞才留；分支名 `maxwell` 整棵樹都一樣 → 不留）。
4. **參數值進檔名時用緊湊寫法**：`_R150`（不是 `_R150um`）、
   `_soff3mm`、`_k22_0p8340`（小數點寫 `p`）。
5. **資料夾名同此**：`R500_gap` 好，`R500_with_air_gap_remesh` 不好。
6. **變體超過 ~5 個** → 開子資料夾分層，不要繼續往檔名後面加後綴。

## 判準（寫完檔名自問）

- 這個檔名**超過 ~25 字元**了嗎？超過就檢查有沒有第 1 條的重複。
- 拿掉某一段之後，**同一個資料夾裡**還分得出來嗎？分得出來就拿掉。

## 觸發片語

- 「新建 / 另存 / 輸出成 …」「檔名叫什麼」「加個變體」
- 你正要寫一個超過 ~25 字元的檔名時 → 停，套上面的判準。

## 何時不適用

- **既有檔案不強制回溯改名**（改名會打斷既有引用與 git 歷史）；
  動到哪個檔、或那條輸出路徑本來就要改時，順手縮短。
- 外部工具規定的檔名（Maxwell 匯出、ANSYS jobname 附檔）。
- 使用者明確要求的長檔名。
