# 改現有檔、別開新檔（強制）

**使用者拍板（2026-07-09）**：**要改腳本/程式，一律改現有那支檔；不要動不動就複製成新的變體檔。能改現有就改，除非真的必要才創新的。**

對應 memory：`feedback_modify_existing_no_new_files.md`
相關規則：`no-structure-change-without-ask`（移/改名/刪/新建資料夾要先問）；產物落點見 `magnetic_sim/ANSYS/main/CLAUDE.md`「資料夾架構地圖」。

---

## 🔒 核心規則

當你要做「某個變體 / 微調 / 換參數 / 換 case」時：

1. **優先就地改現有腳本**，不要複製成 `_v2` / `_LowFwd` / `_HungMatch` / `_new` 之類的近乎重複檔。
2. **首選參數化**：在現有腳本頂部加一個參數（預設＝現值），需要跑不同情況時只改那個參數。同一支檔案涵蓋多情況 > 一堆變體檔。
   - APDL 例：`GLOW_R = 48.0   ! 下極 gap radial；實驗設 45`，下游引用 `GLOW_R`。
   - MATLAB 例：`VARIANT = 'gap200um_mueq'` 這種單一開關（已在用）。
3. **產物 ≠ 腳本**：`.db` / `.dat` / `.mat` / 圖 / `.iges` 是輸出，**輸出到新資料夾或新檔名是允許的**（用來保留對照、不覆蓋舊結果）。本規則只綁**腳本檔**（`.txt` deck、`.py`、`.m`、`.sh`）。
4. **真的要開新腳本**的門檻：邏輯本質不同、或硬塞進現檔會讓它太亂難讀。開之前先自問「參數化現有的行不行？」；不行才開，且**回報/commit 時說明為何非開不可**。

## ❌ 反例（本 session 踩的，別再犯）
gap/no-gap 實驗一路複製新 deck：`MT_Sim_P*_nogap.txt`、`*_samemesh.txt`、`MT_Mesh_Graded_BaseGap_HungMatch.txt` + 一堆 scratchpad 產生器 → 檔案爆炸。**正解＝參數化現有的 gap deck（一個 radial / 材料 / variant 參數），產物落新夾。**

## 觸發片語（任一即套用）
- 「改這支腳本 / 改 deck / 改 fit」
- 「做個變體 / 換位置 / 換 variant / 換參數再跑一次」
- 你正打算複製一支 `*_something.txt` 出來時 → 停，先想能不能改現有的。

## 何時不適用
- 純新功能、跟現有檔邏輯無關的全新東西（本來就該新檔）。
- 產物輸出（db/data/圖/iges）——那不是腳本。
- 使用者明確說「另存一份 / 開新的」。
