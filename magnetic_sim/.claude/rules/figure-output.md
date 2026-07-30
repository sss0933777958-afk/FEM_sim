# 繪圖輸出與迭代規則

使用者 2026-06-26 拍板：「**以後圖片都輸出出來，有要改的話直接改完覆蓋，改到好為止**」。

## 🔒 規則

1. **一律實際輸出圖檔**：每次畫圖都 export 成真正的 `.png`，讓使用者**直接開檔看**。
   - **不要**只用文字描述、**不要**只丟 temp preview 然後「等定案才落地」。
2. **預設輸出到該功能組的 `figures/`**（正式位置、canonical 檔名）。例：
   `matlab/<model>/<功能組>/figures/<canonical>.png`。
3. **要改＝原地覆蓋迭代**：要修改就**直接改腳本 → 重跑 → 覆蓋同一個檔**
   （同檔名，**不另存** `_v2` / `_new` / `_final`），一輪一輪覆寫**到使用者說定案為止**。
4. **每輪輸出後必先自驗，再回報（使用者拍板 2026-07-30）**：Claude **一定要 `Read` 出來的 `.png` 目視、
   並逐項對照使用者這次的「明確要求」勾過一遍**，確認達標才給使用者看/宣稱完成。
   - 對照清單 = 使用者這輪講的每一條（視野範圍、tick、輪廓、sensor、colorbar、標題…）+ 風格要點。
   - **場圖/向量圖要額外驗「方向 / 號誌」**：箭頭指向對不對、跟參考圖/物理是否一致 —— 不可只看風格就說「方向都對」
     （2026-07-30 踩過：Maxwell P1 圖只看風格、漏驗 P2 方向就回報 → 被抓到）。
   - 發現沒達標**先自己修**（原地覆蓋迭代），或**如實講出哪裡沒到 + 原因**，不要含糊帶過。
5. **凌駕**舊 `main/CLAUDE.md` 繪圖規則 #5 的「定案前一律用 preview、不落地最終檔」——
   現在**一律落地 + 覆蓋迭代**。

## 仍守的前置（不變）

- ①畫圖前**先問功能組 + 風格選項**（見 `figure-style.md`）。
- ②**一任務一腳本、原地改**（不另開第二支腳本；一張圖對一支腳本）。
- ③場圖**真實 FEM 節點、不內插**（除非使用者明示，且須在圖說標示為內插）。

## Calibration：figures/ 佈局（2026-07-23 併入共用夾）

三模型 figures 已併入**共用夾** `Calibration_using_FEM_modeling/figures/<model>/<base>/<param>/`：
- `<model>` = long2016_hexapole_halfcut | hung_hexapole | NTU_hexapole；`<base>` = current | voltage。
- `<param>` 三子夾：
  - **`single/`** — 純 fix（USE_BIAS=false）的圖。
  - **`eighteen/`** — 純 18-param bias（USE_BIAS=true）的圖。
  - **`common/`** — 同時比較 fix+bias、或與模型無關的場/幾何/sensor 診斷圖。
- plot 腳本併入 `Calibration_using_FEM_modeling/plot/<model>/<base>/`（**不分 param**）。

plot 腳本的 `figdir` 寫 `fullfile(CAL,'figures',MODEL,BASE,<param>)`；各 `<param>/` 只放 `.png`。
（舊命名 `{single_param,eighteen_param,shared}` → 新 `{single,eighteen,common}`；舊 per-model
`<model>/Calibration_using_FEM_modeling/{current_base,voltage_base}/figures/` 已退役。`results/` 見 `results-pdf-only.md`。）

## 小提醒

- 輸出解析度要讓檔**能被直接開 / 被 `Read` 目視**（PNG 太大時可降 DPI，例如 ~150；
  不影響「實際落地覆蓋」的本意）。
- 覆蓋的是**正式 `figures/` 檔**；迭代過程不堆暫存副本。

## 觸發片語

- 「畫圖 / 出圖 / 場圖 / 畫這張」「改這張圖 / 調一下 / 換個樣式」——啟動本規則：**輸出實檔 → 覆蓋迭代**。

## 相關

- `figure-style.md`（風格 preset；畫圖前先問選項）。
- `main/CLAUDE.md`「🎨 繪圖腳本規則」（流程：先問功能組+風格 → 輸出實檔 → 覆蓋迭代到定案）。
- memory `feedback_field_quiver_style`（「先落 figures/…不要只放 temp」）、`plot_real_nodes`。
