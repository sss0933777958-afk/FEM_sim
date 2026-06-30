# 繪圖風格選項（figure style presets）

`main/` 圖表的**視覺風格 preset 目錄**。每個 preset 是一套可整批套用的軸/字/色階樣式。
**這些是「選項」，不是強制預設**——畫圖流程仍照 `main/CLAUDE.md` 繪圖腳本規則（先確認功能組 → **問風格選項** → **輸出實檔 → 覆蓋迭代到定案**，見 `figure-output.md`）。

## 🔒 強制：畫圖前必先問使用者要哪個風格選項

使用者 2026-06-26 拍板：**任何「要產圖」的任務，動手前必須先問使用者「要用哪個風格選項？」**，把當前清單列出來給挑（目前只有：**①粗體框圖**）。

- 不可自己預設、不可憑記憶猜。
- **問的時機**＝每張圖任務的**開始**（跟「先確認功能組」一起問一次）；同一張圖在 preview 階段來回微調**不必每次重問**，除非使用者要換風格。
- preset 變多時，列出當前所有選項讓使用者選。

---

## 選項①：粗體框圖（bold-framed）

使用者 2026-06-26 拍板的 7 條（場圖 / contour / 一般 2D 圖通用）：

1. **字體加大加粗**：所有文字（軸標題、tick 數字、colorbar 標籤與數字）放大 + 粗體（`FontWeight bold`）。
2. **軸線 + tick 加大加粗**：spine / box line 加粗；tick mark 加長加粗。
3. **外框框出**：`box on`——四邊 spine 全顯示。
4. **移除背景網格**：`grid off`。
5. **tick 數量減半（x、y 兩軸都要）**：取現有 tick 每隔一個。
6. **右邊 colorbar 同樣處理**：字加大加粗、colorbar tick 也減半。
7. **單位用括號 `()`**：`x (mm)` / `z (mm)` / `|B| (T)`（**不是** `[]`）；座標一律用 **mm**。

（補充慣例，沿用既有圖：通常**無標題**、圖上**不標「內插」**字樣——見交叉連結。）

### 具體實作參數（可調，數值沿用既有定案 memory）

**MATLAB**
```matlab
set(ax,'FontSize',16,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]);
box on; grid off;
xt = get(ax,'XTick'); set(ax,'XTick',xt(1:2:end));   % x tick 減半
yt = get(ax,'YTick'); set(ax,'YTick',yt(1:2:end));   % y tick 減半
xlabel('x (mm)'); ylabel('z (mm)');                  % 單位用 ()
cb = colorbar; cb.FontSize = 16; cb.FontWeight = 'bold';
cb.Ticks = cb.Ticks(1:2:end);                         % colorbar tick 減半
cb.Label.String = '|B| (T)'; cb.Label.FontWeight = 'bold';
```

**matplotlib（對應）**
```python
ax.tick_params(width=2, length=7, labelsize=16)
for s in ax.spines.values():
    s.set_visible(True); s.set_linewidth(2)          # 四邊 box + 加粗
ax.grid(False)
ax.set_xticks(ax.get_xticks()[::2])                  # x tick 減半
ax.set_yticks(ax.get_yticks()[::2])                  # y tick 減半
ax.set_xlabel('x (mm)', fontweight='bold', fontsize=16)
ax.set_ylabel('z (mm)', fontweight='bold', fontsize=16)
for lbl in ax.get_xticklabels()+ax.get_yticklabels(): lbl.set_fontweight('bold')
cbar.ax.tick_params(labelsize=16)
cbar.set_ticks(cbar.get_ticks()[::2])                # colorbar tick 減半
cbar.set_label('|B| (T)', fontweight='bold', fontsize=16)
```

匯出建議：DPI 150、figure ~1180px 寬（PNG < 2000px 才能被 Read 目視）。

---

## 觸發片語
- 「畫場圖 / 畫 contour / 畫 quiver / 出圖」——啟動本規則 → **先問要哪個風格選項**。
- 「用選項① / 粗體框圖 / 套那個風格」——直接套對應 preset。

## 之後新增 preset
往本檔 `## 選項②…` 續寫；同步更新 `README.md`、`read-rules-first.md` 清單，以及「畫圖前問選項」時列出的當前清單。

## 相關
- memory `feedback_field_quiver_style`（同款風格 + y=0 場 quiver 專屬坑：source/interp/raw/cap/前端發散）。
- memory `plot_real_nodes`（場圖預設畫真實 FEM 節點、不內插；內插須在圖說標示，除非使用者明示不標）。
- `main/CLAUDE.md`「🎨 繪圖腳本規則」、repo `…/FEM_sim/CLAUDE.md` Figure Production。
