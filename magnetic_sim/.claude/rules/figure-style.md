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
   - **字體分工（使用者拍板 2026-07-27，定案）**：**軸標題 / 刻度標題（axis label、colorbar label）一律用「標準數學字體」= LaTeX `\mathbf`（`Interpreter='latex'`，Computer Modern 粗體）**；如 `$\mathbf{|B|\;(mT)}$` / `$\mathbf{{}^{B}\hat{g}_{I}\;(mT/A)}$`。**刻度「數字」維持 Helvetica 粗體**（`FontWeight='bold'`，**不**套 `TickLabelInterpreter='latex'`）。即：標題數學字體、數字 sans-serif 粗體。⚠ 曾試把標題改 Helvetica（`\mathsf`/tex）**已被否決**——標題一律標準數學字體（且 `ᴮĝ_I` 的 hat/左上標**只有 latex 畫得出來**）。
   - **字體大小統一 = 36（使用者拍板 2026-07-28，所有 paper 圖通用）**：刻度數字 `set(ax,'FontSize',36,'FontWeight','bold')`；軸標題/colorbar 標題同 36。不要各圖各用 28/30——一律 36。（3D box 圖本就 36，2D 圖也統一到 36。）
8. **圖例每則的第一個字首字母大寫**（使用者拍板 2026-07-27）：legend 每一條的**開頭單字**要大寫（`Sampling range ≤ 150 µm`、`Mean = 0.250 mT`——不是 `sampling`/`mean`）。純符號/數學開頭（`|B|`、`ĝ_I`…）不受此限。範例：`plot_err_hist.m`。

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

### 3D 版（黑色粗體框線：**兩種 sanctioned 變體，依幾何選**）

使用者 2026-07-01 拍板、2026-07-02 補充：3D 框體**不是只有一種**。**先看三軸幾何再選框法**——選錯會被 X 掉（踩過）。

| 變體 | 何時用 | 框法 | 立方比 | 範例 |
|---|---|---|---|---|
| **A. 手動框邊** | **同尺度立方幾何**：x,y,z 都是空間、同單位、可等比（例：磁荷位置示意） | `box off` + 手動 `draw_box_edges`（省最遠角 3 邊） | `daspect([1 1 1])` | `plot_charge_positions_3d.m` |
| **B. box on** | **異質軸**：z 是跟 x,y 不同的物理量/單位/尺度（例：surf 高度/山丘圖，x,y=µm、z=gain/iso） | 直接 `box on`（MATLAB 標準 3D 框：外框 + 後方三面邊） | `pbaspect([1 1 1])`（**不能 daspect**，會把小範圍的 z 壓扁） | `plot_svd_gain_iso_3d.m`、`plot_upperP2P5_circuit_3d.m` |

**判準一句話**：**能 `daspect([1 1 1])`（三軸同單位同尺度）→ 用 A 手動框邊；不能（z 是別的量）→ 用 B `box on` + `pbaspect`**。

**變體 A 配方**（同尺度立方）：
```matlab
grid off; box off; daspect([1 1 1]);                 % daspect（不是 axis equal）保住 limits
xlim([-bh bh]); ylim([-bh bh]); zlim([-bh bh]);
view(az,el);                                         % 先設 view（draw_box_edges 用 campos）
set(gca,'FontSize',13,'FontWeight','bold','LineWidth',1.5);
set(gca,'XTick',-2:1:2,'YTick',-2:1:2,'ZTick',-2:1:2);   % 三軸「同刻度」（z 跟 x/y 一致）
draw_box_edges(bh, 3.0);                             % 手動畫框邊（省最近角 3 邊；LineWidth 隨圖，3D paper 圖用 3.0）
xlabel('x (mm)','FontWeight','bold'); ...
```
`draw_box_edges`：立方 12 邊黑，**省略「離相機最近角（`min` campos 距離）」相連 3 邊**（＝前面 3 條，含那條會**橫穿內部**穿過資料的中間直邊）→ 剩 9 邊 = 外框輪廓 + **後方框邊**，全黑等粗（標準 3D 開口箱外觀）。用 `campos` 找最近角（故 view 先設）。
⚠⚠ **2026-07-26 更正（別再犯）**：早期寫「省**最遠角**」是**錯的**——最遠角那 3 邊是**後方**框邊，使用者要它們**保留且加粗**；真正要省的是**最近角（前面）**的 3 邊（其中一條穿過球體/資料變雜線）。實作 = `cp=campos; [~,near]=min(sum((C-cp).^2,2)); 省 near 的 3 邊`。⚠ 若軸異質（z 範圍遠小於 x,y），`campos` 距離會被大範圍軸蓋掉、選錯角 → **這種就別用 A，改用 B**。

**變體 B 配方**（異質軸 surf/山丘）：
```matlab
surf(X,Y,Z,'EdgeColor','none'); shading interp; colormap(jet); caxis([min(Z(:)) max(Z(:))]);
grid off; box on;                                    % ← 標準 box on（外框 + 後方三面邊，粗體黑）
xlim(...); ylim(...); zlim([min(Z(:)) max(Z(:))]);
view(-40,30); pbaspect([1 1 1]);                     % pbaspect 立方 plot box（不用 daspect）
set(gca,'FontSize',16,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]);   % LineWidth 2 = 粗體黑框
ax=gca; ax.Toolbar.Visible='off';                    % 匯出不帶 axes 工具列
```

🔴 **踩過的坑（別重犯）**：
- **不要 `axis equal`**（撐開 limits）；A 用 `daspect([1 1 1])`、B 用 `pbaspect([1 1 1])`。
- **不要 `BoxStyle','full'`**（多畫細內邊）。
- **`box on` vs 手動框邊不可混用**：A 幾何用 `box on` 會**前粗後細/後方 3 邊沒加粗到**（使用者打槍過）；故 A 要 `box off` + 手動 `draw_box_edges`（省最近角、9 邊全粗）；**但 B（異質軸 surf）就是要 `box on`**——別硬套 A 的手動框邊。
- A 的**三軸刻度要一致**（z 常被設 0.5 間距）→ 明設 `ZTick` 跟 x/y 同；B 的 z 是別的量、各軸自然刻度即可。

沿用①：`grid off`、單位 `()`、軸標粗體。

**🔒 仰角<0(從下往上看)的框線(使用者拍板 2026-07-28)**：當 `view` 的**仰角 el<0**（相機在下、看到 box 底面），`draw_box_edges3` 改走此樣式（用 `[~,el]=view(ax)` 判斷、`el>=0` 維持上面「省最近角 3 邊」的標準開口箱）：
- **判準 = 只保留「有刻度數字」的邊，沒刻度的框線不要**（使用者原話：「我不要的是底邊沒有數字刻度的框線」）。
- **頂框(z=max) + 底框(z=min) 手動一律不畫**。有刻度數字的 x/y 底邊（例 x=2/4/6/8、y=−1/0/1）**由 MATLAB 座標軸 ruler 自己畫**（`set(ax,'LineWidth',3)` 讓 ruler 線一樣粗），不需 `draw_box_edges3` 補；沒刻度的底邊（遠端長邊、近端短邊）就此消失＝使用者要的效果。
- **垂直邊只留遠端**：省掉「最近 2 個角」相連的垂直邊（前方 stub）；z 軸 ruler 由 MATLAB 畫（含 0/−1/−2）。
- 實作 = 條件式 `draw_box_edges3`（`if el<0`：`istop||isbot → continue`（底/頂框都不手動畫）、垂直省最近 2 角；`else`：標準省最近角 3 邊）。三支 3D 腳本共用同版：`plot_p1_pole_full.m`（el=−20）、`plot_p2_pole_full.m`（el≥0 不受影響）、`plot_hexapole_sensors_3d.m`（下極 el=−25 套用、上極/合併 el>0 不受影響）。

**3D tick 放置（使用者拍板 2026-07-25；2026-07-26 補字體/框邊，A/B 皆適用）**：
- **刻度數字字體 = 36 粗體**（3D paper box 圖統一；`set(ax,'FontSize',36,'FontWeight','bold')`）。合併多面板圖同此（見 `project_long2016_paper_figures`）。
- **角落 / 三軸相交處不放 tick 與數字**：tick 取 xlim/ylim/zlim 內部、**不含 ±極值端點**（box 三軸在角落交會，端點數字會兩兩擠成一團）。做法 = tick 範圍比 limit 內縮一格（例：`±0.9` 框 → `XTick -0.5:0.5:0.5`，不要到 ±0.9）。
- **tick 數量（含原點）取奇數**：讓 `0` 一定是 tick、且對稱（例 3 個 `[-0.5,0,0.5]`、5 個 `[-500,-250,0,250,500]`）。奇數 + 不含端點 → 乾淨對稱、角落不擠。
- **刻度間距一律等距（平均）（使用者拍板 2026-07-29）**：明設 `XTick`/`YTick` 時，**各 tick 間距必須相同（等差步長）**。做法 = 中心 ± 固定 nice 步長：`s=nice((hi-lo)/N)`（nice ∈ `[1 2 2.5 3 4 5 10]×10^k`）、`ctr=round(mid/s)*s`、`tk=ctr+[-1 0 1]*s`（N=3；見 local `ticks3`）。**不要用 `round(linspace(...))`**——四捨五入會破壞等距（例變成 0/2/6 而非 0/3/6）。R-sweep 趨勢圖（`ell_gain_vs_R`/`ell_gV_vs_R`）上下兩 panel 的 y-tick 皆須等距。適用 2D/3D 各軸。
- **框邊加粗**：變體 A 用 `box off` + `draw_box_edges(bh,3.0)`（省**最近角** 3 邊、留 9 邊全粗；**別用 BoxStyle full、別省最遠角**——見上「別再犯」）。
- 範例：`figures/paper_fig_plot/plot_sphere_lattice_3d.m`（±0.9 框、三軸 `-0.5:0.5:0.5`、font 36、手動框邊省最近角）。

範例圖：A＝`…/fix_dir/figures/charge_positions_P1P2_3d.png`（view −30/−20）；
B＝`…/fix_dir/figures/svd_gain_3d.png`・`svd_iso_3d.png`（surf 山丘，view −40/30）、`…/Hall_sensor_base_fix_dir/figures/circuit_3d_*.png`。

### 方向箭頭(n+ 等)頭看不清 → 自畫箭頭、頭兩翼繞 n 自轉面向相機（方向不動）（使用者拍板 2026-07-28）

當一個方向箭頭（sensor n+、法向…）幾乎**平行視線**、箭頭頭看不清時：**方向(n+)絕不改、視角(view)也不動**，改**自己畫箭頭**（別用 `quiver3`——它的頭方位不能控、平行視線時兩翼側對相機變一條線）：
- 桿 = `plot3` 從 sc 到 `sc+n·Ln`；箭頭頭 = 兩條翼 `plot3`。
- **兩翼展開方向 `w = normalize(cross(n, 視線))`**（⊥ n 且 ⊥ 視線）→ 兩翼**面向相機**、任何視角頭都清楚（＝把箭頭頭「繞 n 自轉」到面向相機，方向 n 完全不變）。
- 翼：`b=-n; d1=cos(a)·b+sin(a)·w; d2=cos(a)·b−sin(a)·w`（`a≈24°`、翼長 `≈0.32·Ln`）。視線 `vd=[sind(az)cosd(el);-cosd(az)cosd(el);sind(el)]`。
- ⚠ 踩過的坑（別再犯）：(1) 別改箭頭**方向**（曾轉離相機/歪出 y=0 平面 → 物理錯、使用者打槍）；(2) 別轉 **view 視角**（會把整個 pole 重定向 → 不是要的）。正解只有「自畫箭頭 + 兩翼 ⊥視線」。
- 實作：local `draw_narrow(ax, sc, n, Ln, col, lw)`。範例：`plot_p1_pole_full.m` / `plot_p2_pole_full.m`（都用真實 n+）。

### 3D 磁路場箭頭（quiver3，定案 2026-07-28）

磁極磁路示意（錐內磁路 + 錐面外/尖端射出）的 3D 場箭頭，**定案畫法**：

- **內部（鐵件，`steel_ids`）+ 尖端射出扇 = raw graded 真實節點**（不內插），`voxpick3` 取每格 `|B|max`（內部 0.30mm、尖端 0.22mm）。
- **磁極表面（錐面）改用重心內插加密**（使用者要「整根均勻、後端也要密、但別像規則格那樣假」，2026-07-28）：**真實貼面節點（自然、優先）+ `scatteredInterpolant`（`linear`＝barycentric）在錐面細格（軸×方位）查詢點填補空隙 → 再用 `vsz≈0.22mm` 體素「每格留一、真實優先」**。這樣整根均勻密、graded 網格超密的尖端不成團、後端粗網格區靠內插填滿，且位置以真實節點為主（非純規則格 → 不假）。半切下極只取鋼側（`z<z_tip`）phi 0..π；全錐上極 phi 0..2π。**所有箭頭(內部+表面)一律限制在磁極軸長度 L 內**（`(P−tip)·axk ≤ L`），超過 L（根部寬端外）不畫。使用者明示**此圖表面內插不必標示**。
- **turbo 依 `|B|(mT)` 分 bin（28 bin）**上色（同 `plot_p2_charge_merged` 右圖）。此類全極示意圖使用者要**不放顏色表**。
- **箭頭長度依 `|B|` 變化，但最大長度不誇張（🔒 定案）**：用壓縮映射 `len = lmin + (lmax−lmin)·(|B|/|B|max)^0.35`，方向 = 單位 `B`×`len`（`quiver3(...,0,'AutoScale','off')` 隱含用實 len）。**`lmax` 取相對圖幅的適度上限**——**最大箭頭不可長到橫跨大半視野**（例：全極 8mm box 用 `lmin=0.15 / lmax=0.55mm`）。**不要**用單位化固定長度（全部等長、看不出大小），也**不要**線性 raw×大比例（最大箭頭爆掉）。
- 範例：`figures/paper_fig_plot/plot_p2_pole_full.m`（全極 P2，`lmin/lmax=0.15/0.55`）、`plot_p2_charge_merged.m` 的 `render_3d`（尖端 zoom，`0.014/0.060`）。

---

## 通用數值標註慣例（圖 + 結果 PDF）

使用者 2026-07-01 拍板，**圖與 emit_mat/emit_labeled_matrix 產的結果 PDF 都適用**：

1. **小指數（|指數|≤1）一律乘進值裡、不抽因子**（使用者 2026-07-10 拍板擴充）：**10^-1 / 10^0 / 10^1 都不抽 `×10^n` 因子，直接顯示原值**（如 `0.4230` 不是 `4.230×10^-1`；`7.327 mT/A` 不是 `7.327×10^0`）。**只有 |指數|≥2 才抽 10^n 因子**（如 `7.154×10^-3`、`1.23×10^5`）。實作：auto-factor helper 的判準 = 矩陣 `if abs(e)>=2` 才抽（原 `if e~=0`）、純量 `if abs(ge)<=1` 印原值（原 `if ge==0`）。範例：`gen_B_matrix.m` 的 B 矩陣（對角 0.42 直接顯示）、各 `emit_model_results.m` 的 `emit_mat/emit_e/emit_scalar_unit`。
2. **無單位 / 無因次不標**：dimensionless 的量（如 K̄_I）**不加**單位標記——不要 `[--]`、`[-]`、`[\text{--}]`；有單位才標（`[mT]`/`[A]`/`µm`/`mT/A`…，圖沿用①的 `()`）。
3. **正值不標 `+` 號**（使用者拍板 2026-07-26）：矩陣 / 數值輸出**只有負值標 `−`、正值不加 `+`**（如 `7.0445` 不是 `+7.0445`；`-1.7016` 保留）。實作：`fprintf` 格式用 `%9.4f`（不是 `%+9.4f`）。範例：`Calibration_using_FEM_modeling/function/emit_tex.m` 的 `emit_mat`/`emit_e`。
4. **矩陣用標準 bmatrix、不用欄位標籤表格**（使用者拍板 2026-07-26）：給使用者的矩陣（K̄_I / ᴮĤ / G / F…）一律 `\begin{bmatrix}` 標準矩陣（隱含 P1~P6 順序），**不要**帶 `P1..P6` 欄/列標籤的 `array{c|cccccc}` 表格。實作：`emit_tex.m` 的 `emit_mat` 已改 bmatrix。
5. **水平軸「起點 + 終點」都要標數字；縱軸首末都不標**（使用者拍板 2026-07-26/27，2D 圖通用；直方圖 / 線圖皆適用）：
   - 一張圖兩軸的端點值**只由水平軸負責標**——**x 軸起點與終點都要有數字**（起點如 `0`/線圖 `40`；終點如直方圖 `1`/線圖 `500`），**縱軸起點與末端都不標**（縱軸只留內縮 tick）。
   - **起始點不進 `XTick`**（原點在角落、不畫刻度線）：`XTick` 只放內部值（如 `[2 4 6 8]`），起始 `0` 用 `text(ax,0,-0.022*ytop,'0',...,'VerticalAlignment','top','Clipping','off')` 在原點下方**補字**（有數字、無 tick mark）。線圖第一點若本就在軸內（如 x 從 40 起），該點可直接留在 `XTick`。
   - **水平軸「終點」也一定要標數字**（使用者拍板 2026-07-27，**無例外**——直方圖也要）：水平軸的**起點與終點都要有數字**（縱軸則首末都不標）。兩種做法依情境：
     - **線圖 / 端點是真 tick**（如掃描 40~500）：起點 + 終點都放進 `XTick`（`[40 100 200 300 400 500]`），正常刻度線。
     - **直方圖 / 端點在角落**（起點在原點、終點在右緣，該處不畫 tick line）：`XTick` 只放內部值，**起點 + 終點各用 `text` 補在該端下方**（無 tick mark）：`xr=xlim(ax); text(ax,xr(1),-0.022*ytop,sprintf('%g',xr(1)),...); text(ax,xr(2),-0.022*ytop,sprintf('%g',xr(2)),...)`（`VerticalAlignment top`、`Clipping off`）。⚠ 之前「直方圖終點不標」是**錯的**，已更正。
   - 連帶：講「tick 數量」時**不含起始/終點**（「x 軸 4 個內部 tick」= `[2 4 6 8]`，另加起點 0 + 終點 1；縱軸「3 個 tick」= 首末不標的 3 個內縮值）。
   - 範例：`plot_err_hist.m`（`XTick=[0.2 0.4 0.6 0.8]` + `text` 補起點 0 + 終點 1、y 三內縮）、`plot_ell_gain_vs_R.m`（x 起點 40 + 終點 500 都在 XTick）。

範例實作：`…/fix_dir/code/function/emit_model_results.m`（K̄_I 無單位、^Bĝ_I 的 `ge==0` 分支）、`Calibration_using_FEM_modeling/function/emit_tex.m`（bmatrix + 正值不標 +）。

## 分布 / 疊圖直方圖 bin 間距（histogram spacing）

使用者 2026-07-07 拍板：分布圖 / 兩模型疊圖一律用**離散 `histogram` bars**（**不要**連續曲線 / KDE / `smoothdata` 平滑——曲線版已被否決）。

- **bin 數固定 `nb = 180`**：`edges = linspace(min(allData), max(allData), nb+1)`（bin 寬 ≈ 資料範圍/180，夠密）。兩組疊圖**共用同一 edges**。
- bars 樣式：`FaceAlpha 0.55`、`EdgeColor 'w'`、`LineWidth 0.3`；mean 用 `xline` 虛線（同色）；legend 用 histogram handles；上方留 headroom（`ylim*1.20`）讓 legend 不壓 bar。
- 範例：`plot_gain_iso_overlay.m` 的 `render_overlay`（gain 𝒞 / iso κ 疊圖）。

## |B| 場 colorbar 標準樣式（使用者拍板 2026-07-29）

**所有磁路 / |B| 場圖的 colorbar 一律照此標準**（source of truth = `plot_circuit_side.m` 的 `style_cbar`）：
- **colormap = `turbo`**；`clim([0 CLIM])`，`CLIM = ceil(max(|B|_mT)/50)*50`（進位到 50mT）。
- **標題** = LaTeX 數學粗體 `cb.Label.String = '$\mathbf{|B|\;(mT)}$'`、`cb.Label.Interpreter='latex'`、`cb.Label.FontSize=36`。
- **刻度數字** `cb.FontSize=36; cb.FontWeight='bold'`。
- **箭頭配色**：nb=28 個 turbo bin、`edges=linspace(0,CLIM,nb+1)`，逐 bin `quiver(...,'Color',cmap(k,:))`（箭頭自帶 truecolor，colormap/clim 只驅動 colorbar 圖例）。
- **貼近 panel**：colorbar 緊靠圖框、不要留大白邊（2D 圖用預設 east 即可；合併圖用下方寬度比例法）。3D 透視圖若預設 colorbar 離框太遠或標題被裁 → 手動 `set(ax,'Position',...)` 縮軸 + `cb.Position` 明確定位，把標題留在圖內。
- 直接呼叫 `style_cbar(cb, 36)`（`plot_circuit_side.m` / `plot_p1p2_poles_3d.m` 都用這支）即符合。

## Colorbar 寬度（paper 合併圖，定案 2026-07-27）

使用者拍板：**多 panel 合併 paper 圖（`plot_circuit_side.m` 的 P1|P2 側視、`plot_flux_arrows_merged.m` 的巨觀|zoom）的 colorbar 寬度，用「佔圖寬固定比例」控制、定案 `CBW_RATIO = 0.009`（= cbw / figW，細長條）**。

- **做法**：manual pixel 佈局，colorbar `Units='pixels'`、`Position=[x y cbw H]`，其中 `cbw` 由 `cbw = CBW_RATIO·base/(1−n·CBW_RATIO)` 解出（`base` = 圖寬扣掉 colorbar 的部分、`n` = 該圖 colorbar 數）。這樣不論圖寬窄，**colorbar 佔圖比例一致**（各圖縮到同寬時視覺同粗）。
- **為何用比例不用絕對 px**：不同圖 figW / 匯出解析度不同，固定 px 會讓縮放後粗細不一；比例法保證一致。
- **改寬度**：只改腳本頂端 `CBW_RATIO` 常數（`plot_circuit_side.m` 與 `plot_flux_arrows_merged.m` 共用同值，改要同步）。
- **P1/P2 共用單一 colorbar**（使用者拍板 2026-07-27，比照 `flux_arrows_merged`：一個 colorbar、shared clim = 兩 panel 全域 max；弱場 panel 顯冷色，如 flux 巨觀 panel 幾乎全藍）——**不是每 panel 各自 colorbar**（那是誤解，已改回）。合併圖靠「兩 panel 共用 box 高度 H + 同 y0」對齊上下邊框。

## 同類比較圖共用色階 / 軸尺度（shared color/axis scale for comparison）

使用者 2026-07-10 拍板：**一組要「互相比較」的同類圖**（如 filled/halfcut/tipcut 三形狀側視場圖、多變體 / 多 case 的 |B| 場圖、同量的多面板）——**必須共用同一 colorbar/`clim`（必要時也共用軸範圍）**。

- **做法**：先把「該組全部圖的資料合併」算一個共用上限（如 `CAP = max_s prctile(bsum_s, 92)`），每張都用同一 `clim([0 CAP])` + 同一 colorbar 範圍。
- **禁止各圖各自 auto-scale**：各自 `prctile`/`caxis auto` 會讓每張 colorbar 最大值不同 → **跨圖比較被誤導**（例：單極 halfcut 近表面峰值真的較高，但獨立色階會讓人誤以為「halfcut 最大值 > filled」純粹因為色階不同）。
- **兩段式實作**：Pass 1 載入全組 + 算共用 CAP；Pass 2 逐張用共用 CAP 渲染。範例：`no_fix_dir/code/plot/plot_singlepole_sideview.m`（3 形狀側視共用 CAP）。
- 何時不適用：**單張獨立圖**（無跨圖比較意圖）可自身 auto-scale；刻意要凸顯各自分布形狀的圖（另註明）。

## 觸發片語
- 「畫場圖 / 畫 contour / 畫 quiver / 出圖」——啟動本規則 → **先問要哪個風格選項**。
- 「用選項① / 粗體框圖 / 套那個風格」——直接套對應 preset。

## 之後新增 preset
往本檔 `## 選項②…` 續寫；同步更新 `README.md`、`read-rules-first.md` 清單，以及「畫圖前問選項」時列出的當前清單。

## 相關
- memory `feedback_field_quiver_style`（同款風格 + y=0 場 quiver 專屬坑：source/interp/raw/cap/前端發散）。
- memory `plot_real_nodes`（場圖預設畫真實 FEM 節點、不內插；內插須在圖說標示，除非使用者明示不標）。
- `main/CLAUDE.md`「🎨 繪圖腳本規則」、repo `…/FEM_sim/CLAUDE.md` Figure Production。
