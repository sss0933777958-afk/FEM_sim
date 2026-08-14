function plot_nmin_ratio_stem(MODE, RSEL, TOLPC, KWIN, USE_BIAS)
% plot_nmin_ratio_stem -- paper 圖：每個取樣範圍所需的點數佔全池的比例（火柴棒圖）
% =========================================================================
%   水平軸 = 取樣半徑 R [um]、縱軸 = N / N_total × 100 [%]。
%   火柴棒（stem）：細桿 + 圓頭，圓頭旁標「取樣點數 / 總格數」（旋轉 90 度避免相撞）。
%
%   MODE（第 1 引數，預設 'conv'）：
%     'conv' — [ADDED 2026-08-13] **等測度網格 + 現行三條件判準**的 N_c
%              （l_hat 穩定 ∧ g_I 穩定 ∧ K_I 符合物理，各持續 10 步）。
%              N_c 與 N_total 都直接讀 data/full_vs_conv_vs_R_maxwell.mat。
%              N_total = 該 R 球內**真實 20 um .fld 格點**數。
%              ⚠ R 從 140 起：R<=120 時 N_c 反而**超過**該 R 內的真實格點數
%                （R=80 佔 980%、R=100 佔 285%、R=120 佔 128%）——等測度點是內插的、
%                不受格點數限制，所以小 R 根本稱不上「減量」。畫進去會讓縱軸失去解析度。
%     'nmin' — 舊研究：隨機巢狀取樣的 N_min（10 組中位數、後 KWIN 點變化率 <= TOLPC%），
%              N_total = 該 R 球內 10 um 內插格的全部點數。
%              資料來自 data/nmin_vs_R_maxwell10_single_random.mat。
%
%   🔑 **不做任何擬合**，全部從快取取值。
%
%   風格 = 選項①粗體框 @ paper scale。
%   輸出 → figures/paper_fig/Section2_E/nmin_ratio_stem_maxwell{10_single,_conv}.png
% =========================================================================
    clc;
    if nargin < 1 || isempty(MODE),  MODE  = 'conv';      end
    if nargin < 3 || isempty(TOLPC), TOLPC = 0.2;         end
    if nargin < 4 || isempty(KWIN),  KWIN  = 7;           end
    if nargin < 5 || isempty(USE_BIAS), USE_BIAS = false; end   % conv 模式：false=single、true=eighteen

    here   = fileparts(fileparts(mfilename('fullpath')));
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section2_E');
    if ~exist(figdir,'dir'); mkdir(figdir); end

    if strcmpi(MODE,'conv')
        D = load(fullfile(here,'data','full_vs_conv_vs_R_maxwell.mat'));
        % [MODIFIED 2026-08-13] 起點 per-model（使用者拍板）：
        %   eighteen 從 80（全段最大才 40.6%，畫得下）；
        %   single 從 **140** —— R<=120 時 N_c 反而超過該 R 內的真實格點數
        %   （R=80 佔 980%、100 佔 285%、120 佔 128%；等測度點是內插的、不受格點數限制），
        %   放進去會把縱軸撐到 1520%，R>=140 的 38%→0.28% 變化全部壓成貼底線的一排。
        if nargin < 2 || isempty(RSEL)
            RSEL = D.R_um;
            if ~USE_BIAS, RSEL = RSEL(RSEL >= 140); end
        end
        [tf, loc] = ismember(RSEL, D.R_um);
        assert(all(tf), 'full_vs_conv 快取缺少 R = %s', num2str(RSEL(~tf)));
        fld = 'n_c1';  mstr = 'single';   lbl = 'Single parameter';
        if USE_BIAS, fld = 'n_c2';  mstr = 'eighteen';  lbl = 'Eighteen parameters'; end
        Nuse = D.(fld)(loc);   pool = D.npts_f(loc);
        outn = sprintf('nmin_ratio_stem_maxwell_conv_%s.png', mstr);
    else
        if nargin < 2 || isempty(RSEL), RSEL = 100:20:500; end
        S = load(fullfile(here,'data','nmin_vs_R_maxwell10_single_random.mat'));
        nR = numel(RSEL);   Nuse = nan(1,nR);   pool = nan(1,nR);
        for t = 1:nR
            i = find(S.RLIST == RSEL(t));   assert(~isempty(i), '快取沒有 R=%d', RSEL(t));
            n = S.nall{i};   E = S.ellall{i};   pool(t) = S.pool(i);
            nd = nan(1, size(E,1));
            for d = 1:size(E,1), nd(d) = flat_one(n, E(d,:), TOLPC, KWIN); end
            [~, k] = min(abs(n - median(nd,'omitnan')));
            Nuse(t) = n(k);
        end
        outn = 'nmin_ratio_stem_maxwell10_single.png';
        lbl  = 'Single parameter';
    end
    pct = Nuse ./ pool * 100;
    fprintf('   R    取樣點數   總格數    佔比%%\n');
    fprintf('  %3d  %8d  %7d   %6.2f\n', [RSEL(:).'; Nuse(:).'; pool(:).'; pct(:).']);
    fprintf('\n  佔比範圍 %.2f ~ %.2f %%（中位數 %.2f%%）\n', min(pct), max(pct), median(pct));

    % ---- 畫圖（火柴棒；選項①粗體框）----
    FS = 36;  LWBOX = 4.0;   col = [0.05 0.10 0.95];
    if strcmpi(MODE,'conv') && USE_BIAS, col = [0.85 0.10 0.10]; end   % 顏色 = 模型
    fig = figure('Color','w','Position',[100 100 1240 860]);
    ax  = axes(fig);  hold(ax,'on');
    hs = stem(ax, RSEL, pct, 'filled', 'Color',col, 'LineWidth',1.8, ...
              'MarkerFaceColor',col, 'MarkerEdgeColor',col, 'MarkerSize',12, ...
              'BaseValue',0, 'Clipping','off');
    hs.BaseLine.LineWidth = 2.0;   hs.BaseLine.Color = [0.35 0.35 0.35];

    box(ax,'on');  grid(ax,'off');
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LWBOX,'TickLength',[.015 .015]);
    ax.Toolbar.Visible = 'off';

    XL = [RSEL(1)-20, RSEL(end)+20];   xlim(ax, XL);      % 火柴棒兩端留一格，桿子不貼框
    set(ax,'XTick', inner_ticks(RSEL(1), RSEL(end), 3));
    % [MODIFIED 2026-08-12] 基準線（0%）**對齊下框**（使用者要求）：ylim 下緣固定 0，
    %   不留下方空白；刻度取奇數個等距、上緣留白 = 間距。
    % [MODIFIED 2026-08-13] 圓頭旁要標「取樣點數/總格數」（旋轉 90 度往上長），8% 裕度放不下
    %   → 上緣放大到 HEAD×。⚠ 刻度必須**跟著填滿放大後的軸**（原本先算 8% 版的刻度、再改
    %   上緣，結果上半截整片沒有刻度、違反規則）——故把上緣一起傳進去算刻度。
    [YL, YT] = axlim_from_zero(max(pct), 1.55);
    ylim(ax, YL);   set(ax,'YTick',YT);

    yoff = YL(1) - 0.022*diff(YL);                         % 水平軸端點只標數字、不畫 tick
    for xv = [RSEL(1) RSEL(end)]
        text(ax, xv, yoff, sprintf('%g',xv), 'HorizontalAlignment','center', ...
             'VerticalAlignment','top','FontSize',FS,'FontWeight','bold','Clipping','off');
    end

    % [RESTORED 2026-08-13] 圓頭旁標「取樣點數/總格數」（使用者要求），旋轉 90 度避免相撞
    for t = 1:numel(RSEL)
        text(ax, RSEL(t), pct(t) + 0.02*diff(YL), sprintf('%d/%d', Nuse(t), pool(t)), ...
             'Rotation',90, 'HorizontalAlignment','left', 'VerticalAlignment','middle', ...
             'FontSize',17, 'FontWeight','bold', 'Color',col, 'Clipping','off');
    end

    xlabel(ax, '$\mathbf{Sampling\;range\;(micro\;meter)}$', 'Interpreter','latex', 'FontSize',FS);
    ylab = '$\mathbf{N_{min}/N_{total}\;(\%)}$';
    if strcmpi(MODE,'conv'), ylab = '$\mathbf{N/N_{total}\;(\%)}$'; end   % [MODIFIED] N_c -> N
    ylabel(ax, ylab, 'Interpreter','latex', 'FontSize',FS);

    lg = legend(ax, hs, {lbl}, ...
                'Interpreter','tex', 'Location','northoutside', 'NumColumns',1);
    lg.FontSize = 24;  lg.FontWeight = 'bold';
    lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = 2.5;
    hold(ax,'off');  drawnow;
    % [MODIFIED 2026-08-12] 圖例寬度**自動貼合內容**：只有一則短標籤時不要拉滿座標框寬
    %   （原標準寫死 width = axp(3)，單則標籤會拉出一條又寬又空的框）。
    %   自然寬度 < 框寬 70% → 保持自然寬度並置中；否則沿用切齊框寬。
    axp = get(ax,'Position');   lgp = get(lg,'Position');
    lgw = lgp(3);   lgh = lgp(4);
    GAPN = 0.022;   newTop = 1 - lgh - GAPN - 0.006;
    axp(4) = newTop - axp(2);   set(ax,'Position',axp);
    if lgw < 0.70*axp(3)
        set(lg, 'Position', [axp(1) + (axp(3)-lgw)/2, newTop + GAPN, lgw, lgh]);
    else
        set(lg, 'Position', [axp(1), newTop + GAPN, axp(3), lgh]);
    end

    out = fullfile(figdir, outn);
    exportgraphics(fig, out, 'Resolution', 200);
    fprintf('\nwrote %s\n', out);
end

% ============================================================================
function k = flat_one(nlist, ell, TOLPC, KWIN)
    dch = [NaN, abs(diff(ell))./ell(2:end)*100];
    k = NaN;
    for a = 1:numel(nlist)-KWIN
        if all(dch(a+1 : a+KWIN) <= TOLPC), k = nlist(a);  return; end
    end
end

% ============================================================================
function tk = inner_ticks(lo, hi, n)
    cand = [1 2 2.5 3 4 5 10];
    s = (hi-lo)/(n+1);   k = floor(log10(s));
    [~, i] = min(abs(cand*10^k - s));   s = cand(i)*10^k;
    ctr = round(((lo+hi)/2)/s)*s;
    tk  = ctr + (-(n-1)/2 : (n-1)/2)*s;
    tk  = tk(tk > lo & tk < hi);
end

% ============================================================================
function [lim, tk] = axlim_from_zero(maxv, HEAD)
% [MODIFIED 2026-08-13] 第 2 引數改成「上緣裕度倍率」（原本是刻度數候選、實際沒在用）。
%   HEAD=1.08 → 照 figure-style「自 0 起的軸上緣只留 8%」；需要放標註時傳更大的倍率，
%   刻度會跟著填滿放大後的軸（奇數個、等距、最大刻度 < 上緣）。
    if nargin < 2 || isempty(HEAD), HEAD = 1.08; end
    cand = [1 2 2.5 3 4 5 10];
    x = maxv/4;   k = floor(log10(x));
    s = cand(find(cand*10^k >= x, 1)) * 10^k;      % 間距由**資料範圍**決定，與上緣無關
    top = HEAD * maxv;
    n = floor((top - 1e-12)/s);
    if mod(n,2) == 0, n = n - 1; end
    n = max(n, 1);
    lim = [0, top];   tk = (1:n)*s;
end

% ============================================================================
function [lim, tk] = axlim_from_zero_old(maxv, nlist) %#ok<DEFNU>
% 下緣固定 0（基準線貼下框）；**上緣與刻度間距解耦**，避免過多留白。
%   [MODIFIED 2026-08-12] 原本強制 上緣 = (n+1)*s（留白 = 間距），資料 39.2 會被撐到
%   上緣 60（填充率僅 65%、上方空掉三分之一）。改為：
%     ① 上緣 = **最小的 nice 值 >= 1.05 × 資料最大值**（1.05 = 資料不貼框的裕度）
%     ② 刻度 = (1:n)*s，n 取奇數、s 取 nice，且最大刻度 < 上緣
%     ③ 在可行組合取「刻度數最多」者（刻度較密、讀數方便）
%   本例：上緣 50、刻度 10/20/30 → 填充率 78%（原 65%）。
    cand = [1 2 2.5 3 4 5 10];
    % ① 刻度間距 s：照資料範圍取 nice（維持一般慣例，約 3~4 格涵蓋資料）
    x = maxv/4;   k = floor(log10(x));
    s = cand(find(cand*10^k >= x, 1)) * 10^k;
    % ② 上緣**只留 8% 裕度**（壓掉多餘空白；limit 值不顯示，不必是整數）
    top = 1.08 * maxv;
    % ③ 刻度 = (1:n)*s，n 取「不超出上緣」的**最大奇數**
    n = floor((top - 1e-12)/s);
    if mod(n,2) == 0, n = n - 1; end
    n = max(n, 1);
    lim = [0, top];   tk = (1:n)*s;
end

% ============================================================================
function [lim, tk] = axlim_auto(lo, hi, nlist)
% 奇數個等距 tick、兩端留白 = tick 間距；資料自 0 起時優先讓最低 tick = 0。
    cand = [1 2 2.5 3 4 5 10];
    mid  = (lo+hi)/2;   rng_ = max(hi-lo, realmin);
    best = {};   bestSpan = inf;
    for n = nlist
        k0 = floor(log10(rng_/(n+1)));
        for k = k0:(k0+4)
            hit = false;
            for c = cand
                s = c*10^k;
                T = { round(mid/s)*s + (-(n-1)/2 : (n-1)/2)*s };
                if lo >= 0 && lo < s, T = [{(0:n-1)*s}, T]; end %#ok<AGROW>
                for it = 1:numel(T)
                    t = T{it};   L = [t(1)-s, t(end)+s];   clr = 0.15*s;
                    if lo >= L(1)+clr && hi <= L(2)-clr
                        if (n+1)*s < bestSpan, bestSpan = (n+1)*s;  best = {L, t}; end
                        hit = true;  break;
                    end
                end
                if hit, break; end
            end
            if hit, break; end
        end
    end
    if isempty(best)
        n = nlist(1);  s = rng_/(n+1);
        t = mid + (-(n-1)/2 : (n-1)/2)*s;   best = {[t(1)-s, t(end)+s], t};
    end
    lim = best{1};   tk = best{2};
end
