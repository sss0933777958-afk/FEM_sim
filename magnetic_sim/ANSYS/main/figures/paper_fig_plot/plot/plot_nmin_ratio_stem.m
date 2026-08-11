function plot_nmin_ratio_stem(RSEL, TOLPC, KWIN)
% plot_nmin_ratio_stem -- paper 圖：每個取樣範圍所需的點數佔全池的比例（火柴棒圖）
% =========================================================================
%   水平軸 = 取樣半徑 R [um]、縱軸 = N_min / pool × 100 [%]。
%   火柴棒（stem）：細桿 + 圓頭，圓頭旁標「N_min/pool」的實際數量（旋轉 90 度避免相撞）。
%
%   N_min(R)：對每組隨機排列，找第一個 n 使**其後 KWIN 點**的相對變化率都 <= TOLPC%；
%             10 組取中位數（判準與 plot_ell_nmin_vs_full 完全一致）。
%   pool(R) ：該 R 球內 10 um 內插格的全部點數。
%
%   🔑 **不做任何擬合**，全部從 data/nmin_vs_R_maxwell10_single_random.mat 取值。
%
%   風格 = 選項①粗體框 @ paper scale。
%   輸出 → figures/paper_fig/Section2_E/nmin_ratio_stem_maxwell10_single.png
% =========================================================================
    clc;
    if nargin < 1 || isempty(RSEL),  RSEL  = 100:20:500; end
    if nargin < 2 || isempty(TOLPC), TOLPC = 0.2;        end
    if nargin < 3 || isempty(KWIN),  KWIN  = 7;          end

    here   = fileparts(fileparts(mfilename('fullpath')));
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section2_E');
    if ~exist(figdir,'dir'); mkdir(figdir); end
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
    pct = Nuse ./ pool * 100;
    fprintf('   R    N_min     pool     佔比%%\n');
    fprintf('  %3d  %6d  %7d   %6.2f\n', [RSEL; Nuse; pool; pct]);
    fprintf('\n  佔比範圍 %.2f ~ %.2f %%（中位數 %.2f%%）\n', min(pct), max(pct), median(pct));

    % ---- 畫圖（火柴棒；選項①粗體框）----
    FS = 36;  LWBOX = 4.0;   col = [0.05 0.10 0.95];
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
    [YL, YT] = axlim_from_zero(max(pct), [5 3]);
    ylim(ax, YL);   set(ax,'YTick',YT);

    yoff = YL(1) - 0.022*diff(YL);                         % 水平軸端點只標數字、不畫 tick
    for xv = [RSEL(1) RSEL(end)]
        text(ax, xv, yoff, sprintf('%g',xv), 'HorizontalAlignment','center', ...
             'VerticalAlignment','top','FontSize',FS,'FontWeight','bold','Clipping','off');
    end

    % [REMOVED 2026-08-12] 圓頭旁的「N_min/pool」數字標註已拿掉（使用者要求）。

    xlabel(ax, '$\mathbf{Sampling\;range\;(micro\;meter)}$', 'Interpreter','latex', 'FontSize',FS);
    ylabel(ax, '$\mathbf{N_{min}/N_{total}\;(\%)}$',         'Interpreter','latex', 'FontSize',FS);

    lg = legend(ax, hs, {'Single parameter'}, ...
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

    out = fullfile(figdir, 'nmin_ratio_stem_maxwell10_single.png');
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
function [lim, tk] = axlim_from_zero(maxv, nlist)
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
