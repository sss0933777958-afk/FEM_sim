function plot_ring_azimuth_profile(CASE)
% plot_ring_azimuth_profile -- 磁極表面「環」上 ‖b‖ 的剖面（Maxwell）
% =========================================================================
%   環 = 垂直極軸的平面在 s0 = 4.4832mm（沿錐面斜距 4.572mm）截磁極，沿軸給寬 0.30mm，
%   切成 0.1 x 0.1mm 小塊；每塊格心沿外法線位移 δ = 0.15mm 後內插取 ‖b‖。
%
%   CASE（對應 data/ring_<CASE>_maxwell.mat，快取自描述）：
%     'P2cone'  P2 完整錐 -> O 形整圈，橫軸 = 方位角（0 = P2 sensor 方位）
%     'P1flat'  P1 削平面弦段，       橫軸 = 面上位置 mm（0 = 中心 = P1 flat sensor 方位）
%     'P1cone'  P1 下半錐面，         橫軸 = 方位角（0 = 正下方 = P1 cone sensor 方位；±90 = 切邊）
%   P1 的 D 形環拆成兩段畫：削平面與錐面是**不同的面**（法線差近 90°），
%   平滑不跨切邊做，否則會把切邊的峰抹掉。
%
%   取樣（詳見快取 .note）：
%     δ = 0.15mm 固定，每格 1 個內插點（格內撒點無效——誤差為系統性，實測 100 點/格
%     只降噪 11%；跨格平均才有效，近 √N）。平滑 = 5 格移動平均。
%     δ=0.15 是 0.1mm 規則格能支撐**真實內插**的最近距離（δ≤0.05 被邊界四面體壓平）。
%     空氣側篩選須留 0.02mm 餘裕 + ‖B‖<40mT，否則會混入鋼件格點。
%   ⚠ 錐面半徑 **R(s) = s·tan(beta)，無 R_TIP 偏移**——已用原始格點分類實測驗證
%     （R = 0.0097 + 0.1985·s，三站位誤差 < 0.005mm）。
%   ⚠ 這是**內插**（Maxwell 匯出 0.1mm 規則格）。
%   ⚠ 高導磁界面上磁通近乎垂直射出（μr=280）→ ‖b‖ ≈ b·n̂⁺。
%
%   風格①粗體框圖；橫軸依 figure-style 2026-08-06：首末點貼框、端點只標數字不畫 tick。
% =========================================================================
    clc;
    if nargin < 1 || isempty(CASE), CASE = 'P2cone'; end

    here   = fileparts(fileparts(mfilename('fullpath')));
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section3_A');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    cachef = fullfile(here, 'data', sprintf('ring_%s_maxwell.mat', CASE));
    assert(exist(cachef,'file')==2, 'cache 不在:%s', cachef);
    C = load(cachef);

    x = C.x;  raw = C.mag;  sm = C.mag_sm5;
    if isfield(C,'periodic') && C.periodic          % 週期環：補一點讓曲線首尾接上框邊
        x = [x; -x(1)];  raw = [raw; raw(1)];  sm = [sm; sm(1)];
    end

    fprintf('%s   δ = %.2f mm\n  %s\n', C.seg, C.delta, C.x0note);
    fprintf('  環段平均 %.3f mT | 全距 原始 %.2f%% -> 5格平滑 %.2f%%\n', mean(C.mag), ...
            (max(C.mag)-min(C.mag))/mean(C.mag)*100, (max(sm)-min(sm))/mean(C.mag)*100);
    [~,i0] = min(abs(C.x));
    fprintf('  x=0 -> %.3f mT | 峰 %.3f @x %+.3f | 谷 %.3f @x %+.3f\n', ...
            C.mag_sm5(i0), max(sm), x(find(sm==max(sm),1)), min(sm), x(find(sm==min(sm),1)));

    FS  = 36;
    COL = [0.85 0.10 0.10];  if isfield(C,'color'), COL = C.color; end
    fig = figure('Color','w','Position',[100 100 1250 900]);  ax = axes(fig);  hold(ax,'on');
    plot(ax, x, sm, '-o', 'Color',COL, 'LineWidth',3.0, ...
         'MarkerSize',8, 'MarkerFaceColor',COL, 'MarkerEdgeColor','k', 'Clipping','off');

    box(ax,'on');  grid(ax,'off');
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',2.5,'TickLength',[.015 .015],'TickDir','out');

    xe = [min(x) max(x)];
    xlim(ax, xe);
    if isfield(C,'xtick')                                  % 複合單位軸：自訂刻度 + 標籤
        set(ax,'XTick', C.xtick, 'XTickLabel', C.xticklab);
    else
        set(ax,'XTick', xticks_odd3(xe));                  % 3 個(奇數)等距、端點不進 tick
    end
    [yr, yt] = ylim_pick(sm);
    ylim(ax, yr);  set(ax,'YTick', yt);

    for xc = C.edge(:).'                                   % 切邊（若有）= 單位切換點
        plot(ax, [xc xc], yr, '--', 'Color',[0.45 0.45 0.45], 'LineWidth',2);
    end

    yoff = yr(1) - 0.022*diff(yr);                         % 端點只標數字、不畫 tick
    if isfield(C,'endlab'), elab = C.endlab; else, elab = {num2str(xe(1),'%.2f'), num2str(xe(2),'%.2f')}; end
    for k = 1:2
        text(ax, xe(k), yoff, elab{k}, 'HorizontalAlignment','center', ...
             'VerticalAlignment','top','FontSize',FS,'FontWeight','bold','Clipping','off');
    end

    xlabel(ax, C.xlab, 'Interpreter','latex', 'FontSize',FS);
    ylabel(ax, '$\mathbf{\|b\|\;(mT)}$', 'Interpreter','latex', 'FontSize',FS);
    ax.Toolbar.Visible = 'off';  hold(ax,'off');

    out = fullfile(figdir, [C.fname '.png']);
    exportgraphics(fig, out, 'Resolution', 150);
    fprintf('wrote %s\n', out);
end

% ============================================================================
function t = xticks_odd3(xe)
% 3 個奇數等距內部 tick（不含端點）：以中點 ± nice 步長
    s = nice_step(diff(xe)/4);   c = round(mean(xe)/s)*s;
    t = c + (-1:1)*s;
    while t(1) <= xe(1) || t(end) >= xe(2), s = s/2; t = c + (-1:1)*s; end
end

% ============================================================================
function [yr, yt] = ylim_pick(v)
% 奇數個等距 nice tick + 兩端留白 = 間距；先試 5 個，太空曠才退回 3 個
    lo = min(v);  hi = max(v);
    for n = [5 3]
        s = nice_step((hi-lo)/(n-1));
        while (n-1)*s < (hi-lo), s = next_step(s); end
        ctr = round((lo+hi)/2 / s)*s;
        yt  = ctr + (-(n-1)/2:(n-1)/2)*s;
        yr  = [yt(1)-s, yt(end)+s];
        if yr(1) <= lo && yr(2) >= hi && (hi-lo)/diff(yr) > 0.35, return; end
    end
end

function s = nice_step(x)
    k = floor(log10(max(x, realmin)));  m = x/10^k;
    cand = [1 2 2.5 5 10];  [~,i] = min(abs(cand - m));  s = cand(i)*10^k;
end

function s2 = next_step(s)
    k = floor(log10(s));  m = round(s/10^k, 4);
    cand = [1 2 2.5 5 10];  i = find(cand > m + 1e-9, 1);
    if isempty(i), s2 = 10^(k+1); else, s2 = cand(i)*10^k; end
end
