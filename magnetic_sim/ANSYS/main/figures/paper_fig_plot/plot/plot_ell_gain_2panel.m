function plot_ell_gain_2panel()
% plot_ell_gain_2panel -- 降取樣（收斂點設計）結果的 l_hat 與 g_I_hat vs 取樣半徑 R
% =========================================================================
%   上下兩個子圖、共用橫軸（取樣半徑 R）：
%     上：l_hat   [um]     single（藍）與 eighteen（紅）兩條疊在一起
%     下：g_I_hat [mT/A]   同上
%   **只畫降取樣（conv）那一組**——即每個 R 各自用「收斂點設計」的少量取樣點校正的結果；
%   全格點那一組（ell_f*/gI_f*）不畫（那是 full_vs_conv_vs_R_*.png 的內容）。
%
%   資料直接讀 plot_full_vs_conv_vs_R.m 的快取 data/full_vs_conv_vs_R_maxwell.mat
%   （欄位 ell_c1/ell_c2、gI_c1/gI_c2、R_um），**不重跑任何擬合**。
%
%   **橫軸（取樣範圍 R）自 0 起**（使用者指定）：資料只有 R = 80~500，
%   故在 R=0 補一個**初始值 0** 的錨點（兩個子圖皆然），曲線自原點拉起。
%   縱軸連帶自 0 起，照 figure-style「自 0 起的軸上緣只留 8% 裕度」。
%
%   風格①粗體框圖：FS 36 粗體、box on、grid off、刻度奇數等距、曲線首末點貼齊左右框邊、
%   端點數字以 text 補、圖例照 figure-style 標準樣式（共用一個、置於最上方）。
%   配色：single 藍 / eighteen 紅（顏色 = 模型）。
%
%   輸出 → figures/paper_fig/Section2_E/ell_gain_2panel_maxwell.png（覆蓋迭代）
% =========================================================================
    clc;
    here   = fileparts(fileparts(mfilename('fullpath')));       % → paper_fig_plot/
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section2_E');
    if ~exist(figdir,'dir'); mkdir(figdir); end

    df = fullfile(here, 'data', 'full_vs_conv_vs_R_maxwell.mat');
    assert(exist(df,'file')==2, '找不到 %s —— 先跑 plot_full_vs_conv_vs_R(true)', df);
    D = load(df);
    R = D.R_um(:).';

    fprintf('R  = %s\n', num2str(R,'%6d'));
    fprintf('l_hat  single  : %.1f ~ %.1f um   (N_c %d ~ %d)\n', min(D.ell_c1), max(D.ell_c1), min(D.n_c1), max(D.n_c1));
    fprintf('l_hat  eighteen: %.1f ~ %.1f um   (N_c %d ~ %d)\n', min(D.ell_c2), max(D.ell_c2), min(D.n_c2), max(D.n_c2));
    fprintf('g_I    single  : %.3f ~ %.3f mT/A\n', min(D.gI_c1), max(D.gI_c1));
    fprintf('g_I    eighteen: %.3f ~ %.3f mT/A\n', min(D.gI_c2), max(D.gI_c2));

    % ---- 樣式 ----
    FS = 36;  LW = 3.0;  MS = 10;  LWBOX = 2.5;
    c1 = [0.05 0.10 0.95];      % single   藍
    c2 = [0.85 0.10 0.10];      % eighteen 紅

    % 手動 axes 定位（**不用 tiledlayout**：它不允許事後設 Position，
    %   圖例就無法照 figure-style 切齊座標框）。
    % [MODIFIED 2026-08-14] 改**左右排列**（使用者指定）：左 = l_hat、右 = g_I，兩格等寬等高，
    %   各自有完整的橫軸刻度與軸標題。
    fig = figure('Color','w','Position',[60 60 2000 900]);
    AXL  = 0.075;   AXW = 0.395;   GAPH = 0.115;   % 左緣 / 每格寬 / 兩格水平間距
    BOT  = 0.165;   TOP1 = 0.885;                  % 下緣（留 xlabel）/ 上緣（留圖例）
    PH   = TOP1 - BOT;
    AXL2 = AXL + AXW + GAPH;                       % 右格左緣
    ax1 = axes(fig, 'Position', [AXL,  BOT, AXW, PH]);   % 左：l_hat
    ax2 = axes(fig, 'Position', [AXL2, BOT, AXW, PH]);   % 右：g_I

    % [MODIFIED 2026-08-14] 每個 R 都是**真的跑過一次校正**（快取 full_vs_conv_vs_R_maxwell.mat，
    %   由 plot_full_vs_conv_vs_R 逐 R 搜尋收斂設計後擬合）。不補任何假錨點。
    %   橫軸自 0 起（使用者指定）；R=0 球內 0 點無法擬合，故 0 ~ 第一個 R 之間沒有資料點。
    % [MODIFIED 2026-08-14] 曲線自 **R=0 的起始點**拉起（使用者指定）：
    %   l_hat 的起始點 = 擬合初值 **500 um**；g_I 的起始點 = **0**。
    %   R>0 的每一點都是真的跑過一次校正（R = 20:20:500）。
    Rp = [0,   R];
    e1 = [500, D.ell_c1];   e2 = [500, D.ell_c2];
    g1 = [0,   D.gI_c1];    g2 = [0,   D.gI_c2];
    % 水平軸：整數刻度（使用者指定）
    XL = [0 R(end)];
    XT = 0:100:500;                          % 0 / 100 / 200 / 300 / 400 / 500

    % ===================== 上：l_hat =====================
    hold(ax1,'on');
    h1 = plot(ax1, Rp, e1, '-o', 'Color',c1, 'LineWidth',LW, 'MarkerSize',MS, ...
              'MarkerFaceColor',c1, 'Clipping','off');
    h2 = plot(ax1, Rp, e2, '-s', 'Color',c2, 'LineWidth',LW, 'MarkerSize',MS, ...
              'MarkerFaceColor','w', 'Clipping','off');
    % [MODIFIED 2026-08-14] 使用者指定：l_hat 縱軸 **自 500 起**，且刻度要把曲線夾在中間
    %   （資料 797~891 → 上下各有 700 / 900 兩個刻度）
    YL1 = [500 930];   YT1 = 500:100:900;   % 上緣不留過多空白（最上刻度 900 之上留 30）
    style_panel(ax1, XL, XT, YL1, YT1, FS, LWBOX);
    ylabel(ax1, '$\mathbf{\hat{\ell}\;(micro\;meter)}$', 'Interpreter','latex', 'FontSize',FS);
    xlabel(ax1, '$\mathbf{Sampling\;range\;(micro\;meter)}$', 'Interpreter','latex', 'FontSize',FS);
    hold(ax1,'off');

    % ===================== 下：g_I_hat =====================
    hold(ax2,'on');
    % [MODIFIED 2026-08-14] g_I 縱軸自 **0** 起（起始點也是 0），刻度用**整數**
    %   （使用者指定：不要小數點）。上緣見下方 2026-08-17 的修訂。
    % [MODIFIED 2026-08-14] **所有點都畫**（使用者指定）—— 含 single 在 R=20/40 的
    %   小球近簡併值（六顆電荷近簡併、最小平方解的 G 下極對角整層翻負）。
    % [MODIFIED 2026-08-17] solve_current 改成 ĝ_I=(6/5)|H_I(1,1)| → R=20/40 由
    %   −40.13 / −14.85 變成 **+40.13 / +14.85**，資料範圍改為 0 ~ 40.13。
    %   縱軸改自 0 起（figure-style「自 0 起的軸上緣只留 8% 裕度」）：
    %   top = 1.08×40.13 = 43.3；刻度取**整數、奇數個、等距**的 8:8:40。
    YL2 = [0 43.4];   YT2 = 8:8:40;
    plot(ax2, Rp, g1, '-o', 'Color',c1, 'LineWidth',LW, 'MarkerSize',MS, ...
         'MarkerFaceColor',c1, 'Clipping','off');
    plot(ax2, Rp, g2, '-s', 'Color',c2, 'LineWidth',LW, 'MarkerSize',MS, ...
         'MarkerFaceColor','w', 'Clipping','off');
    style_panel(ax2, XL, XT, YL2, YT2, FS, LWBOX);
    ylabel(ax2, '$\mathbf{{}^{B}\hat{g}_{I}\;(mT/A)}$', 'Interpreter','latex', 'FontSize',FS);
    xlabel(ax2, '$\mathbf{Sampling\;range\;(micro\;meter)}$', 'Interpreter','latex', 'FontSize',FS);

    % [MODIFIED 2026-08-14] 端點 0 / 500 現在本身就是 tick（水平軸 5 個 tick 含兩端），
    %   不再另用 text 補數字。
    hold(ax2,'off');

    % ===================== 共用圖例（置於最上方）=====================
    lg = legend(ax1, [h1 h2], {'Single parameter', 'Eighteen parameters'}, ...
                'Interpreter','tex', 'Location','northoutside', 'NumColumns',2);
    lg.FontSize = 24;  lg.FontWeight = 'bold';
    lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = 2.5;

    % 圖例橫跨兩格：左緣切齊左格、右緣切齊右格
    drawnow;
    lgp = get(lg,'Position');   lgw = lgp(3);   lgh = lgp(4);
    GAPN = 0.020;   SPAN = AXL2 + AXW - AXL;      % 兩格合起來的總寬
    if lgw < 0.70*SPAN                            % 內容窄 → 保持自然寬度、對整體置中
        set(lg, 'Position', [AXL + (SPAN-lgw)/2, TOP1 + GAPN, lgw, lgh]);
    else
        set(lg, 'Position', [AXL, TOP1 + GAPN, SPAN, lgh]);
    end

    out = fullfile(figdir, 'ell_gain_2panel_maxwell.png');
    exportgraphics(fig, out, 'Resolution', 200);
    fprintf('wrote %s\n', out);
    close(fig);
end

% ============================================================================
function style_panel(ax, XL, XT, YL, YT, FS, LWBOX)
    box(ax,'on');  grid(ax,'off');
    set(ax, 'FontSize',FS, 'FontWeight','bold', 'LineWidth',LWBOX, 'TickLength',[.02 .02]);
    xlim(ax, XL);   set(ax, 'XTick', XT);
    ylim(ax, YL);   set(ax, 'YTick', YT);
    ax.Toolbar.Visible = 'off';
end

% ============================================================================
function [lim, tk] = axlim_from_zero(maxv)
% 自 0 起的軸（figure-style）：下緣固定 0、上緣只留 8% 裕度；
%   刻度間距由**資料範圍**取 nice（不隨上緣壓縮而變小）；刻度數取不超出上緣的最大奇數。
    cand = [1 2 2.5 3 4 5 10];
    x = maxv/4;   k = floor(log10(x));
    s = cand(find(cand*10^k >= x, 1)) * 10^k;
    top = 1.08 * maxv;
    n = floor(top/s);   if mod(n,2)==0, n = n - 1; end
    lim = [0 top];   tk = (1:n)*s;
end

% ============================================================================
function [lim, tk] = axlim_span(lo, hi)
% 涵蓋 lo~hi 的軸：奇數個等距 nice 刻度、兩端留白 = 刻度間距（figure-style）。
    cand = [1 2 2.5 3 4 5 10];
    mid  = (lo+hi)/2;   rng_ = max(hi-lo, realmin);
    for n = [3 5]
        for k = floor(log10(rng_/(n+1))) + (0:4)
            for c = cand
                s   = c*10^k;
                ctr = round(mid/s)*s;
                t   = ctr + (-(n-1)/2 : (n-1)/2)*s;
                L   = [t(1)-s, t(end)+s];
                if lo >= L(1)+0.15*s && hi <= L(2)-0.15*s
                    lim = L;   tk = t;   return;
                end
            end
        end
    end
    s = rng_/4;   t = mid + (-1:1)*s;   lim = [t(1)-s, t(end)+s];   tk = t;
end

% ============================================================================
function tk = ticks_lin_inner(lo, hi, n)
% 線性軸：n 個（奇數）等距 nice 刻度，嚴格落在 (lo,hi) 內部、不含端點。
%   n 是目標個數；湊不出就依序退到 n+2 / n-2 / n+4（都維持奇數）。
    cand = [1 1.5 2 2.5 3 4 5 10];   % 含 1.5：0~500 才湊得出 3 個等距（150/300/450）
    for nn = [n, n+2, n-2, n+4]
        if nn < 1, continue; end
        k0 = floor(log10(max(hi-lo, realmin)/(nn+1)));
        for k = k0:(k0+3)
            for c = cand
                s = c*10^k;
                t = s*ceil((lo + 0.5*s)/s) : s : (hi - 0.5*s);
                if numel(t) == nn, tk = t;  return; end
            end
        end
    end
    tk = lo + (1:n)/(n+1)*(hi-lo);
end
