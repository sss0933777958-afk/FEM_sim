function plot_ell_gain_overlay(SRC, SAMPLING)
% plot_ell_gain_overlay -- paper 圖：ℓ̂(R) 與 ĝ_I(R) 的「single vs eighteen 疊圖」
% =========================================================================
%   把 plot_ell_gain_vs_R 產的兩張（single / eighteen）疊成一張：
%     上 panel = ℓ̂(R) [µm]、下 panel = ᴮĝ_I(R) [mT/A]，每個 panel 兩條曲線。
%   配色沿用專案 overlay 慣例（同 plot_err_hist_overlay）：**藍 = single、紅 = eighteen**
%   （注意：不同於單張圖的「藍=ℓ̂、紅=ĝ_I」——疊圖的顏色改為編碼「模型」而非「物理量」）。
%
%   **只讀 plot_ell_gain_vs_R 的 sweep 快取、不重算**（兩顆快取都要存在）。
%   軸規則（2026-08-03 定案）：兩 panel tick 數一致；縱軸兩端留白 = tick 間距、資料不貼框；
%   水平軸 = 資料範圍本身（80–500），5 個等距 tick。
%   輸出 → figures/paper_fig/Section2_E/ell_gain_vs_R_overlay_<solver>.png（覆蓋迭代）。
% =========================================================================
    clc;
    if nargin < 1, SRC = 'maxwell'; end
    if nargin < 2 || isempty(SAMPLING)
        if strcmpi(SRC,'maxwell'), SAMPLING = 'interp'; else, SAMPLING = 'raw'; end
    end
    here   = fileparts(mfilename('fullpath'));
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section2_E');
    if ~exist(figdir,'dir'); mkdir(figdir); end

    sstr = ''; if strcmpi(SRC,'maxwell'),     sstr = '_maxwell'; end
    istr = ''; if strcmpi(SAMPLING,'interp'), istr = '_interp';  end

    % ---- 讀兩顆 sweep 快取 ----
    f0 = fullfile(here, sprintf('ell_gain_sweep%s%s.mat',      sstr, istr));          % single
    f1 = fullfile(here, sprintf('ell_gain_sweep%s%s_bias.mat', sstr, istr));          % eighteen
    assert(exist(f0,'file')==2, 'single 快取不存在：%s\n先跑 plot_ell_gain_vs_R(false,''%s'')。', f0, SRC);
    assert(exist(f1,'file')==2, 'eighteen 快取不存在：%s\n先跑 plot_ell_gain_vs_R(true,''%s'')。', f1, SRC);
    S = load(f0);   E = load(f1);
    fprintf('loaded\n  single  : %s\n  eighteen: %s\n', f0, f1);

    % ---- 繪圖範圍（與單張圖同：maxwell 從 80µm 起）----
    if strcmpi(SRC,'maxwell'), Rmin = 80; else, Rmin = 40; end
    m0 = S.Rum >= Rmin;   R0 = S.Rum(m0);  ell0 = S.ell_R(m0);  g0 = S.gI_R(m0);
    m1 = E.Rum >= Rmin;   R1 = E.Rum(m1);  ell1 = E.ell_R(m1);  g1 = E.gI_R(m1);

    [~, XT] = axlim_auto(R0(1), R0(end), 5);
    XL = [R0(1), R0(end)];                      % 水平軸 = 資料範圍本身（80–500）

    ell_all = [ell0(:); ell1(:)];   g_all = [g0(:); g1(:)];        % 兩模型聯集決定 y 範圍
    NY = pick_common_n(ell_all, g_all, [3 5]);                     % 兩 panel tick 數一致

    fprintf('  ℓ̂  : single %.1f→%.1f  eighteen %.1f→%.1f um\n', ell0(1), ell0(end), ell1(1), ell1(end));
    fprintf('  ĝ_I : single %.4f→%.4f  eighteen %.4f→%.4f mT/A\n', g0(1), g0(end), g1(1), g1(end));

    % ---- 畫圖 ----
    % [MODIFIED 2026-08-03] 換一組配色（使用者要求「用不同顏色配」）：沿用專案既有的第二組
    %   紫 #7B52AB / 橘 #E69F00（= plot_err_hist_overlay 的 pick_bar_colors 非 R150 那組），
    %   不另發明新色；與藍/紅的直方圖疊圖區隔開。
    cS = [0.482 0.322 0.671];   cE = [0.902 0.624 0.000];   % 紫 = single / 橘 = eighteen
    FS = 30;
    fig = figure('Color','w','Position',[100 80 980 1000]);        % 加高補償外置圖例
    t = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');

    % 上：ℓ̂(R)
    ax1 = nexttile(t);  hold(ax1,'on');
    h0 = plot(ax1, R0, ell0, '-o', 'Color',cS, 'LineWidth',3.0, 'MarkerSize',8, 'MarkerFaceColor',cS);
    h1 = plot(ax1, R1, ell1, '-s', 'Color',cE, 'LineWidth',3.0, 'MarkerSize',8, 'MarkerFaceColor',cE);
    style_panel(ax1, FS);
    xlim(ax1, XL);   set(ax1,'XTick',XT);
    [yl,yt] = axlim_auto(min(ell_all), max(ell_all), NY);   ylim(ax1,yl);  set(ax1,'YTick',yt);
    bold_ticks(ax1);                                    % [ADDED] latex tick 要 \mathbf 才會粗
    ax1.XTickLabel = {};                                % 上 panel x 數字隱藏(共用軸)
    ylabel(ax1, '$\mathbf{\hat{\ell}\;(\mu m)}$', 'Interpreter','latex', 'FontSize',36);
    lg = legend([h0 h1], {'Single parameter','Eighteen parameters'}, ...
                'Interpreter','tex', 'Location','northoutside', 'Orientation','horizontal');
    lg.FontSize = 24;  lg.FontWeight = 'bold';  lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = 2.5;
    hold(ax1,'off');

    % 下：ĝ_I(R)
    ax2 = nexttile(t);  hold(ax2,'on');
    plot(ax2, R0, g0, '-o', 'Color',cS, 'LineWidth',3.0, 'MarkerSize',8, 'MarkerFaceColor',cS);
    plot(ax2, R1, g1, '-s', 'Color',cE, 'LineWidth',3.0, 'MarkerSize',8, 'MarkerFaceColor',cE);
    style_panel(ax2, FS);
    xlim(ax2, XL);   set(ax2,'XTick',XT);
    [yl,yt] = axlim_auto(min(g_all), max(g_all), NY);       ylim(ax2,yl);  set(ax2,'YTick',yt);
    bold_ticks(ax2);                                    % [ADDED] 同上
    ylabel(ax2, '$\mathbf{{}^{B}\hat{g}_{I}\;(mT/A)}$', 'Interpreter','latex', 'FontSize',36);
    xlabel(ax2, '$\mathbf{Sampling\;range\;(\mu m)}$', 'Interpreter','latex', 'FontSize',36);
    hold(ax2,'off');

    out = fullfile(figdir, sprintf('ell_gain_vs_R_overlay%s.png', sstr));
    exportgraphics(fig, out, 'Resolution', 200);
    fprintf('\nwrote %s\n', out);
end

% ============================================================================
function n = pick_common_n(v1, v2, nlist)
% 兩 panel 共用的 tick 數：取「兩邊都可行」且「較差那個 panel 填充率最高」者。
    best = nlist(1);  bestFill = -inf;
    for n_ = nlist
        [L1, ~] = axlim_auto(min(v1), max(v1), n_);
        [L2, ~] = axlim_auto(min(v2), max(v2), n_);
        f = min( (max(v1)-min(v1))/diff(L1), (max(v2)-min(v2))/diff(L2) );
        if f > bestFill, bestFill = f;  best = n_; end
    end
    n = best;
end

% ============================================================================
function [lim, tk] = axlim_auto(lo, hi, nlist)
% 軸範圍/刻度：奇數個等距 tick、兩端留白 = tick 間距、資料與框線留 0.15 格淨空。
%   （與 plot_npts_cost_vs_R 同一份邏輯；改動要兩邊同步。）
    cand = [1 2 2.5 3 4 5 10];
    mid  = (lo+hi)/2;   rng = max(hi-lo, realmin);
    best = {};   bestSpan = inf;
    for n = nlist
        k0 = floor(log10(rng/(n+1)));
        for k = k0:(k0+4)
            hit = false;
            for c = cand
                s   = c*10^k;
                ctr = round(mid/s)*s;
                T = { ctr + (-(n-1)/2 : (n-1)/2)*s };
                if lo >= 0 && lo < s, T = [{(0:n-1)*s}, T]; end %#ok<AGROW>
                for it = 1:numel(T)
                    t   = T{it};
                    L   = [t(1)-s, t(end)+s];
                    clr = 0.15*s;
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
        n = nlist(1);  s = rng/(n+1);
        t = mid + (-(n-1)/2 : (n-1)/2)*s;   best = {[t(1)-s, t(end)+s], t};
    end
    lim = best{1};   tk = best{2};
end

% ============================================================================
function bold_ticks(ax)
% [ADDED 2026-08-03] 讓 LaTeX 刻度數字變粗體。
%   ⚠ `TickLabelInterpreter='latex'` 時 MATLAB **不吃 `FontWeight='bold'`**（永遠 regular），
%   唯一辦法是把每個數字自己包成 $\mathbf{...}$。呼叫時機必須在 set(ax,'XTick'/'YTick') 之後
%   （之後若再改 tick，label 會被 MATLAB 重設回非粗體，得重呼叫一次）。
    xt = get(ax,'XTick');   set(ax,'XTickLabel', arrayfun(@(v) sprintf('$\\mathbf{%g}$', v), xt, 'UniformOutput', false));
    yt = get(ax,'YTick');   set(ax,'YTickLabel', arrayfun(@(v) sprintf('$\\mathbf{%g}$', v), yt, 'UniformOutput', false));
end

% ============================================================================
function style_panel(ax, FS)
% 選項①粗體框 @ paper scale；刻度數字用 LaTeX Computer Modern（2026-08-03 統一）。
    box(ax,'on');  grid(ax,'off');
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',2.5,'TickLength',[.02 .02]);
    set(ax,'TickLabelInterpreter','latex');
    ax.Toolbar.Visible = 'off';
end
