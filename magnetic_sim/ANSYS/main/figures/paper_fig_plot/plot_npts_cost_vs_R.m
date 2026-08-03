function plot_npts_cost_vs_R(USE_BIAS, SRC, SAMPLING)
% plot_npts_cost_vs_R -- paper 圖：取樣點數 N 與擬合 cost J 隨取樣半徑 R 的趨勢
% =========================================================================
%   上 panel = 取樣點數 N(R)（無單位 → 不標單位）
%   下 panel = cost J(R) [mT²]
%       J = Σ_{j=1..6} ‖ S(ℓ̂,Pc)·g_j − b_j ‖²，g_j = (SᵀS)⁻¹Sᵀb_j（變數投影消去電荷強度）
%       即 fitting.m 最小化並回傳的目標函數；求和涵蓋取樣球內所有點×3 分量×6 激發。
%   ⚠ J 是外延量：會隨 N 一起長大（上 panel 正是給讀者看這件事）。
%
%   **本圖只讀 plot_ell_gain_vs_R 產生的 sweep 快取、不重算**（快取需含 J_R；
%   若舊快取沒有 J_R，先跑一次 plot_ell_gain_vs_R 重建）。
%   風格 = 選項①粗體框 @ paper scale，與 ell_gain_vs_R 完全一致（同字體/框/刻度規則）。
%   輸出 → figures/paper_fig/Section2_E/npts_cost_vs_R[_bias][_maxwell].png（覆蓋迭代）。
% =========================================================================
    clc;
    if nargin < 1, USE_BIAS = false; end   % false = fix(無 e)；true = 18-param bias
    if nargin < 2, SRC = 'maxwell';  end   % 'apdl' | 'maxwell'
    if nargin < 3 || isempty(SAMPLING)     % 與 ell_gain_vs_R 同慣例：maxwell 預設 interp
        if strcmpi(SRC,'maxwell'), SAMPLING = 'interp'; else, SAMPLING = 'raw'; end
    end
    here   = fileparts(mfilename('fullpath'));
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section2_E');
    if ~exist(figdir,'dir'); mkdir(figdir); end

    bstr = ''; if USE_BIAS, bstr = '_bias'; end
    sstr = ''; if strcmpi(SRC,'maxwell'),      sstr = '_maxwell'; end
    istr = ''; if strcmpi(SAMPLING,'interp'),  istr = '_interp';  end

    % ---- 讀 sweep 快取（與 plot_ell_gain_vs_R 同一顆）----
    cachef = fullfile(here, sprintf('ell_gain_sweep%s%s%s.mat', sstr, istr, bstr));
    assert(exist(cachef,'file')==2, ...
        'sweep 快取不存在：%s\n先跑 plot_ell_gain_vs_R(%d,''%s'')。', cachef, USE_BIAS, SRC);
    S = load(cachef);
    assert(isfield(S,'J_R'), ...
        '快取缺 J_R（舊版）：%s\n刪掉它並重跑 plot_ell_gain_vs_R(%d,''%s'') 重建。', cachef, USE_BIAS, SRC);
    Rum = S.Rum;  npts_R = S.npts_R;  J_R = S.J_R;
    fprintf('loaded %s\n', cachef);

    % ---- 繪圖範圍（與 ell_gain_vs_R 同：maxwell 從 80µm 起，小 R 格點太少擬合病態）----
    if strcmpi(SRC,'maxwell'), Rmin = 80; else, Rmin = 40; end
    m = Rum >= Rmin;   Rum = Rum(m);  npts_R = npts_R(m);  J_R = J_R(m);
    XT = [100 200 300 400 500];   % 等距(見 figure-style tick 均勻規則);x 軸起點對齊第一 tick(100),避免大留白/壓字

    fprintf('  R=%d: N=%d, J=%.4g mT^2   →   R=%d: N=%d, J=%.4g mT^2\n', ...
            Rum(1), npts_R(1), J_R(1), Rum(end), npts_R(end), J_R(end));

    % ---- 畫圖(2×1 panel,選項①粗體框；與 ell_gain_vs_R 逐項一致)----
    cN = [0.05 0.10 0.95];   cJ = [0.85 0.10 0.10];     % N 亮藍 / J 紅
    FS = 30;
    fig = figure('Color','w','Position',[100 80 980 940]);
    t = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');

    % 上：N(R)
    ax1 = nexttile(t);  hold(ax1,'on');
    plot(ax1, Rum, npts_R, '-o', 'Color',cN, 'LineWidth',3.0, 'MarkerSize',8, 'MarkerFaceColor',cN);
    style_panel(ax1, FS);
    xlim(ax1,[XT(1) 500]);   set(ax1,'XTick',XT);
    [yl,yt] = ylim_auto(npts_R);   ylim(ax1,yl);  set(ax1,'YTick',yt);
    ax1.XTickLabel = {};                                % 上 panel x 數字隱藏(共用軸)
    ylabel(ax1, '$\mathbf{Number\;of\;points}$', 'Interpreter','latex', 'FontSize',36);   % 無單位 → 不標
    hold(ax1,'off');

    % 下：J(R)
    ax2 = nexttile(t);  hold(ax2,'on');
    plot(ax2, Rum, J_R, '-s', 'Color',cJ, 'LineWidth',3.0, 'MarkerSize',8, 'MarkerFaceColor',cJ);
    style_panel(ax2, FS);
    xlim(ax2,[XT(1) 500]);   set(ax2,'XTick',XT);
    [yl,yt] = ylim_auto(J_R);      ylim(ax2,yl);  set(ax2,'YTick',yt);
    ylabel(ax2, '$\mathbf{Cost\;(mT^{2})}$', 'Interpreter','latex', 'FontSize',36);
    xlabel(ax2, '$\mathbf{Sampling\;range\;(\mu m)}$', 'Interpreter','latex', 'FontSize',36);
    hold(ax2,'off');

    out = fullfile(figdir, sprintf('npts_cost_vs_R%s%s.png', bstr, sstr));
    exportgraphics(fig, out, 'Resolution', 200);
    fprintf('\nwrote %s\n', out);
end

% ============================================================================
function [yl, yt] = ylim_auto(v)
% y 範圍 / 刻度：4 個等距內縮 tick（不含端點），照 figure-style「等距、不含端點」。
    lo = min(v(:));  hi = max(v(:));  sp = hi - lo;
    yl = [lo - 0.12*sp, hi + 0.12*sp];
    s  = nice_step(diff(yl)/4);                          % 3 個 tick 的間距(margin≈間距、不擠;見 figure-style tick 規則)
    c  = round(mean(yl)/s)*s;
    yt = c + (-1:1:1)*s;                                 % 3 個等距 tick(短 panel 不塞太多)
    yt = yt(yt > yl(1) & yt < yl(2));
end

function s = nice_step(x)
    k = floor(log10(x));   m = x/10^k;
    cand = [1 2 2.5 3 4 5 10];
    [~,i] = min(abs(cand - m));
    s = cand(i)*10^k;
end

% ============================================================================
function style_panel(ax, FS)
% 選項①粗體框 @ paper scale(tick 由呼叫端明設,不減半)。
    box(ax,'on');  grid(ax,'off');
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',2.5,'TickLength',[.02 .02]);   % 刻度數字 Helvetica 粗體(維持原字體)
    ax.Toolbar.Visible = 'off';
end
