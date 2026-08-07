function plot_upper_boundary_flux(force, WHICH, SRC)
%   [ADDED 2026-08-05] WHICH='p1p2full' + SRC='maxwell'：
%     每極的**完整側表面**累積漏磁 —— P1（半切下極）= 削平面 + 下錐面；P2（上極）= 完整錐面（上+下半）。
%     SRC='apdl'（預設）| 'maxwell'（0.1mm 規則格；⚠ 見 plot_surface_flux 檔頭的表面積分警告）。
%     Maxwell 匯出已是 all-source（cfg.s_source 全 +1）→ 兩極都不翻號；APDL 則下極 −1、上極 +1。
%PLOT_UPPER_BOUNDARY_FLUX  下極「上邊界」累積漏磁 Φ_cum(x) 疊圖：填滿(上半錐面) vs 半切(削平面)。
%   plot_upper_boundary_flux()      % 有快取就用
%   plot_upper_boundary_flux(true)  % 強制重算
%
%   兩條曲線是**同一個下極上邊界**在兩種設計下的樣子：
%     filled  → 上半錐面（lower_filled/coil1，graded）
%     halfcut → 削平面 z = -13.000（graded/coil1，baseline）
%   Φ_cum(x) = 從尖端 (x=0.408) 沿軸積到 x 的漏磁量 [µWb]；終點 x = 15.235（錐底）。
%   兩者相減 = 那一刀對上邊界漏磁的淨效果。
%
%   號誌：long2016 下極 raw FEM 是 **sink**（磁通朝尖端流入）→ 依
%   .claude/rules/charge-model-source-convention.md 對 P1 套 s_source = -1 轉成 all-source，
%   此後 Φ > 0 = 磁通從該面**漏出去**。
%
%   ⚠ 內插：格心 B 由 FEM 鋼件節點內插（見 plot_surface_flux.m 檔頭聲明），非節點原值。
%
%   風格①粗體框圖。輸出 figures/paper_fig/Section3_A/upper_boundary_flux_cum.png
%   快取 figures/paper_fig_plot/upper_boundary_flux.mat

    if nargin < 1 || isempty(force), force = false; end
    % [ADDED 2026-08-05] WHICH = 比哪個面：
    %   'upper'（預設）= 上邊界：filled 上半錐面 vs halfcut 削平面（**兩個不同的面**，設計層級比較）
    %   'lower'        = 下錐面：兩個幾何**同一個面、同一組座標、同一法線**（半切只削上半）
    %                    → 唯一能隔離「磁極本身」的受控對照，無位置/角度混淆
    %   'p1p2'         = P1 削平面（coil1 自激）vs P2 上半錐面（coil5 自激）——同一顆
    %                    baseline 半切六極、各自自激發，隔離「上下極 + 面型」的差異
    if nargin < 2 || isempty(WHICH), WHICH = 'upper'; end
    if nargin < 3 || isempty(SRC),   SRC   = 'apdl';  end        % [ADDED]
    LOW = strcmpi(WHICH,'lower');  P12 = strcmpi(WHICH,'p1p2');  FULL = strcmpi(WHICH,'p1p2full');
    MX  = strcmpi(SRC,'maxwell');
    HERE   = fileparts(fileparts(mfilename('fullpath')));
    FIGDIR = fullfile(fileparts(HERE), 'paper_fig', 'Section3_A');
    ssfx = ''; if MX, ssfx = '_maxwell'; end                     % [ADDED] 來源進快取/檔名
    CACHE  = fullfile(HERE, 'data', sprintf('%s_boundary_flux%s.mat', lower(WHICH), ssfx));
    if ~exist(FIGDIR,'dir'), mkdir(FIGDIR); end

    XT = 0.408;  XE = 15.235;  SGN = -1;        % SGN：下極 sink → all-source
    NB = 300;                                    % 沿軸累積的分箱數

    if ~force && exist(CACHE,'file')
        S = load(CACHE);  fprintf('cache hit: %s\n', CACHE);
    else
        % sgF/sgH = 各曲線的 s_source（下極 sink → -1；上極 source → +1）
        if FULL
            % [ADDED] 整根極的完整側表面：P1 = 削平面 + 下錐面；P2 = 上半錐 + 下半錐
            %   Maxwell: coil→pole identity（p1/p2）、已 all-source（不翻號）
            %   APDL   : P1=coil1、P2=coil5（map [1,3,6,5,2,4]）、下極 raw 為 sink（−1）
            if MX
                c1 = ''; c2 = '';  sgF = +1;  sgH = +1;
            else
                c1 = 'coil1'; c2 = 'coil5';  sgF = +1;  sgH = -1;
            end
            oF2 = plot_surface_flux('cone_upper','graded',[],[],[],2,c2,SRC);   % P2 上半錐
            oF3 = plot_surface_flux('cone_lower','graded',[],[],[],2,c2,SRC);   % P2 下半錐
            oH2 = plot_surface_flux('flat',      'graded',[],[],[],1,c1,SRC);   % P1 削平面
            oH3 = plot_surface_flux('cone_lower','graded',[],[],[],1,c1,SRC);   % P1 下錐面
            oF = merge_surf(oF2, oF3);   oH = merge_surf(oH2, oH3);             % 同極兩面合併
        elseif P12
            oF = plot_surface_flux('cone_upper','graded',[],[],[],2,'coil5');  % P2 上半錐面
            oH = plot_surface_flux('flat',      'graded',[],[],[],1,'coil1');  % P1 削平面
            sgF = +1;  sgH = -1;
        elseif LOW
            oF = plot_surface_flux('cone_lower','lower_filled'); % 填滿：下錐面
            oH = plot_surface_flux('cone_lower','graded');       % 半切：下錐面（同一個面）
            sgF = SGN;  sgH = SGN;
        else
            oF = plot_surface_flux('cone_upper','lower_filled'); % 填滿：上半錐面
            oH = plot_surface_flux('flat',      'graded');       % 半切：削平面
            sgF = SGN;  sgH = SGN;
        end
        % [MODIFIED] 改用**局部軸向 s**（0 = 極尖、S_CONE = 錐底），兩極各自的軸
        ed = linspace(0, oH.S_CONE, NB+1);  xc = (ed(1:end-1)+ed(2:end))/2;
        S = struct('xc',xc, ...
                   'cumF',cumbin(oF.s_cell, oF.dPhi, ed)*sgF, ...
                   'cumH',cumbin(oH.s_cell, oH.dPhi, ed)*sgH, ...
                   'PhiF',oF.Phi*sgF, 'PhiH',oH.Phi*sgH, 'oF',oF, 'oH',oH);
        save(CACHE, '-struct', 'S');
        fprintf('saved %s\n', CACHE);
    end

    %% ---- 畫圖 ----
    col = [0.85 0.33 0.10;    % filled 上半錐面 橘
           0.00 0.45 0.74];   % halfcut 削平面  藍
    fig = figure('Position',[80 80 1250 900],'Color','w');
    ax  = axes(fig);  hold(ax,'on');
    % [MODIFIED] 兩端補上積分域邊界（x=XT 累積為 0、x=XE 為總量），曲線精確跨滿 [0.408, 15.235]
    xp = [0; S.xc(:); S.oH.S_CONE];
    h1 = plot(ax, xp, [0; S.cumF(:); S.PhiF], '-', 'Color', col(1,:), 'LineWidth', 3.5);
    h2 = plot(ax, xp, [0; S.cumH(:); S.PhiH], '-', 'Color', col(2,:), 'LineWidth', 3.5);
    hold(ax,'off');

    % [MODIFIED] 水平軸對齊曲線起訖（尖端 0.408 → 錐底 15.235），不留空白。
    %   端點落在框角、該處不畫 tick line → 依 figure-style 慣例#5 用 text 補數字，
    %   XTick 只放內部等距值（3 個、奇數、間距 4）。
    xl = [0 S.oH.S_CONE];  xt = [4 8 12];
    [yt, yl] = ticks_y(max([S.cumF(:); S.cumH(:)]));
    xlim(ax, xl);  ylim(ax, yl);
    set(ax, 'XTick', xt, 'YTick', yt);
    set(ax, 'FontSize',36, 'FontWeight','bold', 'LineWidth',2, 'TickLength',[.018 .018]);
    box(ax,'on');  grid(ax,'off');
    for e = 1:2
        text(ax, xl(e), -0.022*yl(2), sprintf('%.2f', xl(e)), 'FontSize',36, 'FontWeight','bold', ...
             'HorizontalAlignment','center', 'VerticalAlignment','top', 'Clipping','off');
    end
    xlabel(ax, '$\mathbf{Axial\;position\;(mm)}$',        'Interpreter','latex', 'FontSize',36);
    ylabel(ax, '$\mathbf{\Phi\;(\mu Wb)}$',               'Interpreter','latex', 'FontSize',36);
    if FULL,    LG = {'P2: full cone surface','P1: flat cut + lower cone'};   % [ADDED]
    elseif P12, LG = {'P2: upper cone face','P1: flat cut face'};
    elseif LOW, LG = {'Filled: lower cone face','Half-cut: lower cone face'};
    else,       LG = {'Filled: upper cone face','Half-cut: flat face'};  end
    lg = legend([h1 h2], LG, 'Location','northwest');
    lg.FontSize = 28;  lg.FontWeight = 'bold';  lg.Box = 'off';
    ax.Toolbar.Visible = 'off';

    out = fullfile(FIGDIR, sprintf('%s_boundary_flux_cum%s.png', lower(WHICH), ssfx));   % [MODIFIED] 檔名帶來源
    exportgraphics(fig, out, 'Resolution', 150);
    fprintf('wrote %s\n', out);  close(fig);

    fprintf('\n--- 表面總漏磁（all-source，正 = 漏出；src=%s）---\n', lower(SRC));
    fprintf('  filled  %s : %+.5f uWb  (面積 %.3f mm²)\n', LG{1}, S.PhiF, S.oF.A_ana);
    fprintf('  halfcut %s : %+.5f uWb  (面積 %.3f mm²)\n', LG{2}, S.PhiH, S.oH.A_ana);
    fprintf('  差    filled - halfcut = %+.5f uWb  (%+.1f%%)\n', ...
            S.PhiF-S.PhiH, (S.PhiF/S.PhiH-1)*100);
end

% ---- [ADDED] 同一顆極的兩個面合併成一條曲線（格併格；面積/總量相加）----
function o = merge_surf(a, b)
    o = a;
    o.s_cell = [a.s_cell(:); b.s_cell(:)];
    o.dPhi   = [a.dPhi(:);   b.dPhi(:)];
    o.Phi    = a.Phi + b.Phi;
    o.A_ana  = a.A_ana + b.A_ana;
    o.ncell  = a.ncell + b.ncell;
end

% ---- 沿軸累積（把每格的 dPhi 依 x 分箱後累加）----
function c = cumbin(xv, dv, ed)
    i = discretize(xv, ed);  ok = ~isnan(i);
    c = cumsum(accumarray(i(ok), dv(ok), [numel(ed)-1, 1]));
end

% ---- 縱軸：奇數個等距 tick、首末不標、兩端留白 = 間距 ----
function [tk, lim] = ticks_y(ymax)
    nice = [1 2 2.5 3 4 5 10];
    for e = -4:2
        for q = nice
            s = q*10^e;
            if 4*s >= ymax*1.001, tk = [s 2*s 3*s];  lim = [0 4*s];  return, end
        end
    end
    s = ymax/3;  tk = [s 2*s 3*s];  lim = [0 4*s];
end
