function plot_sensor_normal_gradient(NL)
% plot_sensor_normal_gradient -- Hall sensor 感測體積內「沿 n+ 的場梯度曲線」(Maxwell)
% =========================================================================
%   P1 sensor 貼**削平面**(FACE_lower='flat'、n̂⁺=+z)、P2 sensor 貼**錐面**，
%   兩者皆 SOFF = 4.572mm、離面 0.41mm、圓柱 R 0.15 x H 0.10mm、撒 500 點。
%   把圓柱沿 n̂⁺ 切 NL 層(預設 10)，**每層有幾個點就用幾個點取平均** <‖b‖>，
%   對「名目層中心」作圖 -> 看得出兩個貼附面的沿法線衰減差異。
%
%   ★ 資料一律取自 plot_sensor_B_hist 的 per-pole 快取
%     (`sensor_B_hist_P<k>_maxwell_soff4.572_n500.mat`)，與
%     `plot_sensor_B_hist`(直方圖) 及 `plot_sensor_cyl_B_3d`(3D 箭頭) **是同一批點**
%     -> 三張圖可逐點互相對照。本檔不自建快取、不讀 .fld。
%     快取不存在時先跑 `plot_sensor_B_hist(true)`。
%
%   縱軸 = <‖b‖> [mT]（**與直方圖同一個量**，非 b·n̂⁺ 投影；使用者拍板 2026-08-06 統一）。
%
%   ⚠ 這是**內插**(Maxwell 匯出為 0.1mm 規則格，sensor 圓柱只有 0.15mm，
%     不內插取不到點)，與校正管線 V_METHOD='scattered' 同一套。
%
%   風格:選項①粗體框圖。橫軸依 figure-style 2026-08-06 新標準:
%     首末資料點貼齊左右框邊(無空白)、端點只標數字不畫 tick、內部 tick 奇數等距、
%     分層資料 x 用名目層中心。
% =========================================================================
    clc;
    if nargin < 1, NL = 10; end                      % 層數

    here   = fileparts(fileparts(mfilename('fullpath')));          % .../paper_fig_plot
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section3_A');
    if ~exist(figdir,'dir'); mkdir(figdir); end

    S1 = load_case(here, 1, 'p1flat');               % P1 削平面（P1 自激發）
    S2 = load_case(here, 2, 'p2cone');               % P2 錐面  （P2 自激發）

    % ---- 分層平均:每層有幾點就用幾點；x = 名目層中心 ----
    %   AT 用**名目圓柱高**(build_V_matrix 的 axial_tol = 0.10mm)，不可用撒點實際最大值
    %   ——後者是隨機的(99.79um)，會讓層中心變成 4.98965/14.969… 這種非整數。
    AT  = 100;                                       % [um] = axial_tol 0.10mm
    edg = linspace(0, AT, NL+1);                     % 0,10,…,100 -> 中心 5,15,…,95
    [c1, m1, n1] = layer_mean(S1.a, S1.b, edg);
    [c2, m2, n2] = layer_mean(S2.a, S2.b, edg);

    fprintf('layer   a(um)   P1flat<‖b‖>  nP1    P2cone<‖b‖>  nP2\n');
    for L = 1:NL
        fprintf('%5d %8.1f %12.4f %6d %14.4f %6d\n', L, c1(L), m1(L), n1(L), m2(L), n2(L));
    end
    fprintf('P1 flat : %.4f -> %.4f mT  (drop %.2f%%)\n', m1(1), m1(end), (m1(1)-m1(end))/m1(1)*100);
    fprintf('P2 cone : %.4f -> %.4f mT  (drop %.2f%%)\n', m2(1), m2(end), (m2(1)-m2(end))/m2(1)*100);
    fprintf('volume mean: P1 %.4f  P2 %.4f mT  (P2/P1 = %.4f, %+.2f%%)\n', ...
            mean(S1.b), mean(S2.b), mean(S2.b)/mean(S1.b), (mean(S2.b)/mean(S1.b)-1)*100);

    % ---- 繪圖 ----
    FS  = 36;
    BLU = [0.10 0.35 1.00];      % P1 削平面（與直方圖同色系）
    RED = [0.85 0.10 0.10];      % P2 錐面
    fig = figure('Color','w','Position',[100 100 1150 900]);  ax = axes(fig);  hold(ax,'on');

    % Clipping off:首末點正好落在框邊上,不要被框線切一半
    h1 = plot(ax, c1, m1, '-o', 'Color',BLU, 'LineWidth',3.5, ...
              'MarkerSize',13, 'MarkerFaceColor',BLU, 'MarkerEdgeColor','k', 'Clipping','off');
    h2 = plot(ax, c2, m2, '-s', 'Color',RED, 'LineWidth',3.5, ...
              'MarkerSize',14, 'MarkerFaceColor',RED, 'MarkerEdgeColor','k', 'Clipping','off');

    box(ax,'on');  grid(ax,'off');
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',2.5,'TickLength',[.015 .015],'TickDir','out');

    % 橫軸:首末點貼齊框邊、端點只標數字不畫 tick、內部 tick 奇數等距
    x1 = c1(1);  x2 = c1(end);
    xlim(ax, [x1, x2]);
    sx = round(AT/4/5)*5;                            % 內部 tick 步長（~25um）
    set(ax,'XTick', sx*(1:3));

    % 縱軸:奇數個、等距、兩端留白 = 間距(首末不標)
    [yr, yt] = ylim_pick([m1(:); m2(:)]);
    ylim(ax, yr);  set(ax,'YTick',yt);

    % 端點數字(無 tick mark)
    yoff = yr(1) - 0.022*diff(yr);
    text(ax, x1, yoff, sprintf('%g',x1), 'HorizontalAlignment','center', ...
         'VerticalAlignment','top','FontSize',FS,'FontWeight','bold','Clipping','off');
    text(ax, x2, yoff, sprintf('%g',x2), 'HorizontalAlignment','center', ...
         'VerticalAlignment','top','FontSize',FS,'FontWeight','bold','Clipping','off');

    xlabel(ax, '$\mathbf{Position\;along\;\hat{n}^{+}\;(micro\;meter)}$', ...
           'Interpreter','latex', 'FontSize',FS);
    ylabel(ax, '$\mathbf{\langle\|b\|\rangle\;(mT)}$', ...
           'Interpreter','latex', 'FontSize',FS);

    lg = legend([h1 h2], {'P1 flat cut surface', 'P2 cone surface'}, ...
                'Interpreter','tex', 'Location','northeast');
    lg.FontSize = 26;  lg.FontWeight = 'bold';  lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = 2.5;

    ax.Toolbar.Visible = 'off';  hold(ax,'off');

    out = fullfile(figdir, 'sensor_normal_gradient_p1flat_p2cone_maxwell.png');
    exportgraphics(fig, out, 'Resolution', 150);
    fprintf('wrote %s\n', out);
end

% ============================================================================
function S = load_case(here, kpole, fld)
% 讀 plot_sensor_B_hist 的 per-pole 快取 -> ‖b‖[mT] 與沿 n̂⁺ 位置[um]（同一批點）
    f = fullfile(here,'data', sprintf('sensor_B_hist_P%d_maxwell_soff4.572_n500.mat', kpole));
    assert(exist(f,'file')==2, 'cache 不在:%s\n(先跑 plot_sensor_B_hist(true))', f);
    C  = load(f);
    P  = C.(['P_' fld]);  B = C.(['B_' fld]);  nh = C.(['n_' fld])(:);
    a  = P*nh;                                   % 沿 n̂⁺ 投影（快取座標為 mm）
    S.a = (a - min(a))*1e3;                      % 以最內側點為 0 [um]
    S.b = vecnorm(B, 2, 2);                      % ‖b‖ [mT]
end

% ============================================================================
function [ctr, mu, n] = layer_mean(a, b, edg)
% 每層:落在 [edg(L), edg(L+1)) 的點取平均(最後一層含右端點)。層內點數不強制相等。
%   x 用**名目層中心**(edg 中點)——不用該層撒點的實際平均位置(帶抖動、非等距)。
    NL  = numel(edg) - 1;
    ctr = zeros(NL,1);  mu = zeros(NL,1);  n = zeros(NL,1);
    for L = 1:NL
        k = a >= edg(L) & a < edg(L+1);
        if L == NL, k = k | (a >= edg(end)); end
        n(L)   = nnz(k);
        ctr(L) = (edg(L) + edg(L+1)) / 2;
        mu(L)  = mean(b(k));
    end
end

% ============================================================================
function [yr, yt] = ylim_pick(v)
% 縱軸:3 個等距 nice tick + 兩端留白 = 間距(limit 由 tick 反推,per figure-style)
    lo = min(v);  hi = max(v);
    s  = nice_step((hi-lo)/2);
    while 3*s < (hi - lo), s = next_step(s); end
    ctr = round((lo+hi)/2 / s) * s;
    yt  = ctr + (-1:1)*s;
    yr  = [yt(1)-s, yt(end)+s];
    if yr(1) > lo || yr(2) < hi
        s = next_step(s);  ctr = round((lo+hi)/2 / s) * s;
        yt = ctr + (-1:1)*s;  yr = [yt(1)-s, yt(end)+s];
    end
end

function s = nice_step(x)
    k = floor(log10(max(x, realmin)));  m = x/10^k;
    cand = [1 2 2.5 5 10];
    [~,i] = min(abs(cand - m));
    s = cand(i)*10^k;
end

function s2 = next_step(s)
    k = floor(log10(s));  m = round(s/10^k, 4);
    cand = [1 2 2.5 5 10];
    i = find(cand > m + 1e-9, 1);
    if isempty(i), s2 = 10^(k+1); else, s2 = cand(i)*10^k; end
end
