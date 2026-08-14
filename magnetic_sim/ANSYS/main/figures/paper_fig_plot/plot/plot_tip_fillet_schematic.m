function plot_tip_fillet_schematic()
% plot_tip_fillet_schematic -- SOFF(4.572mm) 的起算點與「倒圓段」處理方式示意
% =========================================================================
%   回答一個問題：SOFF 這條斜距**從哪裡開始量**、40um 倒圓那一段**算不算 / 怎麼算**。
%
%   幾何一律取自 pole_sensor_geometry（sensor 幾何唯一來源，CAD STEP 實測），
%   本檔不自帶任何幾何常數 —— 只做視覺化。以 P1（下極，beta=11.49deg）為例，
%   畫在**沿極軸的子午剖面**上（橫軸 s = 自極尖沿極軸、縱軸 r = 徑向）。
%
%   panel (a) 全景：SOFF 沿錐面母線走到貼附點，再沿外法線抬 AIR 到 sensor 中心。
%   panel (b) 極尖放大：倒圓球、切點 Q、虛擬錐頂 V、以及 SOFF 的兩段拆解
%                       s_ax = t_tan + (SOFF - t_tan)*cos(beta)
%             其中前段 t_tan 被當成**軸向**長度（非弧長）—— 這就是「怎麼處理圓弧」。
%             圖上同時標出該弧的**真實弧長** r_f*(pi/2-beta)，兩者差 22um。
%
%   風格①粗體框圖；純幾何示意圖**不放軸標題**（保留刻度數字）。
%   輸出 → figures/paper_fig/Section3_A/tip_fillet_schematic.png（覆蓋迭代）
% =========================================================================
    clc;
    here   = fileparts(fileparts(mfilename('fullpath')));            % → paper_fig_plot/
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section3_A');
    if ~exist(figdir,'dir'); mkdir(figdir); end

    % ---- 幾何：一律由唯一來源取得（不自帶常數）----
    CALROOT = fullfile(fileparts(fileparts(here)), 'matlab', 'Maxwell');   % here=paper_fig_plot → main/
    addpath(fullfile(CALROOT,'function'), fullfile(CALROOT,'utils'));
    cfg = model_config('long2016_hexapole_halfcut','tip40um');
    [~, ~, geo] = pole_sensor_geometry(cfg);

    K    = 1e3;                       % m → mm
    i    = 1;                         % P1（下極）
    rf   = geo.r_tip(i)*K;            % 倒圓半徑 0.040
    bet  = geo.beta(i);               % 半錐角 11.4916 deg
    tt   = geo.t_tan(i)*K;            % 倒圓的軸向推進 0.0322
    ev   = geo.apex_off(i)*K;         % 虛擬錐頂外移 0.1608
    sax  = geo.s_ax(i)*K;             % 軸向站位 4.4810
    Rr   = geo.R(i)*K;                % 該站位錐面半徑 0.9437
    SOFF = 4.572;                     % 沿貼附面斜距 [mm]（= 0.18 inch）
    AIR  = 0.41;                      % 離面氣隙 [mm]
    arc  = rf*(pi/2 - bet);           % 倒圓段的**真實弧長** 0.0548

    Q  = [tt, rf*cos(bet)];           % 錐面與倒圓的切點
    V  = [-ev, 0];                    % 虛擬錐頂
    Rs = @(s) (ev + s)*tan(bet);      % 真實錐面半徑
    nh = [-sin(bet), cos(bet)];       % 錐面外法線（2D 子午面）

    fprintf('beta=%.4f deg  r_f=%.4f  t_tan=%.4f  arc=%.4f  e_v=%.4f\n', bet*180/pi, rf, tt, arc, ev);
    fprintf('s_ax=%.4f  R=%.4f  |  [D]嚴格弧長 s_ax=%.4f (%+.1f um)  [C]自切點起算 s_ax=%.4f (%+.1f um)\n', ...
        sax, Rr, tt+(SOFF-arc)*cos(bet), (tt+(SOFF-arc)*cos(bet)-sax)*1e3, ...
        tt+SOFF*cos(bet), (tt+SOFF*cos(bet)-sax)*1e3);

    % ---- 樣式 ----
    FS = 36;  LW = 3.0;  FSA = 24;                 % 刻度字 / 線寬 / 標註字
    CST = [0.80 0.80 0.80];                        % 鋼件填色
    CARC= [0.85 0.10 0.10];                        % 倒圓段（紅）
    CSL = [0.05 0.10 0.95];                        % 錐面母線 / SOFF（藍）
    CVA = [0 0.55 0];                              % 虛擬錐頂輔助（綠）

    fig = figure('Color','w','Position',[60 60 1900 980]);
    tl  = tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact');
    tl.OuterPosition = [0 0.205 1 0.795];

    % =====================================================================
    % (a) 全景
    % =====================================================================
    ax = nexttile(tl,1);  hold(ax,'on');
    LV = 5.45;                                           % 只畫到 5.45mm（全錐 14.8mm 太長）
    s  = linspace(tt, LV, 200);
    fill(ax, [flip(s) s], [flip(Rs(s)) -Rs(s)], CST, 'EdgeColor','none');
    ph = linspace(pi/2+bet, 3*pi/2-bet, 80);             % 倒圓封口
    fill(ax, rf+rf*cos(ph), rf*sin(ph), CST, 'EdgeColor','none');
    plot(ax, s, Rs(s), '-', 'Color',[0 0 0], 'LineWidth',LW);
    plot(ax, s,-Rs(s), '-', 'Color',[0 0 0], 'LineWidth',LW);
    plot(ax, rf+rf*cos(ph), rf*sin(ph), '-', 'Color',[0 0 0], 'LineWidth',LW);
    plot(ax, [-0.4 LV], [0 0], '-.', 'Color',[0.45 0.45 0.45], 'LineWidth',1.6);   % 極軸

    % SOFF 沿母線（自切點 Q 起算的那一段）
    Pf = [sax, Rr];
    plot(ax, [Q(1) Pf(1)], [Q(2) Pf(2)], '-', 'Color',CSL, 'LineWidth',6);
    plot(ax, Q(1), Q(2), 'o', 'MarkerSize',11, 'MarkerFaceColor',CARC, 'MarkerEdgeColor','k','LineWidth',1.5);
    plot(ax, 0, 0, 'ko', 'MarkerFaceColor','k', 'MarkerSize',11);

    % sensor：自貼附點沿外法線抬 AIR
    Sc = Pf + AIR*nh;
    plot(ax, [Pf(1) Sc(1)], [Pf(2) Sc(2)], '-', 'Color',[0 0 0], 'LineWidth',2.5);
    rectangle(ax,'Position',[Sc(1)-0.30 Sc(2)-0.075 0.60 0.15], 'Curvature',0.15, ...
              'FaceColor',[1 0.85 0.25], 'EdgeColor','k', 'LineWidth',2.5);
    text(ax, Sc(1), Sc(2)+0.30, 'Hall sensor', 'Interpreter','latex', ...
         'FontSize',FSA, 'FontWeight','bold', 'HorizontalAlignment','center');
    text(ax, Sc(1)+0.40, Sc(2)-0.18, sprintf('$\\mathbf{gap\\;%.2f}$',AIR), 'Interpreter','latex', ...
         'FontSize',FSA, 'HorizontalAlignment','left');

    % 軸向站位 s_ax 的雙箭頭
    yb = -1.30;
    darrow(ax, [0 yb], [sax yb], [0 0 0], 2.5);
    plot(ax, [0 0], [0 yb], ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.5);
    plot(ax, [sax sax], [Rr yb], ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.5);
    text(ax, sax/2, yb-0.20, sprintf('$\\mathbf{s_{ax}=%.3f}$',sax), 'Interpreter','latex', ...
         'FontSize',FSA, 'HorizontalAlignment','center', 'VerticalAlignment','top');
    text(ax, (Q(1)+Pf(1))/2-0.35, (Q(2)+Pf(2))/2+0.42, ...
         sprintf('$\\mathbf{SOFF-t_{tan}=%.3f}$',SOFF-tt), 'Interpreter','latex', ...
         'FontSize',FSA, 'Color',CSL, 'HorizontalAlignment','center');
    text(ax, sax+0.12, Rr/2, sprintf('$\\mathbf{R=%.3f}$',Rr), 'Interpreter','latex', ...
         'FontSize',FSA, 'HorizontalAlignment','left');
    plot(ax, [sax sax], [0 Rr], '--', 'Color',[0.3 0.3 0.3], 'LineWidth',1.8);

    % 放大框
    zb = [-0.19 0.13 -0.095 0.115];
    plot(ax, zb([1 2 2 1 1]), zb([3 3 4 4 3]), '-', 'Color',CVA, 'LineWidth',2.5);
    text(ax, zb(2)+0.05, 0.42, '(b)', 'Interpreter','latex', 'FontSize',FSA, 'Color',CVA, ...
         'HorizontalAlignment','left');
    text(ax, -0.10, -0.50, 'Tip', 'Interpreter','latex', 'FontSize',FSA, ...
         'HorizontalAlignment','right');

    style2d(ax, [-0.45 5.45], [-1.75 1.75], [0 2 4], [-1 0 1], FS, LW);
    title(ax, '(a)', 'Interpreter','latex', 'FontSize',FS, 'FontWeight','bold');

    % =====================================================================
    % (b) 極尖放大
    % =====================================================================
    ax = nexttile(tl,2);  hold(ax,'on');
    s2 = linspace(tt, zb(2), 120);
    fill(ax, [flip(s2) s2], [flip(Rs(s2)) -Rs(s2)], CST, 'EdgeColor','none');
    fill(ax, rf+rf*cos(ph), rf*sin(ph), CST, 'EdgeColor','none');
    plot(ax, s2, Rs(s2), '-', 'Color',CSL, 'LineWidth',LW+1);
    plot(ax, s2,-Rs(s2), '-', 'Color',[0 0 0], 'LineWidth',LW);
    plot(ax, [-0.24 zb(2)], [0 0], '-.', 'Color',[0.45 0.45 0.45], 'LineWidth',1.6);

    % 倒圓段（上半）用紅色突顯 = 真實弧
    pa = linspace(pi, pi/2+bet, 60);
    plot(ax, rf+rf*cos(pa), rf*sin(pa), '-', 'Color',CARC, 'LineWidth',LW+3);
    pb = linspace(pi, 3*pi/2-bet, 60);
    plot(ax, rf+rf*cos(pb), rf*sin(pb), '-', 'Color',[0 0 0], 'LineWidth',LW);

    % 虛擬錐頂 + 母線延伸
    plot(ax, [V(1) Q(1)], [V(2) Q(2)], '--', 'Color',CVA, 'LineWidth',2.5);
    plot(ax, [V(1) Q(1)], [V(2) -Q(2)], '--', 'Color',CVA, 'LineWidth',2.5);
    plot(ax, V(1), V(2), 'o', 'MarkerSize',12, 'MarkerFaceColor','w', ...
         'MarkerEdgeColor',CVA, 'LineWidth',2.5);
    text(ax, V(1), -0.008, 'Virtual apex', 'Interpreter','latex', 'FontSize',FSA, ...
         'Color',CVA, 'HorizontalAlignment','center', 'VerticalAlignment','top');

    % 半錐角（標在角弧外、沿角平分線推出，避開弧線本身）
    aa = linspace(0, bet, 40);  ra = 0.088;
    plot(ax, V(1)+ra*cos(aa), ra*sin(aa), '-', 'Color',CVA, 'LineWidth',2);
    text(ax, V(1)+(ra+0.012)*cos(bet/2), (ra+0.012)*sin(bet/2)+0.006, ...
         sprintf('$\\mathbf{\\beta=%.2f^{\\circ}}$',bet*180/pi), ...
         'Interpreter','latex', 'FontSize',FSA, 'Color',CVA, ...
         'HorizontalAlignment','left', 'VerticalAlignment','bottom');

    % 關鍵點
    plot(ax, 0, 0, 'ko', 'MarkerFaceColor','k', 'MarkerSize',12);
    plot(ax, Q(1), Q(2), 'o', 'MarkerSize',12, 'MarkerFaceColor',CARC, ...
         'MarkerEdgeColor','k', 'LineWidth',1.5);
    plot(ax, [-0.052 -0.004], [0.050 0.004], '-', 'Color','k', 'LineWidth',1.5);   % Tip 引線
    text(ax, -0.056, 0.052, 'Tip', 'Interpreter','latex', 'FontSize',FSA, ...
         'HorizontalAlignment','right', 'VerticalAlignment','bottom');
    text(ax, Q(1)+0.010, Q(2)+0.008, 'Tangent point', 'Interpreter','latex', ...
         'FontSize',FSA, 'Color',CARC, 'HorizontalAlignment','left', 'VerticalAlignment','bottom');

    % t_tan（軸向）雙箭頭 —— 程式實際採用的推進量
    yt = -0.052;
    darrow(ax, [0 yt], [tt yt], CARC, 2.5);
    plot(ax, [tt tt], [Q(2) yt], ':', 'Color',CARC, 'LineWidth',1.5);
    plot(ax, [0 0], [0 yt], ':', 'Color',CARC, 'LineWidth',1.5);
    text(ax, tt+0.014, yt, sprintf('$\\mathbf{t_{tan}=%.4f}$ (axial, used)',tt), ...
         'Interpreter','latex', 'FontSize',FSA, 'Color',CARC, ...
         'HorizontalAlignment','left', 'VerticalAlignment','middle');

    % e_v 雙箭頭
    ye = -0.080;
    darrow(ax, [V(1) ye], [0 ye], CVA, 2.5);
    text(ax, V(1)/2, ye-0.005, sprintf('$\\mathbf{e_{v}=%.4f}$',ev), 'Interpreter','latex', ...
         'FontSize',FSA, 'Color',CVA, 'HorizontalAlignment','center', 'VerticalAlignment','top');

    % 真實弧長（紅弧）與後段母線
    text(ax, 0.004, 0.062, sprintf('$\\mathbf{true\\;arc=%.4f}$',arc), 'Interpreter','latex', ...
         'FontSize',FSA, 'Color',CARC, 'HorizontalAlignment','left');
    text(ax, 0.070, 0.098, '$\\mathbf{SOFF-t_{tan}\\;\\rightarrow}$', ...
         'Interpreter','latex', 'FontSize',FSA, 'Color',CSL, 'HorizontalAlignment','left');

    style2d(ax, zb(1:2), zb(3:4), [-0.1 0 0.1], [-0.05 0 0.05], FS, LW);
    title(ax, '(b)', 'Interpreter','latex', 'FontSize',FS, 'FontWeight','bold');

    % =====================================================================
    % 底部說明：SOFF 的拆解公式 + 三種可能定義的差
    % =====================================================================
    annotation(fig, 'textbox', [0.045 0.015 0.915 0.175], ...
        'String', { ...
        sprintf('$\\mathbf{s_{ax}=t_{tan}+(SOFF-t_{tan})\\cos\\beta=%.4f+%.4f\\times\\cos %.2f^{\\circ}=%.4f\\;mm}$', ...
                tt, SOFF-tt, bet*180/pi, sax), ...
        sprintf(['SOFF is measured {\\bf from the tip}; the fillet is counted as its {\\bf axial} advance ' ...
                 '$t_{tan}=%.4f$, not its arc length $%.4f$.'], tt, arc), ...
        sprintf(['Alternatives: true-arc first segment $\\rightarrow s_{ax}=%.4f$ (%+.1f $\\mu$m)  $\\vert$  ' ...
                 'measuring SOFF from the tangent point $\\rightarrow s_{ax}=%.4f$ (%+.1f $\\mu$m).  All dimensions in mm.'], ...
                tt+(SOFF-arc)*cos(bet), (tt+(SOFF-arc)*cos(bet)-sax)*1e3, ...
                tt+SOFF*cos(bet), (tt+SOFF*cos(bet)-sax)*1e3) }, ...
        'Interpreter','latex', 'FontSize',25, 'FitBoxToText','off', ...
        'HorizontalAlignment','left', 'VerticalAlignment','middle', ...
        'EdgeColor',[0 0 0], 'LineWidth',2.5, 'BackgroundColor',[1 1 1], 'Margin',10);

    out = fullfile(figdir, 'tip_fillet_schematic.png');
    exportgraphics(fig, out, 'Resolution', 150);
    fprintf('wrote %s\n', out);
    close(fig);
end

% ---- 風格①粗體框圖（2D；schematic 不放軸標題）----
function style2d(ax, XL, YL, XT, YT, FS, LW)
    axis(ax,'equal');
    xlim(ax,XL);  ylim(ax,YL);
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LW,'TickLength',[.018 .018], ...
           'XTick',XT,'YTick',YT,'Layer','top');
    box(ax,'on');  grid(ax,'off');
    hold(ax,'off');
end

% ---- 雙向箭頭（資料座標；等比軸故 x/y 同尺度）----
function darrow(ax, p1, p2, col, lw)
    plot(ax, [p1(1) p2(1)], [p1(2) p2(2)], '-', 'Color',col, 'LineWidth',lw);
    d = p2 - p1;  L = hypot(d(1),d(2));  if L == 0, return; end
    u = d/L;  n = [-u(2) u(1)];  h = 0.035*L;  w = 0.35*h;
    for q = [1 -1]
        b = (q>0)*p2 + (q<0)*p1;
        plot(ax, [b(1) b(1)-q*h*u(1)+w*n(1), b(1)-q*h*u(1)-w*n(1), b(1)], ...
                 [b(2) b(2)-q*h*u(2)+w*n(2), b(2)-q*h*u(2)-w*n(2), b(2)], ...
             '-', 'Color',col, 'LineWidth',lw);
    end
end
