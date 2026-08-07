function plot_sensor_ring_schematic()
% plot_sensor_ring_schematic -- 「感測環」幾何示意（討論用；尚未做任何場積分）
% =========================================================================
%   在距極尖 s = 4.572mm(沿錐面斜距)處，用**垂直該極極軸的平面**截磁極，得到一條環；
%   環沿軸向給寬度 dX = 0.30mm(= sensor 圓柱直徑)，再沿周長切成小 tile。
%     P2 上磁極(完整錐) -> 整圈圓環 O 形，周長 2*pi*r
%     P1 下磁極(半切)   -> D 形環：下半圈錐面(pi*r) + 削平面上的一條弦(2r)
%
%   座標：**原點 (0,0,0) = 六極尖共球的球心**，圖上以黑點標示。
%     ⚠ 不使用 "WP" 這個字眼（使用者拍板 2026-08-06；見 .claude/rules/figure-style.md）。
%   [MODIFIED 2026-08-06] 反映真實方位（原本兩根都畫成水平局部框，未反映上極傾角）：
%     P1 下極：極軸**水平**徑向（azimuth 0deg），削平面 = 通過極軸的水平面，保留下半。
%     P2 上極：極軸自 WP 往上斜 **c.upper_incline = 36.59deg**（azimuth 180deg）→ 往 -x 且向上。
%   兩根畫在同一個 panel，傾角差異直接可見。
%
%   本檔**只畫幾何**：確認環的形狀 / 位置 / tile 切法，尚未取任何場值。
%
%   風格①粗體框圖 3D 變體 A（daspect[1 1 1] + box off + 手動框邊省最近角）。
% =========================================================================
    clc;
    here   = fileparts(fileparts(mfilename('fullpath')));
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section3_A');
    if ~exist(figdir,'dir'); mkdir(figdir); end

    % ---- 幾何常數（對齊 config/<model>/tip40um/mt_constants.m）----
    R_POLE = 3.0;    L_CONE = 15.0;
    % [MODIFIED 2026-08-07] 錐面 R(s) = s*tan(beta)，**無 R_TIP 偏移** —— 已用原始格點分類
    %   實測驗證（R = 0.0097 + 0.1985*s；s=1.5~6 誤差 < 0.01mm）。40um 倒圓只修圓極尖、
    %   不把側面往外推（相切球版會外移 0.0408mm，與實測不符）。
    RT2 = 0;
    R_NORM = 0.500;  RXY = R_NORM*sqrt(2/3);  RZ = R_NORM/sqrt(3);  % 極尖球座標
    INC    = upper_incline_deg();                                   % 36.59 deg
    beta   = atan2(R_POLE, L_CONE);                                 % 半錐角 11.31 deg
    SOFF   = 4.572;                                                 % 沿錐面斜距 [mm]
    x0     = SOFF*cos(beta);                                        % 截面的軸向位置 4.4832 mm
    r0     = RT2 + x0*tan(beta);                                    % 截面半徑 0.8966 mm
    LVIS   = 6.5;                                                   % 只畫到 6.5mm（全長 15mm 截面會太小）
    nA     = 18;                                                    % 示意用的徑向格線數（半圓）

    fprintf('beta=%.2f deg  upper_incline=%.2f deg\n', beta*180/pi, INC);
    fprintf('截面 @ 軸向 x0 = %.4f mm（= 沿錐面斜距 %.3f mm）, 半徑 R = %.4f mm\n', x0, SOFF, r0);
    fprintf('P2 整圓截面 : A = %.4f mm²\n', pi*r0^2);
    fprintf('P1 半圓截面 : A = %.4f mm²  （正好一半）\n', pi*r0^2/2);

    % ---- 兩根極的位姿（mm）----
    %   **原點 (0,0,0) = 六極尖共球的球心**（圖上以黑點標示）。
    P1 = make_pole([ RXY; 0; -RZ], [ 1; 0; 0]);                     % azimuth 0deg、下層、極軸水平 +x
    P2 = make_pole([-RXY; 0;  RZ], [-cosd(INC); 0; sind(INC)]);     % azimuth 180deg、上層、上仰 INC
    fprintf('P1 tip = (%.4f, 0, %.4f)   P2 tip = (%.4f, 0, %.4f) mm\n', ...
            P1.tip(1), P1.tip(3), P2.tip(1), P2.tip(3));

    FS  = 32;
    LWBOX = 4.0;   % [ADDED 2026-08-07 使用者拍板] 3D 框線寬（手動框邊與 ruler 共用；原 2.5/2 太細）
    GRY = [0.72 0.74 0.78];   RNG = [0.85 0.10 0.10];   SEN = [0.10 0.60 0.20];
    fig = figure('Color','w','Position',[60 60 1500 950]);  ax = axes(fig);  hold(ax,'on');

    % ---- P1 下極：半切（保留 v<0 側 = z 低於極軸）+ 削平面 + **半圓截面** ----
    draw_cone(ax, P1, beta, RT2, LVIS, pi, 2*pi, GRY, 0.28);
    draw_flat(ax, P1, beta, RT2, LVIS, GRY, 0.34);
    draw_section(ax, P1, beta, RT2, x0, pi, 2*pi, nA, RNG);

    % ---- P2 上極：完整錐 + **整圓截面** ----
    draw_cone(ax, P2, beta, RT2, LVIS, 0, 2*pi, GRY, 0.24);
    draw_section(ax, P2, beta, RT2, x0, 0, 2*pi, 2*nA, RNG);

    % ---- 極軸虛線（凸顯 36.59deg 傾角）----
    for p = [P1 P2]
        e = p.tip + p.ax*LVIS;
        plot3(ax, [p.tip(1) e(1)], [p.tip(2) e(2)], [p.tip(3) e(3)], 'k--', 'LineWidth',2);
    end
    plot3(ax, 0, 0, 0, 'ko', 'MarkerFaceColor','k', 'MarkerSize',10);   % 原點（六極尖共球球心）
    text(ax,  3.2, 0, -2.4, 'P1  (half-cut)', ...
         'FontSize',FS-4,'FontWeight','bold','HorizontalAlignment','center');
    text(ax, -3.6, 0,  4.4, sprintf('P2  (%.2f%c)', INC, char(176)), ...
         'FontSize',FS-4,'FontWeight','bold','HorizontalAlignment','center');

    XL = [-7 7];  YL = [-2.5 2.5];  ZL = [-2.8 5.2];
    xlim(ax,XL); ylim(ax,YL); zlim(ax,ZL);
    daspect(ax,[1 1 1]);  grid(ax,'off');  box(ax,'off');
    view(ax, -28, 14);
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LWBOX);   % ruler 線寬 = 手動框邊，避免粗細不一
    % [MODIFIED 2026-08-07] 刻度修正：
    %   Z：原 [-1 0 1 3] 是 4 個(偶數)且間距 1/1/2(不等距) → 改 [-1 1 3]（奇數 3 個、等距 2、兩端內縮）。
    set(ax,'XTick',[-4 0 4],'YTick',[-1 0 1],'ZTick',[-1 1 3]);
    % [ADDED 2026-08-07] 刻度數字**轉正**：MATLAB 3D 預設把標籤沿軸向旋轉，y 軸那組被轉得最斜、
    %   在此 view 下歪成一團。三軸一律 TickLabelRotation=0（水平），數字才讀得出來。
    ax.XAxis.TickLabelRotation = 0;
    ax.YAxis.TickLabelRotation = 0;
    ax.ZAxis.TickLabelRotation = 0;
    % [MODIFIED 2026-08-07 使用者拍板] 純幾何示意圖**不放軸標題**（xlabel/zlabel）；
    %   刻度數字保留。見 .claude/rules/figure-style.md「示意圖不放軸標題」。
    draw_box_edges3(ax, XL, YL, ZL, LWBOX);
    ax.Toolbar.Visible = 'off';  hold(ax,'off');

    out = fullfile(figdir, 'sensor_ring_schematic.png');
    exportgraphics(fig, out, 'Resolution', 120);
    fprintf('wrote %s\n', out);
end

% ---- 上極傾角（與 mt_constants 同式）----
function d = upper_incline_deg()
    PROT_H=7.0; YOKE_H=2.0; YOKE_MID_R=(84+106)/4; R_NORM=0.5;
    RXY=R_NORM*sqrt(2/3); RZ=R_NORM/sqrt(3);
    end_r = YOKE_MID_R - 11.5;   end_z = YOKE_H + PROT_H + 5;
    tip_z = -PROT_H - 6 + 2*RZ;
    d = atan2(end_z - tip_z, end_r - RXY)*180/pi;
end

% ---- 建一根極的位姿：tip / 極軸 ax / 垂直極軸的正交基 (u,v)；對 P1 v = +z（削平面法線）----
function p = make_pole(tip, ax)
    p.tip = tip;
    p.ax  = ax/norm(ax);
    p.u   = cross([0;0;1], p.ax);   p.u = p.u/norm(p.u);   % 水平、垂直極軸
    p.v   = cross(p.ax, p.u);
end

% ---- 錐面（phi 範圍可指定；半切只畫 v<0 半邊）----
function draw_cone(ax, p, beta, rt, L, p1, p2, col, alp)
    t = linspace(0, L, 60);   ph = linspace(p1, p2, 90);
    [TT, PP] = meshgrid(t, ph);   RR = rt + TT*tan(beta);
    X = p.tip(1) + p.ax(1)*TT + RR.*(p.u(1)*cos(PP) + p.v(1)*sin(PP));
    Y = p.tip(2) + p.ax(2)*TT + RR.*(p.u(2)*cos(PP) + p.v(2)*sin(PP));
    Z = p.tip(3) + p.ax(3)*TT + RR.*(p.u(3)*cos(PP) + p.v(3)*sin(PP));
    surf(ax, X, Y, Z, 'FaceColor',col, 'FaceAlpha',alp, 'EdgeColor','none');
end

% ---- 削平面（通過極軸、法線 = v）----
function draw_flat(ax, p, beta, rt, L, col, alp)
    t = linspace(0, L, 60);   r = rt + t*tan(beta);
    X = [p.tip(1)+p.ax(1)*t - p.u(1)*r; p.tip(1)+p.ax(1)*t + p.u(1)*r];
    Y = [p.tip(2)+p.ax(2)*t - p.u(2)*r; p.tip(2)+p.ax(2)*t + p.u(2)*r];
    Z = [p.tip(3)+p.ax(3)*t - p.u(3)*r; p.tip(3)+p.ax(3)*t + p.u(3)*r];
    surf(ax, X, Y, Z, 'FaceColor',col, 'FaceAlpha',alp, 'EdgeColor','none');
end

% ---- 磁極**內部**的截面（垂直極軸）：整圓（p1..p2 = 0..2pi）或半圓（pi..2pi）----
%   填色面 + 加粗輪廓 + 內部細格線示意「小格加總」。
function draw_section(ax, p, beta, rt, x0, p1, p2, n, col)
    R  = rt + x0*tan(beta);
    ph = linspace(p1, p2, 240);
    ring = p.tip + p.ax*x0 + R*(p.u*cos(ph) + p.v*sin(ph));
    if abs((p2-p1) - 2*pi) > 1e-9                       % 半圓：補上弦（削平面那條直邊）
        poly = [ring, ring(:,1)];
    else
        poly = ring;
    end
    patch(ax, poly(1,:), poly(2,:), poly(3,:), col, 'FaceAlpha',0.35, 'EdgeColor',col, 'LineWidth',3);
    % 內部細格線（示意小格；只畫少量，避免糊掉）
    for f = (1:4)/5
        pr = p.tip + p.ax*x0 + (R*f)*(p.u*cos(ph) + p.v*sin(ph));
        plot3(ax, pr(1,:), pr(2,:), pr(3,:), '-', 'Color',col, 'LineWidth',0.8);
    end
    pe = linspace(p1, p2, n+1);
    for k = 1:numel(pe)
        A = p.tip + p.ax*x0;   B = A + R*(p.u*cos(pe(k)) + p.v*sin(pe(k)));
        plot3(ax, [A(1) B(1)], [A(2) B(2)], [A(3) B(3)], '-', 'Color',col, 'LineWidth',0.8);
    end
end

% ---- 手動框邊：省掉離相機最近角相連的 3 邊（風格①變體 A）----
function draw_box_edges3(ax, XL, YL, ZL, lw)
    C = [XL([1 1 1 1 2 2 2 2]).', YL([1 1 2 2 1 1 2 2]).', ZL([1 2 1 2 1 2 1 2]).'];
    E = [1 2;1 3;1 5;2 4;2 6;3 4;3 7;4 8;5 6;5 7;6 8;7 8];
    cp = campos(ax);  [~,near] = min(sum((C - cp).^2, 2));
    for k = 1:size(E,1)
        if any(E(k,:) == near), continue; end
        plot3(ax, C(E(k,:),1), C(E(k,:),2), C(E(k,:),3), 'k-', 'LineWidth',lw);
    end
end
