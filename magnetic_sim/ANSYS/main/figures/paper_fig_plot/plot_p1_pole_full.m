function plot_p1_pole_full()
% plot_p1_pole_full -- long2016 P1「整根磁極」3D 輪廓(下極:水平軸、半切平頂;透明灰面)
% =========================================================================
%   只畫磁極幾何輪廓(圓角尖端 + 直錐 + 圓柱段),半切平頂(水平);不畫場、不畫虛線。
%   風格沿用 Section2_C/3_A:透明灰 surf + 手動 box 邊(省最近角)、mm、font36 粗、LineWidth3。
%   view(-60,-30);沿軸畫到軸向長度 L(預設 15mm);較長軸 5 tick、較短軸 3 tick(皆內縮不含端點)。
%   輸出 → figures/paper_fig/Section3_A/p1_pole_full.png(覆蓋迭代)。
% =========================================================================
    clc;
    here   = fileparts(mfilename('fullpath'));
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section3_A');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling\common_path');
    c = mt_constants();

    % ---- P1 幾何(下極:水平軸 + 半切平頂;WP frame, mm)----
    th = c.pole_angles(1)*pi/180;                        % P1 → 0°
    tip = [c.pole_tip_x(1); c.pole_tip_y(1); c.pole_tip_z_wp(1)]*1e3;
    a = [cos(th); sin(th); 0];                           % 水平極軸 (+x)
    u = [-sin(th); cos(th); 0];                          % 面內⊥(平頂內、+y)
    v = [0; 0; -1];                                      % 進入保留半邊(向下)→ 平頂在上(z=tip_z)
    rf = c.POLE_TIP_R*1e3;  beta = atan2(c.POLE_R, c.POLE_CONE_LEN);  Rcyl = c.POLE_R*1e3;

    L = 8;                                               % 軸向長度(mm)(與六極圖一致)

    fig = figure('Color','w','Position',[60 40 1180 1000]);
    ax  = axes(fig);  hold(ax,'on');
    [X,Y,Z] = draw_pole_half(ax, tip,a,u,v, rf,beta,Rcyl, L, true, [0.72 0.74 0.78], 0.16);
    fprintf('P1 tip=(%.3f,%.3f,%.3f)  bbox x[%.2f %.2f] y[%.2f %.2f] z[%.2f %.2f]\n', ...
        tip, min(X(:)),max(X(:)), min(Y(:)),max(Y(:)), min(Z(:)),max(Z(:)));
    % ---- sensor 感測區域圓柱（tip400 config sensor_pos/n；paper P1）----
    draw_cyl(ax, [4.8047;0;-1.6189], [-0.1961;0;-0.9806], 0.15, 0.1, [0.10 0.60 0.20], 0.95);  % tip40 sensor;R0.15 H0.1、底面=表面外0.41mm、沿 n+ 長

    XL=[0 9]; YL=[-2 2]; ZL=[-2.5 0];                    % 圓整框(L=8)
    grid(ax,'off'); box(ax,'off'); daspect(ax,[1 1 1]);
    xlim(ax,XL); ylim(ax,YL); zlim(ax,ZL);
    view(ax,-60,-20);  camlight(ax,'headlight'); lighting(ax,'gouraud'); material(ax,'dull');
    set(ax,'FontSize',36,'FontWeight','bold','LineWidth',3.0,'TickLength',[.02 .02]);
    set(ax,'XTick',[2 4 6 8],'YTick',[-1 0 1],'ZTick',[-2 -1 0]);  % 內縮 tick(L=8 重新框)
    ax.XTickLabelRotation=0; ax.YTickLabelRotation=0; ax.ZTickLabelRotation=0;
    draw_box_edges3(ax, XL,YL,ZL, 3.0);
    ax.Clipping='off';  ax.Toolbar.Visible='off';  hold(ax,'off');

    out = fullfile(figdir,'p1_pole_full.png');
    exportgraphics(fig, out, 'Resolution', 130);
    fprintf('wrote %s  (L=%.2fmm)\n', out, L);
end

% ============================================================================
function draw_cyl(ax, cen, ndir, R, H, col, alp)
% 在 cen 畫一根圓柱(sensor 感測區域)：軸沿 ndir、半徑 R、高 H(沿軸 ±H/2)。單位 mm。
    ndir = ndir/norm(ndir);
    tmp = [1;0;0]; if abs(ndir.'*tmp)>0.9, tmp=[0;1;0]; end
    e1 = cross(ndir,tmp); e1=e1/norm(e1);  e2 = cross(ndir,e1);
    th = linspace(0,2*pi,40);  circ = e1*cos(th) + e2*sin(th);
    c1 = cen;  c2 = cen + H*ndir;                       % 底面在 cen、沿 +ndir(n+) 長 H
    P1 = c1 + R*circ;  P2 = c2 + R*circ;
    surf(ax, [P1(1,:);P2(1,:)],[P1(2,:);P2(2,:)],[P1(3,:);P2(3,:)], ...
         'FaceColor',col,'FaceAlpha',alp,'EdgeColor','none');
    patch(ax, P1(1,:),P1(2,:),P1(3,:), col,'FaceAlpha',alp,'EdgeColor','none');
    patch(ax, P2(1,:),P2(2,:),P2(3,:), col,'FaceAlpha',alp,'EdgeColor','none');
end

% ============================================================================
function [X,Y,Z] = draw_pole_half(ax, tip, a, u, v, rf, beta, Rcyl, L, half, col, alp)
% 透明磁極:圓角尖端 + 切線直錐(到 Rcyl)+ 圓柱段(到軸向 L)。half=true → 半切平頂(平頂在 u-a 面)。
    ts = rf*(1+sin(beta));  rt = rf*cos(beta);
    L_cone = ts + (Rcyl - rt)/tan(beta);
    psi = linspace(0, pi/2+beta, 16);
    axoff = rf*(1-cos(psi));   radm = rf*sin(psi);
    Lc = min(L, L_cone);
    tc = linspace(ts, Lc, 44);
    axoff = [axoff, tc(2:end)];   radm = [radm, rt + (tc(2:end)-ts)*tan(beta)];
    if L > L_cone
        tcyl = linspace(L_cone, L, 20);
        axoff = [axoff, tcyl(2:end)];   radm = [radm, Rcyl*ones(1,numel(tcyl)-1)];
    end
    if half, phi = linspace(0, pi, 40); else, phi = linspace(0, 2*pi, 72); end
    N = numel(axoff);  X=zeros(N,numel(phi)); Y=X; Z=X;
    for j = 1:numel(phi)
        rad = u*cos(phi(j)) + v*sin(phi(j));   P = tip + a*axoff + rad*radm;
        X(:,j)=P(1,:).'; Y(:,j)=P(2,:).'; Z(:,j)=P(3,:).';
    end
    surf(ax, X,Y,Z, 'FaceColor',col, 'FaceAlpha',alp, 'EdgeColor','none');       % 曲面(半)
    if half
        % 平頂面(通過軸的直徑平面)：+u 邊 → -u 邊
        eP = tip + a*axoff + u*radm;   eM = tip + a*axoff - u*radm;
        patch(ax, [eP(1,:) fliplr(eM(1,:))],[eP(2,:) fliplr(eM(2,:))],[eP(3,:) fliplr(eM(3,:))], ...
              col,'EdgeColor','none','FaceAlpha',alp);
        % 末端半圓盤
        rimc = tip + a*axoff(end) + (u*cos(phi)+v*sin(phi))*radm(end);
        patch(ax, rimc(1,:),rimc(2,:),rimc(3,:), col,'EdgeColor','none','FaceAlpha',alp);
    else
        rimc = tip + a*axoff(end) + (u*cos(phi)+v*sin(phi))*radm(end);  cen = tip + a*axoff(end);
        patch(ax, [cen(1) rimc(1,:)],[cen(2) rimc(2,:)],[cen(3) rimc(3,:)], col,'EdgeColor','none','FaceAlpha',alp);
    end
end

% ============================================================================
function draw_box_edges3(ax, XL, YL, ZL, lw)
% 長方框 12 邊全黑等粗;省「離相機最近角」相連 3 邊(前面開口)→ 9 邊。
    [Xc,Yc,Zc]=ndgrid(XL,YL,ZL); C=[Xc(:) Yc(:) Zc(:)];
    E=[]; for i=1:8, for j=i+1:8
        if nnz(abs(C(i,:)-C(j,:))>1e-9)==1, E=[E; i j]; end
    end, end
    cp=campos(ax); dd=sum((C-cp).^2,2); [~,near]=min(dd);
    E=E(~(E(:,1)==near | E(:,2)==near), :);
    for k=1:size(E,1)
        p1=C(E(k,1),:); p2=C(E(k,2),:);
        plot3(ax,[p1(1) p2(1)],[p1(2) p2(2)],[p1(3) p2(3)],'k-','LineWidth',lw);
    end
end
