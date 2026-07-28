function plot_p2_pole_full()
% plot_p2_pole_full -- long2016 P2「整根磁極」3D 輪廓(透明灰面,風格同 Section2_C 右圖)
% =========================================================================
%   只畫磁極幾何輪廓(圓角尖端 + 直錐 + 圓柱段),不畫場箭頭、不畫虛線圈。
%   風格沿用 plot_p2_charge_merged 右 3D panel:透明灰 surf + 手動 box 邊(省最近角)、
%   view(65,25)、mm、font36 粗、LineWidth3、各軸 5 內縮 tick 標數字。沿軸畫到軸向長度 L(預設 15.3mm，進圓柱段)。
%   輸出 → figures/paper_fig/Section3_A/p2_pole_full.png(覆蓋迭代)。
% =========================================================================
    clc;
    here   = fileparts(mfilename('fullpath'));
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section3_A');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling\common_path');
    c = mt_constants();

    % ---- P2 幾何(WP frame, mm)----
    inc=c.upper_incline; th=c.pole_angles(2)*pi/180; psi0=atan2(c.R_norm_z,c.R_norm_xy); Rn=c.R_norm*1e3;
    tip = Rn*[cos(psi0)*cos(th); cos(psi0)*sin(th); sin(psi0)];   % 尖端位置
    axk = [cos(inc)*cos(th); cos(inc)*sin(th); sin(inc)];         % 極軸(尖端→根部)
    u = cross(axk,[0;0;1]); u=u/norm(u);  v = cross(axk,u);       % 面內⊥ 基底
    rf = c.POLE_TIP_R*1e3;  beta = atan2(c.POLE_R, c.POLE_CONE_LEN);  Rcyl = c.POLE_R*1e3;

    L = 8;                                                         % 沿軸畫的軸向長度(mm)(與六極圖一致)

    % ---- 畫圖 ----
    fig = figure('Color','w','Position',[60 40 1180 1000]);
    ax  = axes(fig);  hold(ax,'on');
    [X,Y,Z] = draw_pole_full(ax, tip,axk,u,v, rf,beta,Rcyl, L, [0.72 0.74 0.78], 0.16); %#ok<ASGLU>
    % ---- sensor 感測區域圓柱（tip400 config sensor_pos/n；paper P2）----
    draw_cyl(ax, [-3.1454;0;3.9774], [0.7420;0;0.6704], 0.15, 0.1, [0.10 0.60 0.20], 0.95);  % tip40 sensor;R0.15 H0.1、底面=表面外0.41mm、沿 n+ 長

    XL=[-8.5 0]; YL=[-2 2]; ZL=[0 7];                              % 圓整框(L=8)
    grid(ax,'off'); box(ax,'off'); daspect(ax,[1 1 1]);
    xlim(ax,XL); ylim(ax,YL); zlim(ax,ZL);
    view(ax,65,25);  camlight(ax,'headlight'); lighting(ax,'gouraud'); material(ax,'dull');
    set(ax,'FontSize',36,'FontWeight','bold','LineWidth',3.0,'TickLength',[.02 .02]);
    set(ax,'XTick',[-7 -5 -3 -1],'YTick',[-1 0 1],'ZTick',[2 4 6]);  % 內縮 tick(L=8 重新框)
    ax.XTickLabelRotation=0; ax.YTickLabelRotation=0; ax.ZTickLabelRotation=0;        % 刻度數字不隨軸旋轉(保持水平)
    draw_box_edges3(ax, XL,YL,ZL, 3.0);
    ax.Clipping='off';  ax.Toolbar.Visible='off';  hold(ax,'off');

    out = fullfile(figdir,'p2_pole_full.png');
    exportgraphics(fig, out, 'Resolution', 130);
    fprintf('wrote %s  (L=%.2fmm, cone→cyl @ %.2fmm)\n', out, L, ...
        rf*(1+sin(beta)) + (Rcyl-rf*cos(beta))/tan(beta));
end

% ============================================================================
function draw_cyl(ax, cen, ndir, R, H, col, alp)
% 在 cen 畫一根圓柱(sensor 感測區域)：軸沿 ndir、半徑 R、高 H(沿軸 ±H/2)。單位 mm。
    ndir = ndir/norm(ndir);
    tmp = [1;0;0]; if abs(ndir.'*tmp)>0.9, tmp=[0;1;0]; end
    e1 = cross(ndir,tmp); e1=e1/norm(e1);  e2 = cross(ndir,e1);
    th = linspace(0,2*pi,40);  circ = e1*cos(th) + e2*sin(th);     % 3×40
    c1 = cen;  c2 = cen + H*ndir;                       % 底面在 cen、沿 +ndir(n+) 長 H
    P1 = c1 + R*circ;  P2 = c2 + R*circ;
    surf(ax, [P1(1,:);P2(1,:)],[P1(2,:);P2(2,:)],[P1(3,:);P2(3,:)], ...
         'FaceColor',col,'FaceAlpha',alp,'EdgeColor','none');            % 側面
    patch(ax, P1(1,:),P1(2,:),P1(3,:), col,'FaceAlpha',alp,'EdgeColor','none');   % 兩端蓋
    patch(ax, P2(1,:),P2(2,:),P2(3,:), col,'FaceAlpha',alp,'EdgeColor','none');
end

% ============================================================================
function [X,Y,Z] = draw_pole_full(ax, tip, a, u, v, rf, beta, Rcyl, L, col, alp)
% 透明磁極整根:圓角尖端 + 切線直錐(到半徑 Rcyl)+ 圓柱段(到軸向 L)。單位 mm。
    ts = rf*(1+sin(beta));  rt = rf*cos(beta);
    L_cone = ts + (Rcyl - rt)/tan(beta);            % 錐→柱轉換軸向位置
    % 圓角尖端
    psi = linspace(0, pi/2+beta, 16);
    axoff = rf*(1-cos(psi));   radm = rf*sin(psi);
    % 直錐(到 min(L,L_cone))
    Lc = min(L, L_cone);
    tc = linspace(ts, Lc, 44);
    axoff = [axoff, tc(2:end)];   radm = [radm, rt + (tc(2:end)-ts)*tan(beta)];
    % 圓柱段(若軸向超過錐段)
    if L > L_cone
        tcyl = linspace(L_cone, L, 20);
        axoff = [axoff, tcyl(2:end)];   radm = [radm, Rcyl*ones(1,numel(tcyl)-1)];
    end
    phi = linspace(0, 2*pi, 72);
    N = numel(axoff);  X=zeros(N,numel(phi)); Y=X; Z=X;
    for j = 1:numel(phi)
        rad = u*cos(phi(j)) + v*sin(phi(j));   P = tip + a*axoff + rad*radm;
        X(:,j)=P(1,:).'; Y(:,j)=P(2,:).'; Z(:,j)=P(3,:).';
    end
    surf(ax, X,Y,Z, 'FaceColor',col, 'FaceAlpha',alp, 'EdgeColor','none');
    rim = tip + a*axoff(end) + (u*cos(phi)+v*sin(phi))*radm(end);  cen = tip + a*axoff(end);
    patch(ax, [cen(1) rim(1,:)],[cen(2) rim(2,:)],[cen(3) rim(3,:)], col, 'EdgeColor','none','FaceAlpha',alp);
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
