%% plot_pole_positions.m -- 極上候選 Hall sensor 位置概略圖（下極 P1 / 上極 P2 側視）
% =========================================================================
% x-z 子午面（WP 框, mm）；畫「錐段 + 圓柱段(shank)」外形 + 底面中心 + 3 個候選 sensor 位置：
%   錐面 soff=2.05mm（近 tip）/ 錐面 soff=4.572mm（nominal）/ 圓柱 a=18mm（shank）。
% 幾何全來自 mt_constants + build_sensor_geometry；錐半角 β、圓柱半徑 RSHANK、底面 a=36mm。
% 風格 ①粗體框圖。輸出實檔 → figures/shared/pole_sensor_positions.png（覆蓋迭代）。
% =========================================================================
clear; clc;
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');   % mt_constants
here = fileparts(mfilename('fullpath'));  TREE = fileparts(fileparts(here));   % .../voltage_base
addpath(fullfile(TREE,'code','main_function'));                                % build_sensor_geometry
cnst = mt_constants();
[sensor_pos,~] = build_sensor_geometry(cnst);        % nominal 4.572mm，WP 框 [m]

psi0 = atan2(cnst.R_norm_z, cnst.R_norm_xy);          % 仰角 ψ ≈ 35.26°
bta  = atan2(cnst.POLE_R, cnst.POLE_CONE_LEN);        % 半錐角 β ≈ 11.31°
inc  = cnst.upper_incline;                            % 上極錐軸傾角 ≈ 36.59°
Rn   = cnst.R_norm;  RSH = 3.047e-3;  AIR = 0.41e-3;  % 工作半徑 / 圓柱半徑 / standoff
CONE_A = cnst.POLE_CONE_LEN;  ABASE = 36e-3;          % 錐軸向長 15mm / 底面 a=36mm
dirv = @(e,a)[cos(e)*cos(a); cos(e)*sin(a); sin(e)];
q2   = @(v)[v(1); v(3)]*1e3;                          % 3D(WP,m) → 2D(x,z, mm)

col = [.95 .55 0; 0 .45 .2; .6 .1 .65; .1 .5 .75];   % 2.05 橘 / 4.572 綠 / shank18 紫 / shank25 藍

fig_dir = fullfile(TREE,'figures','shared');  if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
fig = figure('Color','w','Units','inches','Position',[1 1 12.5 5.8]);
tl  = tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact');

for panel = 1:2
    if panel==1, ip=1; low=true;  nm='P1 (lower, half-cut)';
    else,        ip=2; low=false; nm='P2 (upper)'; end
    th = cnst.pole_angles(ip)*pi/180;
    if low
        psi=-psi0; ax3=dirv(0,th);   e2=dirv(-bta,th);      nh=dirv(-bta-pi/2,th);
    else
        psi=+psi0; ax3=dirv(inc,th); e2=dirv(inc+bta,th);   nh=dirv(inc+bta+pi/2,th);
    end
    tip = Rn*dirv(psi,th);                               % 極尖 WP [m]
    nr  = nh - dot(nh,ax3)*ax3;  nr = nr/norm(nr);       % 圓柱面徑向外法線（去軸向分量）

    % --- 3 候選 sensor 位置（WP m）---
    s205 = tip + 2.05e-3*e2 + AIR*nh;                    % 錐面 soff=2.05
    s457 = sensor_pos(:,ip);                             % 錐面 soff=4.572（nominal）
    sshk = tip + 18e-3*ax3 + (RSH+AIR)*nr;               % 圓柱 a=18mm
    sshk2= tip + 25e-3*ax3 + (RSH+AIR)*nr;               % 圓柱 a=25mm
    basec= tip + ABASE*ax3;                              % 底面中心

    % --- 外形（2D x-z, mm）---
    T2=q2(tip); a2=[ax3(1);ax3(3)]; a2=a2/norm(a2); pp=[-a2(2);a2(1)]; R=RSH*1e3;
    Cend=T2+CONE_A*1e3*a2;  Base=T2+ABASE*1e3*a2;
    if low   % 半切：平頂(沿軸) + 底面 cap + 圓柱底 + 錐底斜邊
        poly=[T2, Base, Base-R*pp, Cend-R*pp];
    else     % 全錐+全圓柱：上下兩翼
        poly=[T2, Cend+R*pp, Base+R*pp, Base-R*pp, Cend-R*pp];
    end

    ax = nexttile(tl); hold(ax,'on');
    patch('Parent',ax,'XData',poly(1,:),'YData',poly(2,:), ...
          'FaceColor',[.82 .82 .86],'EdgeColor','k','LineWidth',2);
    plot(ax,[0 T2(1)],[0 T2(2)],'-','Color',[.6 .6 .6],'LineWidth',1.2);          % WP→tip
    plot(ax,0,0,'p','MarkerSize',15,'MarkerFaceColor',[1 .85 .1],'MarkerEdgeColor','k');
    text(ax,0,0,' WP','FontWeight','bold','FontSize',12);
    hb = plot(ax,Base(1),Base(2),'x','Color',[.15 .15 .15],'MarkerSize',13,'LineWidth',3);  % 底面中心
    P205=q2(s205); P457=q2(s457); PSHK=q2(sshk); PSHK2=q2(sshk2);
    h1=plot(ax,P205(1),P205(2),'o','MarkerSize',11,'MarkerFaceColor',col(1,:),'MarkerEdgeColor','k','LineWidth',1);
    h2=plot(ax,P457(1),P457(2),'o','MarkerSize',11,'MarkerFaceColor',col(2,:),'MarkerEdgeColor','k','LineWidth',1);
    h3=plot(ax,PSHK(1),PSHK(2),'o','MarkerSize',11,'MarkerFaceColor',col(3,:),'MarkerEdgeColor','k','LineWidth',1);
    h4=plot(ax,PSHK2(1),PSHK2(2),'o','MarkerSize',11,'MarkerFaceColor',col(4,:),'MarkerEdgeColor','k','LineWidth',1);
    hold(ax,'off');
    daspect(ax,[1 1 1]);
    title(ax,nm,'FontSize',14,'FontWeight','bold');
    xlabel(ax,'x (mm)','FontWeight','bold'); ylabel(ax,'z (mm)','FontWeight','bold');
    set(ax,'FontSize',15,'FontWeight','bold','LineWidth',2,'TickLength',[.02 .02]);
    box(ax,'on'); grid(ax,'off');
    xt=get(ax,'XTick'); set(ax,'XTick',xt(1:2:end));
    yt=get(ax,'YTick'); set(ax,'YTick',yt(1:2:end));
    ax.Toolbar.Visible='off';
    if panel==2
        lg=legend([h1 h2 h3 h4 hb], {'cone soff = 2.05 mm','cone soff = 4.572 mm (nominal)', ...
                  'shank a = 18 mm','shank a = 25 mm','base center'}, 'Location','northeast');
        lg.FontSize=12; lg.FontWeight='bold';
    end
end

out = fullfile(fig_dir,'pole_sensor_positions.png');
exportgraphics(fig, out, 'Resolution', 150);
fprintf('已輸出 %s\n', out);
