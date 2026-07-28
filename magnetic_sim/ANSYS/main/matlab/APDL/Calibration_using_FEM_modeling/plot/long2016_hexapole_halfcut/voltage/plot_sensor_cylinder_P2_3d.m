%% plot_sensor_cylinder_P2_3d.m -- P2 sensor 實際圓柱 3D 示意（long2016 半切，真實 WP 框）
% P2(上極,方位180°,full cone)：sensor 在錐面 4.572mm，畫實際取樣圓柱(R=0.15mm、H=0.10mm 沿 n+)
% + n+ 箭頭；WP 中心用黑十字（無字）。輪廓仿 NTU flux_box 半透明綠風格。view(60,35.73)。
% 3D 風格 = 圖規則「變體 A」：daspect([1 1 1]) + box off + 手動 draw_box_edges + 三軸同刻度。
clear; clc;
here = fileparts(mfilename('fullpath'));
CAL  = fileparts(fileparts(here));                                   % ...\voltage_base
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');   % mt_constants
addpath(fullfile(CAL,'code','main_function'));                       % build_sensor_geometry
cnst = mt_constants();
outdir = fullfile(CAL,'figures','shared'); if ~exist(outdir,'dir'), mkdir(outdir); end
MM = 1e3;

MODE='shank';   % 'cone' = sensor 在錐面 4.572mm；'shank' = sensor 在圓柱 shank a=18mm
colP=[0.10 0.55 0.20]; colS=[0.80 0.25 0.15]; redn=[0.85 0.10 0.10];

% ---- P2 幾何（真實 WP 框）----
i2=2; th=cnst.pole_angles(i2)*pi/180;
psi0=atan2(cnst.R_norm_z,cnst.R_norm_xy); beta=atan2(cnst.POLE_R,cnst.POLE_CONE_LEN); inc=cnst.upper_incline;
dir=@(e,a)[cos(e)*cos(a);cos(e)*sin(a);sin(e)];
tip=(cnst.R_norm*dir(psi0,th))*MM;                 % 極尖 [mm]（magic-angle 位置）
ax_c=dir(inc,th);                                  % 極/柱軸 = FEM 實際傾角 inc_up（離 WP）
u=cross(ax_c,[0;0;1]); if norm(u)<1e-6, u=cross(ax_c,[0;1;0]); end
u=u/norm(u); v=cross(ax_c,u); v=v/norm(v);         % 錐/柱面兩正交基底（⊥ 軸）
CONE=cnst.POLE_CONE_LEN*MM; RSH=3.047; RadOff=3.047+0.41;   % 錐長15、柱半徑、離軸

% ---- sensor 位置/法線 ----
if strcmp(MODE,'shank')
    nhp=dir(inc+beta+pi/2,th); ur=nhp-dot(nhp,ax_c)*ax_c; ur=ur/norm(ur);   % 柱面徑向外法線（同 main.m shank_Vgi）
    a18=18; c2=tip+a18*ax_c+RadOff*ur; n2=ur;       % a=18mm 柱面 sensor [mm]
    Lcyl=19;                                        % 柱段畫到 19mm
    fn='sensor_cylinder_P2_shank18_3d.png';
    XL=[-16.5 1]; YL=[-3.6 3.6]; ZL=[-0.6 16.5]; TICK=5; qL=2.5; cr=1.1;
else
    [sp,sn]=build_sensor_geometry(cnst); c2=sp(:,i2)*MM; n2=sn(:,i2)/norm(sn(:,i2));   % 錐面 4.572mm
    fn='sensor_cylinder_P2_3d.png';
    XL=[-5.4 0.6]; YL=[-1.0 1.0]; ZL=[-0.3 4.5]; TICK=1; qL=0.8; cr=0.35;
end
Rsen=0.15; Hsen=0.10;                              % 實際圓柱 R/H [mm]
Lp=CONE; if strcmp(MODE,'cone'), Lp=5.5; end       % cone 模式只畫 5.5mm；shank 畫整根錐
Rc=Lp*tan(beta);

fig=figure('Color','w','Position',[100 100 1200 900]); axh=axes(fig); hold(axh,'on');

% ---- P2 錐面（半透明綠 + 深綠邊）----
ang=linspace(0,2*pi,60); s=linspace(0,Lp,32);
[SS,AA]=meshgrid(s,ang); Rr=SS*tan(beta);
CX=tip(1)+SS*ax_c(1)+Rr.*(cos(AA)*u(1)+sin(AA)*v(1));
CY=tip(2)+SS*ax_c(2)+Rr.*(cos(AA)*u(2)+sin(AA)*v(2));
CZ=tip(3)+SS*ax_c(3)+Rr.*(cos(AA)*u(3)+sin(AA)*v(3));
surf(axh,CX,CY,CZ,'FaceColor',colP,'FaceAlpha',0.13,'EdgeColor','none');
bc=tip+Lp*ax_c; rim=bc+Rc*(cos(ang).*u+sin(ang).*v);
plot3(axh,rim(1,:),rim(2,:),rim(3,:),'-','Color',colP*0.55,'LineWidth',1.3);
for a0=0:pi/3:2*pi-0.01
    ed=bc+Rc*(cos(a0)*u+sin(a0)*v);
    plot3(axh,[tip(1) ed(1)],[tip(2) ed(2)],[tip(3) ed(3)],'-','Color',colP*0.55,'LineWidth',1.0);
end
plot3(axh,tip(1),tip(2),tip(3),'s','MarkerSize',8,'MarkerFaceColor','k','MarkerEdgeColor','k');

% ---- shank 圓柱段（CONE..Lcyl，半徑 RSH）----
if strcmp(MODE,'shank')
    ss2=linspace(CONE,Lcyl,20); [SS2,AA2]=meshgrid(ss2,ang);
    SX=tip(1)+SS2*ax_c(1)+RSH*(cos(AA2)*u(1)+sin(AA2)*v(1));
    SY=tip(2)+SS2*ax_c(2)+RSH*(cos(AA2)*u(2)+sin(AA2)*v(2));
    SZ=tip(3)+SS2*ax_c(3)+RSH*(cos(AA2)*u(3)+sin(AA2)*v(3));
    surf(axh,SX,SY,SZ,'FaceColor',colP,'FaceAlpha',0.13,'EdgeColor','none');
    ec=tip+Lcyl*ax_c; cap=ec+RSH*(cos(ang).*u+sin(ang).*v);        % 端面圈
    plot3(axh,cap(1,:),cap(2,:),cap(3,:),'-','Color',colP*0.55,'LineWidth',1.3);
    patch(axh,cap(1,:),cap(2,:),cap(3,:),colP,'FaceAlpha',0.13,'EdgeColor','none');
end

% ---- sensor 實際圓柱（紅橘半透明）----
u2=cross(n2,[0;0;1]); if norm(u2)<1e-6, u2=cross(n2,[0;1;0]); end
u2=u2/norm(u2); v2=cross(n2,u2); v2=v2/norm(v2);
aa=linspace(0,2*pi,40);
ring=@(t) (c2+t*Hsen*n2) + Rsen*(cos(aa).*u2+sin(aa).*v2);   % t=0 底、t=1 頂（沿 n+）
R0=ring(0); R1=ring(1);
patch(axh,R0(1,:),R0(2,:),R0(3,:),colS,'FaceAlpha',0.55,'EdgeColor',colS*0.6,'LineWidth',1.0);
patch(axh,R1(1,:),R1(2,:),R1(3,:),colS,'FaceAlpha',0.55,'EdgeColor',colS*0.6,'LineWidth',1.0);
surf(axh,[R0(1,:);R1(1,:)],[R0(2,:);R1(2,:)],[R0(3,:);R1(3,:)],'FaceColor',colS,'FaceAlpha',0.55,'EdgeColor','none');

% ---- n+ 箭頭 ----
nt=c2+n2*qL;
quiver3(axh,c2(1),c2(2),c2(3),n2(1)*qL,n2(2)*qL,n2(3)*qL,0,'Color',redn,'LineWidth',2.6,'MaxHeadSize',2.2);
text(axh,nt(1)+0.03*diff(XL),nt(2),nt(3)+0.03*diff(ZL),'n_+','Color',redn,'FontWeight','bold','FontSize',15);

% ---- WP 十字（原點,黑,無字）----
plot3(axh,[-cr cr],[0 0],[0 0],'-k','LineWidth',2);
plot3(axh,[0 0],[-cr cr],[0 0],'-k','LineWidth',2);
plot3(axh,[0 0],[0 0],[-cr cr],'-k','LineWidth',2);

% ---- 3D 風格「變體 A」----
daspect(axh,[1 1 1]); view(axh,60,35.73); grid(axh,'off'); box(axh,'off');
xlim(axh,XL); ylim(axh,YL); zlim(axh,ZL);
set(axh,'FontSize',15,'FontWeight','bold','LineWidth',1.5,'TickLength',[.018 .018]);
tk=@(L)(ceil(L(1)/TICK)*TICK):TICK:(floor(L(2)/TICK)*TICK);
set(axh,'XTick',tk(XL),'YTick',tk(YL),'ZTick',tk(ZL));
xlabel(axh,'x (mm)','FontWeight','bold'); ylabel(axh,'y (mm)','FontWeight','bold'); zlabel(axh,'z (mm)','FontWeight','bold');
axh.Toolbar.Visible='off';
draw_box_edges(axh,XL,YL,ZL,1.5);

out=fullfile(outdir,fn); exportgraphics(fig,out,'Resolution',150);
fprintf('P2 sensor 中心=(%.3f,%.3f,%.3f) mm ; n+=(%.3f,%.3f,%.3f) ; saved %s\n', c2, n2, out);

%% ---- local: 手動框邊（省離相機最遠角相連 3 邊）----
function draw_box_edges(ax,xl,yl,zl,lw)
    [X,Y,Z]=ndgrid(xl,yl,zl); C=[X(:) Y(:) Z(:)];
    E=[]; for i=1:8, for j=i+1:8, if sum(C(i,:)~=C(j,:))==1, E=[E;i j]; end; end; end
    cp=campos(ax); [~,near]=min(vecnorm(C-cp,2,2));   % 省「離相機最近角」3 前緣 → 打開前面、不橫穿幾何（正仰角適用）
    for e=1:size(E,1)
        if any(E(e,:)==near), continue; end
        plot3(ax,C(E(e,:),1),C(E(e,:),2),C(E(e,:),3),'-','Color','k','LineWidth',lw);
    end
    xlim(ax,xl); ylim(ax,yl); zlim(ax,zl);
end
