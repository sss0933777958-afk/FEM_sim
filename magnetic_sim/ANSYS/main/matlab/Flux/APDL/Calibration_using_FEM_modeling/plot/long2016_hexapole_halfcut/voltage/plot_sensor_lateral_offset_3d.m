%% plot_sensor_lateral_offset_3d.m -- P1 sensor 偏移 3D 幾何示意圖（long2016 半切 P1，下極 θ=0）
% 產 3 張(figures/shared/)：
%   ('cone','lat' ,-30): 錐面 nominal(soff=4.572mm) 繞錐軸 -30° 側向偏移 → sensor_lateral_offset_P1_3d.png
%   ('cone','long',-2 ): 錐面 nominal 沿錐面滑 soff(-2mm,向 tip) 縱向偏移 → sensor_long_offset_P1_3d.png
%   ('shank','none')   : 圓柱 shank a=18mm，只畫幾何+nominal(無偏移)   → sensor_lateral_offset_P1_shank18_3d.png
% 極為藍半透明(仿 NTU 風格 FaceAlpha≈0.12+col*0.6 深邊)；偏移弧不加文字標註。
% 3D 風格 = 圖規則「變體 A」：daspect([1 1 1]) + box off + 手動 draw_box_edges(省最遠角 3 邊) + 三軸同刻度。
% 無場、無 colorbar。幾何來自 mt_constants + φ/滑動公式(同 main.m)。
clear; clc;
here = fileparts(mfilename('fullpath'));
CAL  = fileparts(fileparts(here));                            % ...\voltage_base
% [MODIFIED 2026-08-08] 脫離 backup（規則 no-backup-data）→ live config + utils/pole_sensor_geometry。
CALROOT = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling';
addpath(fullfile(CALROOT,'function'), fullfile(CALROOT,'utils'), fullfile(CALROOT,'common_path'));
cnst = model_config('long2016_hexapole_halfcut','tip40um');
outdir = fullfile(CALROOT,'figures','long2016_hexapole_halfcut','voltage','common'); if ~exist(outdir,'dir'), mkdir(outdir); end

render_offset('cone', 'lat',  -30, cnst, outdir);
render_offset('cone', 'long', -2,  cnst, outdir);
render_offset('shank','none',  0,  cnst, outdir);

%% ================= per-mode render =================
function render_offset(GEO, OFF, AMT, cnst, outdir)
    colP=[0.30 0.55 0.95]; green=[0.10 0.55 0.20]; redn=[0.85 0.10 0.10]; arcc=[0.85 0.45 0.05];
    % [MODIFIED 2026-08-08] sensor 幾何改由 utils/pole_sensor_geometry 供給（含 psi 繞極軸旋轉）；
    %   本檔原本的 e2c(p)/ncf(p) 是同一套 psi 慣例，但用舊構法（名目 beta、無真實錐體半徑外推）。
    ell=cnst.R_norm; psi0=atan2(cnst.R_norm_z,cnst.R_norm_xy);
    [~,~,gg] = pole_sensor_geometry(cnst);  beta = gg.beta(1);   % P1 真實半錐角（CAD）
    dirv=@(e,a)[cos(e)*cos(a);cos(e)*sin(a);sin(e)];
    MM=1e3; ah=[1;0;0]; vh=[0;1;0]; e1=dirv(-psi0,0); T=(ell*e1)*MM;
    CONE=gg.cone_len(1)*MM; RSH=3.047; Rc=CONE*tan(beta);      % 真實錐長（CAD）
    soffC=4.572e-3;   % AIR 由 pole_sensor_geometry 內部處理（預設 0.41mm）
    % psi 版位置/法線：直接取共用幾何（psi=0 即 canonical sensor 方位）
    Pg =@(p,so) sg_pos(cnst, so, p);   Ng =@(p,so) sg_nrm(cnst, so, p);
    RadOff=3.047+0.41; a0=18; nr0=[0;0;-1]; wf=@(p) cos(p)*nr0+sin(p)*vh;

    isShank=strcmp(GEO,'shank'); DRAW_OFF=~strcmp(OFF,'none'); Pf=[]; nf=[]; tt=[];
    if isShank
        Pnom=T(:)+a0*ah+RadOff*wf(0); n_nom=wf(0);
        Lcyl=22; Lp=CONE; XL=[-1 24]; YL=[-4.5 4.5]; ZL=[-5.6 0.8]; TICK=5;
        fn='sensor_lateral_offset_P1_shank18_3d.png'; ttl='P1 Hall sensor @ shank a=18mm（3D 示意）';
    else
        Pnom=Pg(0,soffC)*MM; n_nom=Ng(0,soffC);
        Lp=6; XL=[-0.6 6.4]; YL=[-2 2]; ZL=[-2.8 0.4]; TICK=2;
        if strcmp(OFF,'lat')
            tt=linspace(0,AMT,40)*pi/180;
            Pf=@(p)Pg(p,soffC)*MM;  nf=@(p)Ng(p,soffC);
            fn='sensor_lateral_offset_P1_3d.png'; ttl='P1 Hall sensor 側向偏移（3D 示意）';
        else % long：沿錐面滑 soff
            tt=linspace(0,AMT,40)*1e-3;
            Pf=@(d)Pg(0,soffC+d)*MM;  nf=@(d)Ng(0,soffC);
            fn='sensor_long_offset_P1_3d.png'; ttl='P1 Hall sensor 縱向偏移（3D 示意）';
        end
    end
    Pnom=Pnom(:);

    fig=figure('Color','w','Position',[100 100 1250 820]); ax=axes(fig); hold(ax,'on');

    % ---- P1 半切極（藍半透明）----
    ang=linspace(pi,2*pi,46);
    s=linspace(0,Lp,30); [SS,AA]=meshgrid(s,ang); Rr=SS*tan(beta);
    surf(ax,T(1)+SS,Rr.*cos(AA),T(3)+Rr.*sin(AA),'FaceColor',colP,'FaceAlpha',0.12,'EdgeColor','none');
    Rfar=Lp*tan(beta);
    if ~isShank
        patch(ax,[T(1) T(1)+Lp T(1)+Lp],[0 -Rfar Rfar],[T(3) T(3) T(3)],colP,'FaceAlpha',0.12,'EdgeColor','none');
        capx=T(1)+Lp; capR=Rfar; toppoly=[T(1) 0; T(1)+Lp -Rfar; T(1)+Lp Rfar];
    else
        xs=linspace(CONE,Lcyl,22); [XX,AA2]=meshgrid(xs,ang);
        surf(ax,T(1)+XX,RSH*cos(AA2),T(3)+RSH*sin(AA2),'FaceColor',colP,'FaceAlpha',0.12,'EdgeColor','none');
        patch(ax,[T(1) T(1)+CONE T(1)+Lcyl T(1)+Lcyl T(1)+CONE],[0 -Rc -RSH RSH Rc],T(3)*ones(1,5),colP,'FaceAlpha',0.12,'EdgeColor','none');
        capx=T(1)+Lcyl; capR=RSH; toppoly=[T(1) 0; T(1)+CONE -Rc; T(1)+Lcyl -RSH; T(1)+Lcyl RSH; T(1)+CONE Rc];
    end
    patch(ax,capx+0*ang,capR*cos(ang),T(3)+capR*sin(ang),colP,'FaceAlpha',0.12,'EdgeColor','none');
    plot3(ax,capx+0*ang,capR*cos(ang),T(3)+capR*sin(ang),'-','Color',colP*0.6,'LineWidth',1.2);
    plot3(ax,[toppoly(:,1);toppoly(1,1)],[toppoly(:,2);toppoly(1,2)],T(3)*ones(size(toppoly,1)+1,1),'-','Color',colP*0.6,'LineWidth',1.2);
    plot3(ax,[T(1) capx],[0 0],[T(3) T(3)-capR],'-','Color',colP*0.6,'LineWidth',1.2);

    % ---- WP + tip ----
    plot3(ax,0,0,0,'p','MarkerSize',17,'MarkerFaceColor',[1 .84 0],'MarkerEdgeColor','k');
    plot3(ax,[0 T(1)],[0 0],[0 T(3)],'-','Color',[.5 .5 .5],'LineWidth',1.5);
    plot3(ax,T(1),T(2),T(3),'s','MarkerSize',9,'MarkerFaceColor','k','MarkerEdgeColor','k');
    text(ax,-0.02*diff(XL),0,0.12,'WP','FontWeight','bold','FontSize',13,'HorizontalAlignment','right');
    text(ax,T(1),0,T(3)+0.05*abs(diff(ZL)),'  tip','FontWeight','bold','FontSize',12);
    text(ax,T(1)+0.35*Lp,0,T(3)-0.12*abs(diff(ZL)),'P1 pole','Color',colP*0.7,'FontWeight','bold','FontSize',13);

    % ---- nominal sensor（綠盤 + 紅 n+）----
    rd=0.05*diff(XL); qL=0.09*diff(XL); if isShank, qL=1.4; rd=1.1; end   % shank 箭頭/盤縮短防超界
    draw_disc(ax,Pnom,n_nom,rd,green,0.92,1.6);
    plot3(ax,Pnom(1),Pnom(2),Pnom(3),'o','MarkerSize',7,'MarkerFaceColor',green,'MarkerEdgeColor','k');
    quiver3(ax,Pnom(1),Pnom(2),Pnom(3),n_nom(1)*qL,n_nom(2)*qL,n_nom(3)*qL,0,'Color',redn,'LineWidth',2.6,'MaxHeadSize',2.2);
    nt=Pnom+n_nom*qL; text(ax,nt(1)-0.02*diff(XL),nt(2),nt(3)-0.03*abs(diff(ZL)),'n_+','Color',redn,'FontWeight','bold','FontSize',15,'HorizontalAlignment','right');
    text(ax,Pnom(1)-0.03*diff(XL),Pnom(2),Pnom(3)+0.055*abs(diff(ZL)),'sensor','Color',green*0.75,'FontWeight','bold','FontSize',12,'HorizontalAlignment','right');

    % ---- 偏移弧/線 + 端點淡化 sensor（無文字標註）----
    if DRAW_OFF
        arc=cell2mat(arrayfun(@(t)Pf(t),tt,'UniformOutput',false));
        plot3(ax,arc(1,:),arc(2,:),arc(3,:),'-','Color',arcc,'LineWidth',3.0);
        ce=Pf(tt(end)); ne=nf(tt(end));
        draw_disc(ax,ce,ne,rd*0.9,green,0.30,1.0);
        quiver3(ax,ce(1),ce(2),ce(3),ne(1)*qL*0.8,ne(2)*qL*0.8,ne(3)*qL*0.8,0,'Color',[.95 .55 .55],'LineWidth',1.6,'MaxHeadSize',2.2);
    end

    % ---- 3D 風格「變體 A」----
    daspect(ax,[1 1 1]); view(ax,30,-20); grid(ax,'off'); box(ax,'off');
    xlim(ax,XL); ylim(ax,YL); zlim(ax,ZL);
    set(ax,'FontSize',15,'FontWeight','bold','LineWidth',1.5,'TickLength',[.018 .018]);
    tk=@(L)(ceil(L(1)/TICK)*TICK):TICK:(floor(L(2)/TICK)*TICK);
    set(ax,'XTick',tk(XL),'YTick',tk(YL),'ZTick',tk(ZL));
    xlabel(ax,'x (mm)','FontWeight','bold'); ylabel(ax,'y (mm)','FontWeight','bold'); zlabel(ax,'z (mm)','FontWeight','bold');
    title(ax,ttl,'FontSize',15,'FontWeight','bold'); ax.Toolbar.Visible='off';
    draw_box_edges(ax,XL,YL,ZL,1.5);

    out=fullfile(outdir,fn); exportgraphics(fig,out,'Resolution',150);
    fprintf('[%s/%s] nominal=(%.3f,%.3f,%.3f) mm ; saved %s\n', GEO,OFF, Pnom(1),Pnom(2),Pnom(3), out);
end

%% ---- local: 3D 圓盤(⊥ n) ----
function draw_disc(ax,c,n,rd,col,alp,lw)
    c=c(:); n=n(:)/norm(n); t1=cross(n,[0;1;0]); if norm(t1)<1e-6, t1=cross(n,[1;0;0]); end
    t1=t1/norm(t1); t2=cross(n,t1); a=linspace(0,2*pi,44);
    P=c + rd*(t1*cos(a)+t2*sin(a));
    patch(ax,P(1,:),P(2,:),P(3,:),col,'FaceAlpha',alp,'EdgeColor',col*0.6,'LineWidth',lw);
end

%% ---- local: 手動框邊（矩形 box 12 邊，省離相機最遠角相連 3 邊 → 9 邊）----
function draw_box_edges(ax,xl,yl,zl,lw)
    [X,Y,Z]=ndgrid(xl,yl,zl); C=[X(:) Y(:) Z(:)];
    E=[]; for i=1:8, for j=i+1:8, if sum(C(i,:)~=C(j,:))==1, E=[E;i j]; end; end; end
    cp=campos(ax); [~,far]=max(vecnorm(C-cp,2,2));
    for e=1:size(E,1)
        if any(E(e,:)==far), continue; end
        plot3(ax,C(E(e,:),1),C(E(e,:),2),C(E(e,:),3),'-','Color','k','LineWidth',lw);
    end
    xlim(ax,xl); ylim(ax,yl); zlim(ax,zl);
end

%% ---- local：包一層 pole_sensor_geometry（P1、可帶 psi 與 soff）----
function P = sg_pos(cnst, soff, psi)
    S = pole_sensor_geometry(cnst, struct('soff_lower', soff, 'psi', psi));
    P = S(:,1);
end
function N = sg_nrm(cnst, soff, psi)
    [~, Nn] = pole_sensor_geometry(cnst, struct('soff_lower', soff, 'psi', psi));
    N = Nn(:,1);
end
