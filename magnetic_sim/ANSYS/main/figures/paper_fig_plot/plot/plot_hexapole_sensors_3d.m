function plot_hexapole_sensors_3d(layer)
% plot_hexapole_sensors_3d -- long2016 六極整根輪廓(實心) + 6 感測器圓柱(綠) + n+ 法向箭頭(紅)
% =========================================================================
%   6 極:下極(P1/P3/P6)水平半切、上極(P2/P4/P5)傾斜全錐;實心灰面。
%   sensor:tip40 幾何重算(tip + 4.572·e2 + (gap+0.41)·n̂;gap=rf(1−sinβ))、綠圓柱 R0.15/H0.1 沿 n+。
%   n+ 紅箭頭(quiver3)。WP 原點黑十字。view 仿 #220、daspect 立方、手動 box 邊、tick 奇數。
%   layer='upper' → 只畫上極 view(70,30)；'lower' → 只畫下極 view(70,-25);nargin<1 兩張都出。
%   輸出 → figures/paper_fig/Section3_A/hexapole_sensors_{upper,lower}_3d.png(覆蓋迭代)。
% =========================================================================
    if nargin < 1
        plot_hexapole_sensors_3d('upper');  plot_hexapole_sensors_3d('lower');  plot_hexapole_sensors_3d('merge');  return;
    end
    clc;
    here   = fileparts(fileparts(mfilename('fullpath')));
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section3_A');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    % [MODIFIED 2026-08-08] 脫離 backup（規則 no-backup-data）→ live config，
    %   否則 cnst 沒有 pole_cone_slope，pole_sensor_geometry 會靜默回退成名目 beta。
    CAL = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling';
    addpath(fullfile(CAL,'function'), fullfile(CAL,'utils'), fullfile(CAL,'common_path'));
    c = model_config('long2016_hexapole_halfcut','tip40um');

    % =============================== 合併(上|下)===============================
    if strcmp(layer,'merge')
        fig = figure('Color','w','Units','pixels','Position',[20 30 2200 1050]);
        t = tiledlayout(fig,1,2,'TileSpacing','loose','Padding','compact');
        ax1 = nexttile(t);  draw_layer(ax1, 'upper', c);
        ax2 = nexttile(t);  draw_layer(ax2, 'lower', c);
        out = fullfile(figdir,'hexapole_sensors_merged_3d.png');
        exportgraphics(fig, out, 'Resolution', 120);
        fprintf('wrote %s\n', out);  return;
    end

    % =============================== 單張 upper/lower ===========================
    fig = figure('Color','w','Position',[40 30 1250 1150]);  ax = axes(fig);
    draw_layer(ax, layer, c);
    fn = sprintf('hexapole_sensors_%s_3d.png', layer);
    out = fullfile(figdir,fn);
    exportgraphics(fig, out, 'Resolution', 120);
    fprintf('wrote %s\n', out);
end

% ============================================================================
function draw_layer(ax, layer, c)
% 把單一 layer(upper/lower)的 6-半極 + sensor + 框 + 視角 + tick 畫進 ax。
    inc  = c.upper_incline;
    rf   = c.POLE_TIP_R*1e3;  Rcyl = c.POLE_R*1e3;                 % mm
    SOFF = 4.572;                                                   % mm（沿貼附面斜面距離）
    % [MODIFIED 2026-08-08] sensor 幾何改由 utils/pole_sensor_geometry 供給（唯一來源）。
    %   舊寫法 sp = tip + SOFF*e2 + (gap+AIR)*nh（gap = rf(1-sin beta) 加在**法線**上）已作廢。
    [SPOS, SNOR, GEO] = pole_sensor_geometry(c, struct('soff_upper',SOFF*1e-3,'soff_lower',SOFF*1e-3));
    dir3 = @(el,az)[cos(el)*cos(az); cos(el)*sin(az); sin(el)];
    L    = 8;
    SEP  = 0.6;   % [USER] 圖示用:每根極沿自身軸往外推,分開中間匯聚的尖端(mm)
    colP = [0.62 0.65 0.70];  colS = [0.10 0.60 0.20];
    switch layer
        case 'upper', sel=~c.pole_is_lower(:); VW=[70 30];
        case 'lower', sel= c.pole_is_lower(:); VW=[70 -25];
    end

    hold(ax,'on');  allP = [];  ord = find(sel).';
    for k = 1:numel(ord)
        i = ord(k);
        th  = c.pole_angles(i)*pi/180;
        tip = [c.pole_tip_x(i); c.pole_tip_y(i); c.pole_tip_z_wp(i)]*1e3;   % mm
        if c.pole_is_lower(i)
            a=dir3(0,th); u=[-sin(th);cos(th);0]; v=[0;0;-1]; half=true;
        else
            a=dir3(inc,th); u=cross(a,[0;0;1]); u=u/norm(u); v=cross(a,u); half=false;
        end
        beta = GEO.beta(i);                      % per-pole 真實半錐角（CAD STEP）
        tip = tip + SEP*a;   % [USER] 沿自身軸外推,分開中間匯聚尖端(sensor sp 依 tip 一起移)
        [X,Y,Z] = draw_pole(ax, tip,a,u,v, rf,beta,Rcyl, L, half, colP, 1.0);
        allP = [allP; X(:) Y(:) Z(:)]; %#ok<AGROW>
        sp = SPOS(:,i)*1e3 + SEP*a;   sn = SNOR(:,i);   % 位置/法線由共用幾何供給（sp 同步 SEP 外推）
        draw_cyl(ax, sp, sn, 0.15, 0.10, colS, 0.95);
        allP = [allP; sp.']; %#ok<AGROW>
    end

    m = 0.6;  Nz = 5;  if strcmp(layer,'lower'), Nz = 3; end
    [XL,XT] = axis_setup(allP(:,1), 5, m);
    [YL,YT] = axis_setup(allP(:,2), 5, m);
    [ZL,ZT] = axis_setup(allP(:,3), Nz, m);
    grid(ax,'off'); box(ax,'off'); daspect(ax,[1 1 1]);
    xlim(ax,XL); ylim(ax,YL); zlim(ax,ZL);
    view(ax,VW(1),VW(2));  camlight(ax,'headlight'); lighting(ax,'gouraud'); material(ax,'dull');
    set(ax,'FontSize',36,'FontWeight','bold','LineWidth',3.0,'TickLength',[.015 .015]);
    set(ax,'XTick',XT,'YTick',YT,'ZTick',ZT);
    ax.XTickLabelRotation=0; ax.YTickLabelRotation=0; ax.ZTickLabelRotation=0;
    draw_box_edges3(ax, XL,YL,ZL, 3.0);
    ax.Clipping='off';  ax.Toolbar.Visible='off';  hold(ax,'off');
end

% ============================================================================
function [lim, tk] = axis_setup(vals, N, margin)
% 回傳軸範圍 + N 個「含 0、等距、nice step」的 tick(0-錨定,step 為 nice 值的倍數 → 必含 0)。
    emin = min(vals)-margin;  emax = max(vals)+margin;
    if emin > 0, emin = -0.06*(emax-emin); end          % 確保 0 落在範圍
    if emax < 0, emax =  0.06*(emax-emin); end
    span = emax - emin;
    P = 10^floor(log10(span/N));
    cand = sort([P/2 P P*1.5 P*2 P*2.5 P*3 P*5 P*10]);
    best = cand(1);  bestsc = inf;
    for s = cand
        n = floor(emax/s) - ceil(emin/s) + 1;           % s 的倍數落在範圍內的個數
        sc = abs(n-N)*2 + (ceil(emin/s)>0 || floor(emax/s)<0)*100;   % 命中 N + 必含 0
        if sc < bestsc, bestsc = sc; best = s; end
    end
    s = best;  tk = (ceil(emin/s):floor(emax/s))*s;     % 倍數(含 0)
    if numel(tk) > N                                     % 多了 → 以 0 為中心裁到 N
        [~,i0] = min(abs(tk));  lo = i0-floor((N-1)/2);  hi = lo+N-1;
        if lo<1, lo=1; hi=N; end
        if hi>numel(tk), hi=numel(tk); lo=hi-N+1; end
        tk = tk(lo:hi);
    end
    lim = [min(emin, tk(1)-0.35*s), max(emax, tk(end)+0.35*s)];
end

% ============================================================================
function draw_cyl(ax, cen, ndir, R, H, col, alp)
    ndir = ndir/norm(ndir);
    tmp = [1;0;0]; if abs(ndir.'*tmp)>0.9, tmp=[0;1;0]; end
    e1 = cross(ndir,tmp); e1=e1/norm(e1);  e2 = cross(ndir,e1);
    th = linspace(0,2*pi,40);  circ = e1*cos(th) + e2*sin(th);
    c1 = cen;  c2 = cen + H*ndir;
    P1 = c1 + R*circ;  P2 = c2 + R*circ;
    surf(ax, [P1(1,:);P2(1,:)],[P1(2,:);P2(2,:)],[P1(3,:);P2(3,:)], ...
         'FaceColor',col,'FaceAlpha',alp,'EdgeColor','none');
    patch(ax, P1(1,:),P1(2,:),P1(3,:), col,'FaceAlpha',alp,'EdgeColor','none');
    patch(ax, P2(1,:),P2(2,:),P2(3,:), col,'FaceAlpha',alp,'EdgeColor','none');
end

% ============================================================================
function [X,Y,Z] = draw_pole(ax, tip, a, u, v, rf, beta, Rcyl, L, half, col, alp)
% 磁極:圓角尖端+切線直錐(到 Rcyl)+圓柱段(到 L)；half=true→下極半切平頂。
    ts = rf*(1+sin(beta));  rt = rf*cos(beta);  L_cone = ts + (Rcyl-rt)/tan(beta);
    psi = linspace(0, pi/2+beta, 16);  axoff = rf*(1-cos(psi));  radm = rf*sin(psi);
    Lc = min(L, L_cone);  tc = linspace(ts, Lc, 40);
    axoff = [axoff, tc(2:end)];  radm = [radm, rt + (tc(2:end)-ts)*tan(beta)];
    if L > L_cone
        tcyl = linspace(L_cone, L, 16);  axoff = [axoff, tcyl(2:end)];  radm = [radm, Rcyl*ones(1,numel(tcyl)-1)];
    end
    if half, phi = linspace(0, pi, 40); else, phi = linspace(0, 2*pi, 64); end
    N = numel(axoff);  X=zeros(N,numel(phi)); Y=X; Z=X;
    for j = 1:numel(phi)
        rad = u*cos(phi(j)) + v*sin(phi(j));  P = tip + a*axoff + rad*radm;
        X(:,j)=P(1,:).'; Y(:,j)=P(2,:).'; Z(:,j)=P(3,:).';
    end
    surf(ax, X,Y,Z, 'FaceColor',col,'FaceAlpha',alp,'EdgeColor','none');
    if half
        eP = tip + a*axoff + u*radm;  eM = tip + a*axoff - u*radm;
        patch(ax, [eP(1,:) fliplr(eM(1,:))],[eP(2,:) fliplr(eM(2,:))],[eP(3,:) fliplr(eM(3,:))], col,'EdgeColor','none','FaceAlpha',alp);
        rimc = tip + a*axoff(end) + (u*cos(phi)+v*sin(phi))*radm(end);
        patch(ax, rimc(1,:),rimc(2,:),rimc(3,:), col,'EdgeColor','none','FaceAlpha',alp);
    else
        rimc = tip + a*axoff(end) + (u*cos(phi)+v*sin(phi))*radm(end);  cen = tip + a*axoff(end);
        patch(ax, [cen(1) rimc(1,:)],[cen(2) rimc(2,:)],[cen(3) rimc(3,:)], col,'EdgeColor','none','FaceAlpha',alp);
    end
end

% ============================================================================
function draw_box_edges3(ax, XL, YL, ZL, lw)
% 框線(使用者拍板 2026-07-28)：
%   仰角<0(從下往上看)：頂框不畫；底框只留近端 2 條(去掉最遠 2 條)；垂直邊只留遠端(省最近 2 角)。
%   仰角>=0：標準開口箱(12 邊省「最近角」相連 3 邊 → 9 邊)。
    [Xc,Yc,Zc]=ndgrid(XL,YL,ZL); C=[Xc(:) Yc(:) Zc(:)];
    E=[]; for i=1:8, for j=i+1:8
        if nnz(abs(C(i,:)-C(j,:))>1e-9)==1, E=[E; i j]; end
    end, end
    cp=campos(ax); dd=sum((C-cp).^2,2); [~,order]=sort(dd);  near=order(1); near2=order(1:2);
    [~,el]=view(ax); zmin=min(ZL); zmax=max(ZL);
    isbot=@(k) C(E(k,1),3)<zmin+1e-9 & C(E(k,2),3)<zmin+1e-9;
    istop=@(k) C(E(k,1),3)>zmax-1e-9 & C(E(k,2),3)>zmax-1e-9;
    for k=1:size(E,1)
        a=E(k,1); b=E(k,2); p1=C(a,:); p2=C(b,:);
        if el<0
            if istop(k) || isbot(k), continue; end                               % 頂框+底框：手動都不畫(刻度底邊由 ruler 提供)
            if any(a==near2)||any(b==near2), continue; end                       % 垂直邊：省最近兩角
        else
            if a==near || b==near, continue; end                                 % 標準：省最近角 3 邊
        end
        plot3(ax,[p1(1) p2(1)],[p1(2) p2(2)],[p1(3) p2(3)],'k-','LineWidth',lw);
    end
end
