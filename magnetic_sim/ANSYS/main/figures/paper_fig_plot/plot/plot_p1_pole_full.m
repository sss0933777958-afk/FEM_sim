function plot_p1_pole_full()
% plot_p1_pole_full -- long2016 P1「整根磁極」3D 輪廓(下極:水平軸、半切平頂;透明灰面)
% =========================================================================
%   只畫磁極幾何輪廓(圓角尖端 + 直錐 + 圓柱段),半切平頂(水平);不畫場、不畫虛線。
%   風格沿用 Section2_C/3_A:透明灰 surf + 手動 box 邊(省最近角)、mm、font36 粗、LineWidth3。
%   view(-60,-30);沿軸畫到軸向長度 L(預設 15mm);較長軸 5 tick、較短軸 3 tick(皆內縮不含端點)。
%   輸出 → figures/paper_fig/Section3_A/p1_pole_full.png(覆蓋迭代)。
% =========================================================================
    clc;
    here   = fileparts(fileparts(mfilename('fullpath')));
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section3_A');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    % [MODIFIED 2026-08-08] 脫離 backup（規則 no-backup-data）→ live 樹。
    CALROOT = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling';
    addpath(fullfile(CALROOT,'function'), fullfile(CALROOT,'utils'), fullfile(CALROOT,'common_path'));
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling\common_path');
    c = model_config('long2016_hexapole_halfcut','tip40um');

    % ---- P1 幾何(下極:水平軸 + 半切平頂;WP frame, mm)----
    th = c.pole_angles(1)*pi/180;                        % P1 → 0°
    tip = [c.pole_tip_x(1); c.pole_tip_y(1); c.pole_tip_z_wp(1)]*1e3;
    a = [cos(th); sin(th); 0];                           % 水平極軸 (+x)
    u = [-sin(th); cos(th); 0];                          % 面內⊥(平頂內、+y)
    v = [0; 0; -1];                                      % 進入保留半邊(向下)→ 平頂在上(z=tip_z)
    rf = c.POLE_TIP_R*1e3;  beta = atan2(c.POLE_R, c.POLE_CONE_LEN);  Rcyl = c.POLE_R*1e3;

    L = 8;                                               % 軸向長度(mm)(與六極圖一致)
    XL=[0 9]; YL=[-2 2]; ZL=[-2.5 0];                    % 圓整框(L=8)

    fig = figure('Color','w','Position',[60 40 1180 1000]);
    ax  = axes(fig);  hold(ax,'on');
    [X,Y,Z] = draw_pole_half(ax, tip,a,u,v, rf,beta,Rcyl, L, true, [0.72 0.74 0.78], 0.16);
    fprintf('P1 tip=(%.3f,%.3f,%.3f)  bbox x[%.2f %.2f] y[%.2f %.2f] z[%.2f %.2f]\n', ...
        tip, min(X(:)),max(X(:)), min(Y(:)),max(Y(:)), min(Z(:)),max(Z(:)));
    % ---- sensor 感測區域圓柱（tip400 config sensor_pos/n；paper P1）----
    draw_cyl(ax, [4.8047;0;-1.6189], [-0.1961;0;-0.9806], 0.15, 0.1, [0.10 0.60 0.20], 0.95);  % tip40 sensor;R0.15 H0.1、底面=表面外0.41mm、沿 n+ 長
    % ---- 磁路箭頭：raw graded 場、turbo 依 |B|(mT)、長度∝|B|、voxel 抽稀（同 figure-style「3D 磁路場箭頭」定案）----
    draw_circuit_arrows(ax, c, tip, a, u, v, rf, beta, Rcyl, L, XL, YL, ZL);
    % ---- sensor n+ 紅箭頭（方向 = 真實 n+；箭頭頭繞 n+ 自轉使兩翼面向相機、頭好辨識）----
    sc=[4.8047;0;-1.6189]; nrm=[-0.1961;0;-0.9806]; Ln=1.3;
    draw_narrow(ax, sc, nrm, Ln, [0.85 0.10 0.10], 4.0);

    grid(ax,'off'); box(ax,'off'); daspect(ax,[1 1 1]);
    xlim(ax,XL); ylim(ax,YL); zlim(ax,ZL);
    view(ax,-45,-20);  camlight(ax,'headlight'); lighting(ax,'gouraud'); material(ax,'dull');
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
function draw_circuit_arrows(ax, c, tip, axk, u, v, rf, beta, Rcyl, L, XL, YL, ZL)
% 載 coil1(=P1,下極) graded 場 → 內部鐵節點(raw voxel 抽稀) + 尖端射出(raw voxel)
%   + 磁極表面(錐面殼) 用「重心(barycentric)內插」在細格上加密查詢點；turbo 依 |B|、長度 ∝ |B|。
    here = fileparts(fileparts(mfilename('fullpath')));
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling\common_path');
    ddir = ansys_path('long2016_hexapole_halfcut','data','graded','coil1');
    d = import_ansys_data(ddir, 'all', 'coil1');
    fprintf('P1 = coil1 (graded/all): matched=%d, |B|max=%.4f T\n', numel(d.x), max(d.bsum));
    S = load(fullfile(here, 'data', 'steel_ids.mat'));  iron = ismember(d.node_id, S.steel_ids);

    s    = 1 - 2*c.pole_is_lower(1);                    % P1 lower → -1（翻成全 source）
    zoff = -c.SPH_OFST*1e3;
    x=d.x*1e3; y=d.y*1e3; z=d.z*1e3+zoff;               % WP frame, mm
    bx=s*d.bx; by=s*d.by; bz=s*d.bz; bmag=d.bsum;
    coneL = c.POLE_CONE_LEN*1e3;

    % --- P1 錐座標(軸向 t、離軸 rp、線性錐半徑 rcone) ---
    vv = [x-tip(1), y-tip(2), z-tip(3)];
    t  = vv*axk;  rp = sqrt(max(sum(vv.^2,2)-t.^2,0));  rn = sqrt(sum(vv.^2,2));
    rcone = rf + t*(Rcyl-rf)/coneL;
    inbox = x>=XL(1)&x<=XL(2) & y>=YL(1)&y<=YL(2) & z>=ZL(1)&z<=ZL(2);

    % --- 內部鐵(raw voxel 0.30) + 尖端射出扇(raw voxel 0.22) ---
    sel_i = voxpick3(inbox & bmag>1e-4 & iron, 0.30, x,y,z,bmag, XL,YL,ZL);
    sel_t = voxpick3(inbox & bmag>1e-4 & ~iron & t<0.30 & rn<0.45, 0.22, x,y,z,bmag, XL,YL,ZL);
    sr=[sel_i;sel_t];  Pr=[x(sr) y(sr) z(sr)];  Br=[bx(sr) by(sr) bz(sr)];

    % --- 磁極表面：真實貼面節點(自然、優先) + 重心內插(scatteredInterpolant linear)細格填空隙 → 整根均勻加密 ---
    zt = tip(3);                                                          % 平頂 z(半切在此以下=鋼側)
    shell = inbox & ~iron & bmag>1e-6 & rp>rcone-0.03 & rp<rcone+0.30 & t>0.4 & t<L & z<zt+0.05;
    P0=[x(shell) y(shell) z(shell)];  B0=[bx(shell) by(shell) bz(shell)];
    src = inbox & ~iron & bmag>1e-6 & rp<rcone+1.8;                       % 建三角化的 air 節點
    Fx=scatteredInterpolant(x(src),y(src),z(src),bx(src),'linear','none');
    Fy=Fx; Fy.Values=by(src);   Fz=Fx; Fz.Values=bz(src);
    NT=72;  NP=24;  delta=0.12;                                           % 細格(軸×方位) / 往外偏移(mm)。P1 半切 → phi 0..pi
    tt=linspace(0.4,L,NT);  pp=linspace(0.02*pi,0.98*pi,NP);
    [TT,PP]=ndgrid(tt,pp);  TT=TT(:);  PP=PP(:);
    rc=rf+TT*(Rcyl-rf)/coneL;
    radx=u(1)*cos(PP)+v(1)*sin(PP);  rady=u(2)*cos(PP)+v(2)*sin(PP);  radz=u(3)*cos(PP)+v(3)*sin(PP);
    Qx=tip(1)+axk(1)*TT+(rc+delta).*radx;  Qy=tip(2)+axk(2)*TT+(rc+delta).*rady;  Qz=tip(3)+axk(3)*TT+(rc+delta).*radz;
    BQ=[Fx(Qx,Qy,Qz) Fy(Qx,Qy,Qz) Fz(Qx,Qy,Qz)];
    okq=all(~isnan(BQ),2) & Qx>=XL(1)&Qx<=XL(2)&Qy>=YL(1)&Qy<=YL(2)&Qz>=ZL(1)&Qz<=ZL(2) & sqrt(sum(BQ.^2,2))>1e-6;
    Pq=[P0; Qx(okq) Qy(okq) Qz(okq)];  Bq=[B0; BQ(okq,:)];               % 真實在前
    vsz=0.22;  key=floor((Pq-min(Pq,[],1))/vsz);                          % 體素每格留一(真實優先、空隙用內插填)
    [~,iu]=unique(key,'rows','stable');  Pq=Pq(iu,:);  Bq=Bq(iu,:);

    % --- 合併 raw(內部+尖端) + interp(表面) 一起畫 ---
    Pall=[Pr; Pq];  Ball=[Br; Bq];
    tall=(Pall-tip.')*axk;  in=tall<=L;                                  % 超過磁極軸長度 L 的部分不畫箭頭
    Pall=Pall(in,:);  Ball=Ball(in,:);
    Bm=sqrt(sum(Ball.^2,2));  Bm_mT=Bm*1e3;  CLIM=ceil(max(Bm_mT)/50)*50;
    lmin=0.15; lmax=0.55; bmax=max(Bm);                 % 長度 ∝ |B|(壓縮 ^0.35、最大不誇張)
    len=lmin+(lmax-lmin).*(Bm./bmax).^0.35;   bn=Bm+eps;
    U=Ball(:,1)./bn.*len;  V=Ball(:,2)./bn.*len;  W=Ball(:,3)./bn.*len;
    nb=28; edc=linspace(0,CLIM,nb+1); cmap=turbo(nb);
    for k=1:nb
        if k<nb, m=Bm_mT>=edc(k)&Bm_mT<edc(k+1); else, m=Bm_mT>=edc(k); end
        if any(m), quiver3(ax,Pall(m,1),Pall(m,2),Pall(m,3),U(m),V(m),W(m),0,'Color',cmap(k,:),'LineWidth',1.4,'MaxHeadSize',0.5); end
    end
end

% ============================================================================
function sel = voxpick3(mask, vsz, x, y, z, bmag, XL, YL, ZL)
% voxel 降取樣：每個 vsz 立方格取 |B|max 的真實節點(不內插)。
    ex=XL(1):vsz:XL(2); if ex(end)<XL(2)-1e-9, ex=[ex XL(2)]; end
    ey=YL(1):vsz:YL(2); if ey(end)<YL(2)-1e-9, ey=[ey YL(2)]; end
    ez=ZL(1):vsz:ZL(2); if ez(end)<ZL(2)-1e-9, ez=[ez ZL(2)]; end
    ix=discretize(x,ex); iy=discretize(y,ey); iz=discretize(z,ez);
    good=mask&~isnan(ix)&~isnan(iy)&~isnan(iz);
    ncy=numel(ey)-1; ncz=numel(ez)-1;
    cid=(ix-1)*ncy*ncz+(iy-1)*ncz+iz;
    vi=find(good); if isempty(vi), sel=[]; return; end
    [~,~,g]=unique(cid(vi)); sel=zeros(max(g),1);
    for k=1:max(g), cn=vi(g==k); [~,b]=max(bmag(cn)); sel(k)=cn(b); end
end

% ============================================================================
function draw_narrow(ax, sc, n, Ln, col, lw)
% 紅色方向箭頭：桿 + 箭頭頭。方向 = 真實 n(不改)；箭頭頭兩翼取「⊥ n 且 ⊥ 視線」方向
%   → 頭繞 n 自轉面向相機、任何視角都看得清（解決 quiver3 平行視線時頭側對變一條線）。
    sc=sc(:);  n=n(:)/norm(n);  tip=sc+n*Ln;
    [az,el]=view(ax);
    vd=[sind(az)*cosd(el); -cosd(az)*cosd(el); sind(el)];   % 視線
    w=cross(n,vd);  if norm(w)<1e-6, w=cross(n,[0;0;1]); end;  if norm(w)<1e-6, w=cross(n,[0;1;0]); end
    w=w/norm(w);                                            % ⊥ n 且 ⊥ 視線 → 兩翼面向相機
    hl=0.32*Ln;  a=24*pi/180;  b=-n;
    d1=cos(a)*b+sin(a)*w;  d2=cos(a)*b-sin(a)*w;
    plot3(ax,[sc(1) tip(1)],[sc(2) tip(2)],[sc(3) tip(3)],'-','Color',col,'LineWidth',lw);
    plot3(ax,[tip(1) tip(1)+d1(1)*hl],[tip(2) tip(2)+d1(2)*hl],[tip(3) tip(3)+d1(3)*hl],'-','Color',col,'LineWidth',lw);
    plot3(ax,[tip(1) tip(1)+d2(1)*hl],[tip(2) tip(2)+d2(2)*hl],[tip(3) tip(3)+d2(3)*hl],'-','Color',col,'LineWidth',lw);
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
% 框線(使用者拍板 2026-07-28)：
%   仰角<0(從下往上看)：頂框 + 底框都不手動畫(有刻度數字的 x/y 底邊由 MATLAB 座標軸 ruler 提供，
%     沒刻度的底邊使用者不要)；只補畫「遠端垂直邊」(省最近 2 角)給 3D 立體感。
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
            if istop(k) || isbot(k), continue; end                               % 頂框+底框：手動都不畫
            if any(a==near2)||any(b==near2), continue; end                       % 垂直邊：省最近兩角
        else
            if a==near || b==near, continue; end                                 % 標準：省最近角 3 邊
        end
        plot3(ax,[p1(1) p2(1)],[p1(2) p2(2)],[p1(3) p2(3)],'k-','LineWidth',lw);
    end
end
