function plot_p2_pole_full()
% plot_p2_pole_full -- long2016 P2「整根磁極」3D 輪廓(透明灰面,風格同 Section2_C 右圖)
% =========================================================================
%   只畫磁極幾何輪廓(圓角尖端 + 直錐 + 圓柱段),不畫場箭頭、不畫虛線圈。
%   風格沿用 plot_p2_charge_merged 右 3D panel:透明灰 surf + 手動 box 邊(省最近角)、
%   view(65,25)、mm、font36 粗、LineWidth3、各軸 5 內縮 tick 標數字。沿軸畫到軸向長度 L(預設 15.3mm，進圓柱段)。
%   輸出 → figures/paper_fig/Section3_A/p2_pole_full.png(覆蓋迭代)。
% =========================================================================
    clc;
    here   = fileparts(fileparts(mfilename('fullpath')));
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

    XL=[-8.5 0]; YL=[-2 2]; ZL=[0 7];                              % 圓整框(L=8)

    % ---- 畫圖 ----
    fig = figure('Color','w','Position',[60 40 1180 1000]);
    ax  = axes(fig);  hold(ax,'on');
    [X,Y,Z] = draw_pole_full(ax, tip,axk,u,v, rf,beta,Rcyl, L, [0.72 0.74 0.78], 0.16); %#ok<ASGLU>
    % ---- sensor 感測區域圓柱（tip400 config sensor_pos/n；paper P2）----
    draw_cyl(ax, [-3.1454;0;3.9774], [0.7420;0;0.6704], 0.15, 0.1, [0.10 0.60 0.20], 0.95);  % tip40 sensor;R0.15 H0.1、底面=表面外0.41mm、沿 n+ 長
    % ---- 磁路箭頭：raw graded 場、turbo 依 |B|(mT)、長度∝|B|、voxel 抽稀（沿用 p2_charge_merged 右圖法）----
    draw_circuit_arrows(ax, c, tip, axk, u, v, rf, beta, Rcyl, L, XL, YL, ZL);
    % ---- sensor n+ 紅箭頭（方向 = 真實 n+；箭頭頭繞 n+ 自轉使兩翼面向相機、頭好辨識）----
    sc=[-3.1454;0;3.9774]; nrm=[0.7420;0;0.6704]; Ln=1.3;
    draw_narrow(ax, sc, nrm, Ln, [0.85 0.10 0.10], 4.0);

    grid(ax,'off'); box(ax,'off'); daspect(ax,[1 1 1]);
    xlim(ax,XL); ylim(ax,YL); zlim(ax,ZL);
    view(ax,45,25);  camlight(ax,'headlight'); lighting(ax,'gouraud'); material(ax,'dull');
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
function draw_circuit_arrows(ax, c, tip, axk, u, v, rf, beta, Rcyl, L, XL, YL, ZL)
% 載 coil5(=P2) graded 場 → 內部鐵(raw voxel) + 尖端射出(raw voxel)
%   + 磁極表面(全錐面) 真實貼面節點(自然、優先) + 重心內插細格填空隙 → 整根均勻加密。turbo 依 |B|、長度 ∝ |B|。
    here = fileparts(fileparts(mfilename('fullpath')));
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling\common_path');
    ddir = ansys_path('long2016_hexapole_halfcut','data','graded','coil5');
    d = import_ansys_data(ddir, 'all', 'coil5');
    fprintf('P2 = coil5 (graded/all): matched=%d, |B|max=%.4f T\n', numel(d.x), max(d.bsum));
    S = load(fullfile(here, 'data', 'steel_ids.mat'));  iron = ismember(d.node_id, S.steel_ids);

    s    = 1 - 2*c.pole_is_lower(2);                    % P2 upper → +1（全 source）
    zoff = -c.SPH_OFST*1e3;
    x=d.x*1e3; y=d.y*1e3; z=d.z*1e3+zoff;               % WP frame, mm
    bx=s*d.bx; by=s*d.by; bz=s*d.bz; bmag=d.bsum;
    coneL = c.POLE_CONE_LEN*1e3;

    vv = [x-tip(1), y-tip(2), z-tip(3)];
    t  = vv*axk;  rp = sqrt(max(sum(vv.^2,2)-t.^2,0));  rn = sqrt(sum(vv.^2,2));
    rcone = rf + t*(Rcyl-rf)/coneL;
    inbox = x>=XL(1)&x<=XL(2) & y>=YL(1)&y<=YL(2) & z>=ZL(1)&z<=ZL(2);

    % --- 內部鐵(raw voxel 0.30) + 尖端射出扇(raw voxel 0.22) ---
    sel_i = voxpick3(inbox & bmag>1e-4 & iron, 0.30, x,y,z,bmag, XL,YL,ZL);
    sel_t = voxpick3(inbox & bmag>1e-4 & ~iron & t<0.30 & rn<0.45, 0.22, x,y,z,bmag, XL,YL,ZL);
    sr=[sel_i;sel_t];  Pr=[x(sr) y(sr) z(sr)];  Br=[bx(sr) by(sr) bz(sr)];

    % --- 磁極表面(全錐)：真實貼面節點(優先) + 重心內插細格填空隙 → 整根均勻加密 ---
    shell = inbox & ~iron & bmag>1e-6 & rp>rcone-0.03 & rp<rcone+0.30 & t>0.4 & t<L;
    P0=[x(shell) y(shell) z(shell)];  B0=[bx(shell) by(shell) bz(shell)];
    src = inbox & ~iron & bmag>1e-6 & rp<rcone+1.8;
    Fx=scatteredInterpolant(x(src),y(src),z(src),bx(src),'linear','none');
    Fy=Fx; Fy.Values=by(src);   Fz=Fx; Fz.Values=bz(src);
    NT=72;  NP=44;  delta=0.12;                          % 細格(軸×方位，全錐 phi 0..2pi) / 往外偏移(mm)
    tt=linspace(0.4,L,NT);  pp=linspace(0,2*pi,NP+1);  pp(end)=[];
    [TT,PP]=ndgrid(tt,pp);  TT=TT(:);  PP=PP(:);
    rc=rf+TT*(Rcyl-rf)/coneL;
    radx=u(1)*cos(PP)+v(1)*sin(PP);  rady=u(2)*cos(PP)+v(2)*sin(PP);  radz=u(3)*cos(PP)+v(3)*sin(PP);
    Qx=tip(1)+axk(1)*TT+(rc+delta).*radx;  Qy=tip(2)+axk(2)*TT+(rc+delta).*rady;  Qz=tip(3)+axk(3)*TT+(rc+delta).*radz;
    BQ=[Fx(Qx,Qy,Qz) Fy(Qx,Qy,Qz) Fz(Qx,Qy,Qz)];
    okq=all(~isnan(BQ),2) & Qx>=XL(1)&Qx<=XL(2)&Qy>=YL(1)&Qy<=YL(2)&Qz>=ZL(1)&Qz<=ZL(2) & sqrt(sum(BQ.^2,2))>1e-6;
    Pq=[P0; Qx(okq) Qy(okq) Qz(okq)];  Bq=[B0; BQ(okq,:)];
    vsz=0.22;  key=floor((Pq-min(Pq,[],1))/vsz);
    [~,iu]=unique(key,'rows','stable');  Pq=Pq(iu,:);  Bq=Bq(iu,:);

    % --- 合併 raw(內部+尖端) + interp(表面) 一起畫 ---
    Pall=[Pr; Pq];  Ball=[Br; Bq];
    tall=(Pall-tip.')*axk;  in=tall<=L;                                  % 超過磁極軸長度 L 的部分不畫箭頭
    Pall=Pall(in,:);  Ball=Ball(in,:);
    Bm=sqrt(sum(Ball.^2,2));  Bm_mT=Bm*1e3;  CLIM=ceil(max(Bm_mT)/50)*50;
    lmin=0.15; lmax=0.55; bmax=max(Bm);                 % 長度 ∝ |B|(壓縮 ^0.35、最大不誇張)
    len=lmin+(lmax-lmin).*(Bm./bmax).^0.35;  bn=Bm+eps;
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
% 框線(使用者拍板 2026-07-28)：
%   仰角<0(從下往上看)：頂框不畫；底框只留「離相機近的 2 條邊」(去掉最遠 2 條)；垂直邊只留遠端(省最近 2 角)。
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
