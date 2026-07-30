function plot_p1p2_poles_3d()
% plot_p1p2_poles_3d -- long2016 P1+P2 磁路「側視圖」(x-z, y=0) + sensor + n+ + WP 十字
% =========================================================================
%   側視(x-z, y=0)單一 panel：P1(下極半切、水平+x)與 P2(上極傾斜全錐)共尖端於 WP(原點)。
%   磁路 = coil1(=P1)激發、每 (x,z) 格取最近 y=0 真實節點(不內插)、turbo 依 |B|(mT)、限 R=7mm 球內
%   （P3–P6 azimuth 60/120/240/300° 不在 y=0 面 → 自動只剩 P1&P2）。磁極 y=0 截面輪廓灰填+深邊。
%   colorbar 樣式同 circuit_side（style_cbar：粗體 + $\mathbf{|B|\;(mT)}$）。奇數等距 tick、無軸標題、字體 36。
%   輸出 → figures/paper_fig/Section3_A/p1p2_poles_3d.png（覆蓋迭代）。
% =========================================================================
    clc;
    here   = fileparts(mfilename('fullpath'));
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section3_A');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling\common_path');
    c = mt_constants();

    rf = c.POLE_TIP_R*1e3;  beta = atan2(c.POLE_R, c.POLE_CONE_LEN);
    L  = 8;  R_SPH = 7;                                   % 磁極畫長 / WP 球殼半徑 (mm)
    XL = [-8 8];  ZL = [-2.5 6.5];                        % x-z 視野
    XT = [-6 -3 0 3 6];  ZT = [-2 0 2 4 6];              % 奇數、等距、含原點、內縮
    FS = 36;
    GREEN=[0.10 0.60 0.20];  RED=[0.85 0.10 0.10];

    % ---- coil1 y=0 切面場（回 grid-sample 真實節點 + 縮放好的箭頭分量）----
    [Xs,Zs,Uq,Wq,Bm_mT,CLIM] = field_slice(c, XL, ZL, R_SPH);

    % ---- P1 / P2 的 y=0 截面輪廓 ----
    [p1x,p1z] = pole_outline(c, 1, rf, beta, L);
    [p2x,p2z] = pole_outline(c, 2, rf, beta, L);

    % ---- sensor / n+（tip40 定案，x-z 分量）----
    spP1=[4.8047;-1.6189];  nP1=[-0.1961;-0.9806];
    spP2=[-3.1454;3.9774];  nP2=[ 0.7420; 0.6704];

    fig = figure('Color','w','Position',[60 60 1280 760]);
    ax  = axes(fig);  hold(ax,'on');

    % 1) 磁極輪廓（灰填 + 深邊，畫在最底）
    patch(ax, p1x, p1z, [0.82 0.84 0.88], 'FaceAlpha',0.30, 'EdgeColor',[0.28 0.30 0.36], 'LineWidth',2.2);
    patch(ax, p2x, p2z, [0.82 0.84 0.88], 'FaceAlpha',0.30, 'EdgeColor',[0.28 0.30 0.36], 'LineWidth',2.2);

    % 2) 磁路箭頭（turbo 依 |B|mT 分 bin；線寬隨場漸粗）
    nb=28; edges=linspace(0,CLIM,nb+1); cmap=turbo(nb); lw=linspace(0.5,2.2,nb);
    for k=1:nb
        if k<nb, m=Bm_mT>=edges(k)&Bm_mT<edges(k+1); else, m=Bm_mT>=edges(k); end
        if any(m), quiver(ax,Xs(m),Zs(m),Uq(m),Wq(m),0,'Color',cmap(k,:),'LineWidth',lw(k),'MaxHeadSize',0.35); end
    end

    % 3) sensor 感測區域（實際尺寸圓柱側視矩形 2R×H，綠）+ n+ 紅箭頭
    draw_sensor2d(ax, spP1, nP1, 0.15, 0.10, GREEN);
    draw_sensor2d(ax, spP2, nP2, 0.15, 0.10, GREEN);
    draw_narrow2d(ax, spP1, nP1, 1.0, RED, 3.5);
    draw_narrow2d(ax, spP2, nP2, 1.0, RED, 3.5);

    % 4) WP（原點）十字：screen-space marker（MarkerSize22 / LineWidth3）
    plot(ax, 0, 0, '+', 'MarkerSize',22, 'Color','k', 'LineWidth',3);

    axis(ax,'equal');  xlim(ax,XL);  ylim(ax,ZL);
    box(ax,'on');  grid(ax,'off');
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',3.0,'TickLength',[.018 .018]);
    set(ax,'XTick',XT,'YTick',ZT);                       % 奇數等距、無軸標題（不設 xlabel/ylabel）
    colormap(ax,turbo);  clim(ax,[0 CLIM]);
    cb = colorbar(ax);  style_cbar(cb, FS);
    ax.Toolbar.Visible='off';  hold(ax,'off');

    out = fullfile(figdir,'p1p2_poles_3d.png');
    exportgraphics(fig, out, 'Resolution', 150);
    fprintf('wrote %s  (CLIM=%d mT)\n', out, CLIM);
end

% ============================================================================
function [Xs,Zs,Uq,Wq,Bm_mT,CLIM] = field_slice(c, XL, ZL, R_SPH)
% coil1(=P1) y=0 切面：每 (x,z) 格取最近 y=0 的真實節點(不內插)、限 R=7mm 球內。回箭頭分量(mm)。
    ddir = ansys_path('long2016_hexapole_halfcut','data','graded','coil1');
    d = import_ansys_data(ddir, 'all', 'coil1');
    fprintf('P1 = coil1 (graded/all): matched=%d, |B|max=%.4f T\n', numel(d.x), max(d.bsum));
    s = 1 - 2*c.pole_is_lower(1);                         % P1 lower → -1（全 source）
    zoff = -c.SPH_OFST*1e3;
    x=d.x*1e3; y=d.y*1e3; z=d.z*1e3+zoff;                 % mm, WP frame
    rWP = sqrt(x.^2+y.^2+z.^2);

    cell = 0.20;                                          % (x,z) 格邊 (mm)
    inwin = x>=XL(1)&x<=XL(2) & z>=ZL(1)&z<=ZL(2) & rWP<R_SPH;
    gi = find(inwin);
    xe = XL(1):cell:XL(2);  ze = ZL(1):cell:ZL(2);
    ix = discretize(x(gi),xe);  iz = discretize(z(gi),ze);
    ok = ~isnan(ix)&~isnan(iz);  gi=gi(ok); ix=ix(ok); iz=iz(ok);
    nz = numel(ze)-1;  cid=(ix-1)*nz+iz;  ay=abs(y(gi));
    [uc,~,g] = unique(cid);  sel=zeros(numel(uc),1);
    for k=1:numel(uc), idx=find(g==k); [~,b]=min(ay(idx)); sel(k)=gi(idx(b)); end
    keep = d.bsum(sel)>1e-4;  sel=sel(keep);

    Xs=x(sel);  Zs=z(sel);
    Bx=s*d.bx(sel);  Bz=s*d.bz(sel);  Bm=d.bsum(sel);
    arrow_max=0.30;  bmax=max(Bm);                        % 箭頭長度 ∝ |B|^0.25 (mm)
    bxz=hypot(Bx,Bz);  bxz(bxz==0)=1e-12;
    scl=arrow_max.*(Bm./bmax).^0.25./bxz;
    Uq=Bx.*scl;  Wq=Bz.*scl;  Bm_mT=Bm*1e3;
    CLIM=ceil(max(Bm_mT)/50)*50;
end

% ============================================================================
function [px, pz] = pole_outline(c, pidx, rf, beta, L)
% 磁極 y=0 截面外框(x,z)：下極水平半切平頂 / 上極傾斜全錐（龍飛模型）。
    th = c.pole_angles(pidx)*pi/180;
    if c.pole_is_lower(pidx)
        a=[cos(th);sin(th);0];  v=[0;0;-1];  half=true;
    else
        inc=c.upper_incline;  a=[cos(inc)*cos(th);cos(inc)*sin(th);sin(inc)];
        u=cross(a,[0;0;1]); u=u/norm(u);  v=cross(a,u);  half=false;
    end
    tip=[c.pole_tip_x(pidx); c.pole_tip_y(pidx); c.pole_tip_z_wp(pidx)]*1e3;
    ts=rf*(1+sin(beta));  rt=rf*cos(beta);
    psi=linspace(0,pi/2+beta,16);  nt=30;  t=linspace(ts,L,nt);
    axoff=[rf*(1-cos(psi)),  t(2:end)];
    radm =[rf*sin(psi),      rt+(t(2:end)-ts)*tan(beta)];
    edgeP = tip + a.*axoff + v.*radm;
    if half
        axline = tip + a.*axoff;  poly=[axline, fliplr(edgeP)];
    else
        edgeM = tip + a.*axoff - v.*radm;  poly=[edgeP, fliplr(edgeM)];
    end
    px=poly(1,:);  pz=poly(3,:);
end

% ============================================================================
function draw_sensor2d(ax, sp, n, R, H, col)
% sensor 圓柱側視矩形(x-z)：底面在 sp、沿 +n 長 H、半寬 R。綠填黑邊。
    sp=sp(:);  n=n(:)/norm(n);  td=[-n(2);n(1)];
    base=sp;  top=sp+H*n;
    rc=[base+R*td, base-R*td, top-R*td, top+R*td];
    patch(ax, rc(1,:), rc(2,:), col, 'EdgeColor','k', 'LineWidth',1.6);
end

% ============================================================================
function draw_narrow2d(ax, sp, n, Ln, col, lw)
% 紅 n+ 箭頭(x-z)：桿 + 小三角頭（對齊 n）。
    sp=sp(:);  n=n(:)/norm(n);  tip=sp+Ln*n;  td=[-n(2);n(1)];
    hl=0.32*Ln;  a=24*pi/180;
    d1=cos(a)*(-n)+sin(a)*td;  d2=cos(a)*(-n)-sin(a)*td;
    plot(ax,[sp(1) tip(1)],[sp(2) tip(2)],'-','Color',col,'LineWidth',lw);
    plot(ax,[tip(1) tip(1)+d1(1)*hl],[tip(2) tip(2)+d1(2)*hl],'-','Color',col,'LineWidth',lw);
    plot(ax,[tip(1) tip(1)+d2(1)*hl],[tip(2) tip(2)+d2(2)*hl],'-','Color',col,'LineWidth',lw);
end

% ============================================================================
function style_cbar(cb, FS)
% colorbar 標準樣式（同 circuit_side）：軸標題 LaTeX 數學粗體、刻度數字粗體。
    cb.FontSize=FS;  cb.FontWeight='bold';
    cb.Label.Interpreter='latex';  cb.Label.String='$\mathbf{|B|\;(mT)}$';  cb.Label.FontSize=FS;
end
