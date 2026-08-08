function plot_p2_charge_merged()
% plot_p2_charge_merged -- long2016 P2 尖端「電荷匯聚→射出」2D + 3D 合併(raw 真實節點)
% =========================================================================
%   左 = 2D 尖端俯視(p2_charge:匯聚→射出,x −0.5~−0.1 / y −0.2~0.2 mm)、
%   右 = 3D 磁極 + 尖端射出(p2_charge_3d:錐內匯聚柱 + 尖端放射狀射出)。
%   同激發 coil5=P2(上極 source s=+1)、'wp'、**共用 |B| 色階 + 單一 colorbar**。
%   規則:mm、無軸數字、粗框 LineWidth3、font 36、colorbar |B|(mT) LaTeX 粗體;
%   2D box on、3D box off + 手動 draw_box_edges(省最近角)。
%   輸出 → figures/paper_fig/Section2_C/p2_charge_merged.png(覆蓋迭代)。
% =========================================================================
    clc;
    here   = fileparts(fileparts(mfilename('fullpath')));
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section2_C');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    CAL = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling';
    addpath(fullfile(CAL,'function'));  addpath(fullfile(CAL,'common_path'));
    % [MODIFIED 2026-08-08] 脫離 backup（規則 no-backup-data）→ live 樹。
    CALROOT = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling';
    addpath(fullfile(CALROOT,'function'), fullfile(CALROOT,'utils'), fullfile(CALROOT,'common_path'));
    c = model_config('long2016_hexapole_halfcut','tip40um');

    cname='coil5';  pidx=2;
    ddir = ansys_path('long2016_hexapole_halfcut','data','graded',cname);
    d = import_ansys_data(ddir, 'wp', cname);
    fprintf('coil5=P2 (graded/wp): matched=%d, |B|max=%.4f T\n', numel(d.x), max(d.bsum));
    S = load(fullfile(here, 'data', 'steel_ids.mat'));  iron = ismember(d.node_id, S.steel_ids);

    s = 1 - 2*c.pole_is_lower(pidx);                 % P2 → +1
    bx=s*d.bx; by=s*d.by; bz=s*d.bz; bmag=d.bsum;
    CLIM = ceil(max(bmag)*1e3/50)*50;                % 共用色階(≈750 mT)
    fprintf('shared CLIM = %d mT\n', CLIM);

    fig = figure('Color','w','Position',[30 30 1760 900]);
    t = tiledlayout(fig,1,8,'TileSpacing','compact','Padding','compact');   % 無空白格,colorbar 緊接

    % ===== 左:2D 尖端俯視 =====
    ax1 = nexttile(t,1,[1 4]);  hold(ax1,'on');
    render_2d(ax1, bx,by,bmag, d, CLIM);
    ax1.Toolbar.Visible='off';  hold(ax1,'off');

    % ===== 右:3D 磁極 + 尖端射出 =====
    ax2 = nexttile(t,5,[1 4]);  hold(ax2,'on');
    render_3d(ax2, bx,by,bz,bmag, d, iron, c, CLIM);
    ax2.Toolbar.Visible='off';  hold(ax2,'off');

    % ===== 共用 colorbar =====
    colormap(fig,turbo);  clim(ax1,[0 CLIM]);  clim(ax2,[0 CLIM]);
    cb = colorbar(ax2);  cb.Layout.Tile='east';
    cb.FontSize=36;  cb.FontWeight='bold';
    cb.Label.Interpreter='latex';  cb.Label.String='$\mathbf{|B|\;(mT)}$';  cb.Label.FontSize=36;

    out = fullfile(figdir,'p2_charge_merged.png');
    exportgraphics(fig, out, 'Resolution', 105);
    fprintf('wrote %s\n', out);
end

% ============================================================================
function render_2d(ax, bx,by,bmag, d, CLIM)
% 2D 尖端俯視:每格取 |B|max 真實節點、|B|^0.25 長度、turbo 依 |B|(mT)。單位 mm。
    x=d.x*1e6; y=d.y*1e6; zwp=d.z*1e3+12.711;        % µm / mm
    XR=[-500 -100]; YR=[-200 200];                   % µm(zoom tip 前端)
    inbox = x>=XR(1)&x<=XR(2) & y>=YR(1)&y<=YR(2) & abs(zwp)<0.6;
    ncx=70; ncy=70;
    xe=linspace(XR(1),XR(2),ncx+1); ye=linspace(YR(1),YR(2),ncy+1);
    ix=discretize(x,xe); iy=discretize(y,ye);
    good=inbox&~isnan(ix)&~isnan(iy);
    cid=(ix-1)*ncy+iy; vi=find(good);
    [~,~,g]=unique(cid(vi)); sel=zeros(max(g),1);
    for k=1:max(g), cn=vi(g==k); [~,b]=max(bmag(cn)); sel(k)=cn(b); end
    keep=bmag(sel)>1e-4; sel=sel(keep);
    Xs=x(sel);Ys=y(sel);BXs=bx(sel);BYs=by(sel);Bm=bmag(sel); Bm_mT=Bm*1e3;
    arrow_max=15; bmax=max(Bm);
    bxy=hypot(BXs,BYs); bxy(bxy==0)=1e-12;
    scl=arrow_max.*(Bm./bmax).^0.25./bxy;
    Uq=BXs.*scl; Vq=BYs.*scl;
    Xs=Xs*1e-3; Ys=Ys*1e-3; Uq=Uq*1e-3; Vq=Vq*1e-3;  % mm
    nb=28; edges=linspace(0,CLIM,nb+1); cmap=turbo(nb); lw=linspace(0.5,2.0,nb);
    for k=1:nb
        if k<nb, m=Bm_mT>=edges(k)&Bm_mT<edges(k+1); else, m=Bm_mT>=edges(k); end
        if any(m), quiver(ax,Xs(m),Ys(m),Uq(m),Vq(m),0,'Color',cmap(k,:),'LineWidth',lw(k),'MaxHeadSize',0.3); end
    end
    axis(ax,'equal'); xlim(ax,[-0.5 -0.1]); ylim(ax,[-0.2 0.2]);
    box(ax,'on'); grid(ax,'off');
    set(ax,'FontSize',36,'FontWeight','bold','LineWidth',3.0,'TickLength',[.02 .02]);
    set(ax,'XTick',[-0.4 -0.2],'YTick',[-0.1 0 0.1],'XTickLabel',{},'YTickLabel',{});
end

% ============================================================================
function render_3d(ax, bx,by,bz,bmag, d, iron, c, CLIM)
% 3D:P2 透明磁極 + 錐內匯聚柱(鐵件)+ 尖端 apex 射出扇(air);raw 真實節點、voxel 降取樣。
    x=d.x*1e3; y=d.y*1e3; zwp=d.z*1e3+12.711;        % mm
    Rn=c.R_norm*1e3; r0=c.POLE_TIP_R*1e3;
    beta=atan2(c.POLE_R,c.POLE_CONE_LEN); inc=c.upper_incline; psi0=atan2(c.R_norm_z,c.R_norm_xy);
    th=c.pole_angles(2)*pi/180;
    axk=[cos(inc)*cos(th);cos(inc)*sin(th);sin(inc)];
    tip2=Rn*[cos(psi0)*cos(th);cos(psi0)*sin(th);sin(psi0)];
    u=cross(axk,[0;0;1]); u=u/norm(u); v=cross(axk,u);
    Lpole=0.60;
    XL=[-0.92 -0.34]; YL=[-0.20 0.20]; ZL=[0.14 0.70];
    inbox=x>=XL(1)&x<=XL(2)&y>=YL(1)&y<=YL(2)&zwp>=ZL(1)&zwp<=ZL(2);
    relx=x-tip2(1); rely=y-tip2(2); relz=zwp-tip2(3);
    proj=relx*axk(1)+rely*axk(2)+relz*axk(3);
    rnorm=sqrt(relx.^2+rely.^2+relz.^2);
    inb=inbox&bmag>1e-4&(iron|(~iron&proj<0.10&rnorm<0.10));   % 極內柱 + 尖端射出扇
    sel=voxpick3(inb,0.033,x,y,zwp,bmag,XL,YL,ZL);
    Xs=x(sel);Ys=y(sel);Zs=zwp(sel); BXs=bx(sel);BYs=by(sel);BZs=bz(sel); Bm=bmag(sel); Bm_mT=Bm*1e3;
    lmin=0.014; lmax=0.060; bmax=max(Bm);
    len=lmin+(lmax-lmin).*(Bm./bmax).^0.35;
    bn=sqrt(BXs.^2+BYs.^2+BZs.^2)+eps;
    Uq=BXs./bn.*len; Vq=BYs./bn.*len; Wq=BZs./bn.*len;
    draw_pole_t(ax, tip2,axk,u,v,r0,beta,Lpole,false,[0.72 0.74 0.78],0.14);
    nb=28; edc=linspace(0,CLIM,nb+1); cmap=turbo(nb);
    for k=1:nb
        if k<nb, m=Bm_mT>=edc(k)&Bm_mT<edc(k+1); else, m=Bm_mT>=edc(k); end
        if any(m), quiver3(ax,Xs(m),Ys(m),Zs(m),Uq(m),Vq(m),Wq(m),0,'Color',cmap(k,:),'LineWidth',1.4,'MaxHeadSize',0.5); end
    end
    % ---- 磁極內某段沿錐面輪廓的黑色粗虛線圈 ----
    projR = 0.35;  Rring = r0 + projR*tan(beta);  cen = tip2 + projR*axk;   % 該處錐半徑/中心
    thr = linspace(0,2*pi,120);  ring = cen + Rring*(u*cos(thr)+v*sin(thr));
    plot3(ax, ring(1,:),ring(2,:),ring(3,:), 'k--','LineWidth',5.0);
    grid(ax,'off'); box(ax,'off'); daspect(ax,[1 1 1]);
    xlim(ax,XL); ylim(ax,YL); zlim(ax,ZL);
    view(ax,65,25);  camlight(ax,'headlight'); lighting(ax,'gouraud'); material(ax,'dull');
    set(ax,'FontSize',36,'FontWeight','bold','LineWidth',3.0,'TickLength',[.02 .02]);
    set(ax,'XTick',[-0.8 -0.6 -0.4],'YTick',[-0.1 0 0.1],'ZTick',[0.2 0.4 0.6], ...
           'XTickLabel',{},'YTickLabel',{},'ZTickLabel',{});
    draw_box_edges3(ax, XL,YL,ZL, 3.0);
    ax.Clipping='off';
end

% ============================================================================
function sel = voxpick3(mask, vsz, x, y, z, bmag, XL, YL, ZL)
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

% ============================================================================
function draw_pole_t(ax, tip, a, u, v, rf, beta, L, half, col, alp)
% 透明磁極:rf 圓角 + 切線直錐 → L;half=true 半切平頂。單位 mm。
    ts=rf*(1+sin(beta)); rt=rf*cos(beta);
    psi=linspace(0,pi/2+beta,16); nt=30; t=linspace(ts,L,nt);
    axoff=[rf*(1-cos(psi)), t(2:end)]; radm=[rf*sin(psi), rt+(t(2:end)-ts)*tan(beta)];
    if half, phi=linspace(0,pi,44); else, phi=linspace(0,2*pi,64); end
    N=numel(axoff); X=zeros(N,numel(phi)); Y=X; Z=X;
    for j=1:numel(phi)
        rad=u*cos(phi(j))+v*sin(phi(j)); P=tip+a*axoff+rad*radm;
        X(:,j)=P(1,:).'; Y(:,j)=P(2,:).'; Z(:,j)=P(3,:).';
    end
    surf(ax,X,Y,Z,'FaceColor',col,'FaceAlpha',alp,'EdgeColor','none');
    rim=tip+a*axoff(end)+(u*cos(phi)+v*sin(phi))*radm(end); cen=tip+a*axoff(end);
    patch(ax,[cen(1) rim(1,:)],[cen(2) rim(2,:)],[cen(3) rim(3,:)],col,'EdgeColor','none','FaceAlpha',alp);
end
