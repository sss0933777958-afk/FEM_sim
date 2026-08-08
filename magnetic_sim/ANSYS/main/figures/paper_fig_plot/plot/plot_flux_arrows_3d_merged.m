function plot_flux_arrows_3d_merged()
% plot_flux_arrows_3d_merged -- long2016 半切六極 coil5=P2 3D 磁路箭頭:兩視角併一張
% =========================================================================
%   左 = view(65,25)、右 = view(105,25),同一 3D 場(graded / coil5=P2 上極 source s=+1)。
%   共用 |B| 色階 + 單一 colorbar(shared clim)。箭頭畫在磁極內(真實 FEM 鐵件節點 steel_ids)
%   + WP 尖端球;透明磁極輪廓(解析重建);3D voxel 降取樣(每體素 |B|max 真實節點,不內插)。
%   樣式:box on + daspect + µm(3D 粗體框 B 變體)、粗框、tick -500:500:500(奇數含 0 去端點)、
%   colorbar |B| (mT) 用 LaTeX 數學字型 + 粗體。P1..P6 per-view 置於清空區。
%   輸出 → figures/paper_fig/flux_arrows_3d_merged.png(覆蓋迭代)。
% =========================================================================
    clc;
    here   = fileparts(fileparts(mfilename('fullpath')));
    figdir = fullfile(fileparts(here), 'paper_fig');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    CAL = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling';
    addpath(fullfile(CAL,'function'));  addpath(fullfile(CAL,'common_path'));
    % [MODIFIED 2026-08-08] 脫離 backup（規則 no-backup-data）→ live 樹。
    CALROOT = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling';
    addpath(fullfile(CALROOT,'function'), fullfile(CALROOT,'utils'), fullfile(CALROOT,'common_path'));
    c = model_config('long2016_hexapole_halfcut','tip40um');

    COIL = 5;  cname = sprintf('coil%d',COIL);
    apdl_to_paper_idx = [1,3,6,5,2,4];  pidx = apdl_to_paper_idx(COIL);
    pname = sprintf('P%d', pidx);
    BOX = 900;                                         % 半視窗 [µm]

    % ---- 讀 graded WP 場 + 指紋 ----
    ddir = ansys_path('long2016_hexapole_halfcut','data','graded',cname);
    d = import_ansys_data(ddir, 'wp', cname);
    fprintf('\n==== %s = %s (3D merged) ====\n  matched=%d, |B|max=%.4f T\n', cname, pname, numel(d.x), max(d.bsum));

    % ---- 真實鐵件節點(mesh MAT2)----
    S = load(fullfile(here, 'data', 'steel_ids.mat'));
    iron = ismember(d.node_id, S.steel_ids);
    fprintf('  iron nodes(real): %d\n', nnz(iron));

    % ---- source 版 + WP 框座標(µm)----
    s = 1 - 2*c.pole_is_lower(pidx);
    bx = s*d.bx;  by = s*d.by;  bz = s*d.bz;
    bmag = d.bsum;
    x = d.x*1e6;  y = d.y*1e6;  zwp = d.z*1e6 + 12711;

    inbox = abs(x)<BOX & abs(y)<BOX & abs(zwp)<BOX;
    RWP = 550;  r3d = sqrt(x.^2+y.^2+zwp.^2);
    inb = inbox & (iron | r3d<RWP);

    sel = voxpick(inb, 70, x,y,zwp,bmag,BOX);
    Xs=x(sel);Ys=y(sel);Zs=zwp(sel); BXs=bx(sel);BYs=by(sel);BZs=bz(sel); Bm=bmag(sel);
    Bm_mT = Bm*1e3;
    fprintf('  quiver3 nodes(iron + WP): %d\n', numel(Xs));

    lmin=18; lmax=80;  bmax=max(Bm);
    len = lmin + (lmax-lmin).*(Bm./bmax).^0.35;
    bn = sqrt(BXs.^2+BYs.^2+BZs.^2)+eps;
    Uq=BXs./bn.*len; Vq=BYs./bn.*len; Wq=BZs./bn.*len;

    CLIM = ceil(max(Bm_mT)/50)*50;
    fprintf('  shared CLIM = %d mT\n', CLIM);
    % µm → mm(位置+箭頭同除，比例不變)
    SCN = struct('X',Xs/1e3,'Y',Ys/1e3,'Z',Zs/1e3,'U',Uq/1e3,'V',Vq/1e3,'W',Wq/1e3, ...
                 'Bm_mT',Bm_mT,'CLIM',CLIM,'BOX',BOX/1e3);

    FS = 36;
    fig = figure('Color','w','Position',[30 30 1900 940]);
    % 9 欄格:兩圖各佔 4 欄(1-4、5-8),第 9 欄留空 = 與 colorbar 隔開;colorbar 照原本 east 自動樣式
    t = tiledlayout(fig,1,9,'TileSpacing','compact','Padding','compact');

    % ---- 左:view(65,25)（tiles 1-4、不標 P、無軸標題）----
    ax1 = nexttile(t,1,[1 4]);  hold(ax1,'on');
    draw_scene(ax1, c, SCN, FS);
    view(ax1,65,25);  hold(ax1,'off');

    % ---- 右:view(105,25)（tiles 5-8、不標 P、無軸標題）----
    ax2 = nexttile(t,5,[1 4]);  hold(ax2,'on');
    draw_scene(ax2, c, SCN, FS);
    view(ax2,105,25);  hold(ax2,'off');

    % ---- 單一共用 colorbar ----
    colormap(fig,turbo);  clim(ax1,[0 CLIM]);  clim(ax2,[0 CLIM]);
    cb = colorbar(ax2);  cb.Layout.Tile='east';        % 原本自動樣式(不改),靠右邊空白格與圖隔開
    cb.FontSize=36;  cb.FontWeight='bold';
    cb.Label.Interpreter='latex';  cb.Label.String='$\mathbf{|B|\;(mT)}$';  cb.Label.FontSize=36;

    out = fullfile(figdir,'flux_arrows_3d_merged.png');
    exportgraphics(fig, out, 'Resolution', 100);
    fprintf('wrote %s\n', out);
end

% ============================================================================
function draw_scene(ax, c, S, FS)
% 把整個 3D 場景畫進指定 ax:透明磁極 + 箭頭 + WP 十字 + 框/刻度/軸標題。
    draw_poles_transparent(ax, c);
    nb=28; edc=linspace(0,S.CLIM,nb+1); cmap=turbo(nb);
    for k=1:nb
        if k<nb, m=S.Bm_mT>=edc(k)&S.Bm_mT<edc(k+1); else, m=S.Bm_mT>=edc(k); end
        if any(m), quiver3(ax,S.X(m),S.Y(m),S.Z(m),S.U(m),S.V(m),S.W(m),0, ...
                          'Color',cmap(k,:),'LineWidth',1.2,'MaxHeadSize',0.5); end
    end
    plot3(ax,0,0,0,'k+','MarkerSize',14,'LineWidth',2.4);
    grid(ax,'off'); box(ax,'on'); daspect(ax,[1 1 1]);
    xlim(ax,[-S.BOX S.BOX]); ylim(ax,[-S.BOX S.BOX]); zlim(ax,[-S.BOX S.BOX]);
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',3.0,'TickLength',[.018 .018]);
    set(ax,'XTick',[-0.5 0 0.5],'YTick',[-0.5 0 0.5],'ZTick',[-0.5 0 0.5]);   % mm、奇數含 0 去端點
    ax.Toolbar.Visible='off';                                                  % 無 x/y/z 軸標題
end

% ============================================================================
function h = draw_plabels(ax, plab, pabs, tip, dh)
    h = gobjects(6,1);
    for k=1:6
        if ~isempty(pabs{k}), t=pabs{k}; else, t=tip(:,k)+340*dh(:,k); end
        h(k) = text(ax,t(1),t(2),t(3),plab{k},'FontSize',22,'FontWeight','bold','Color',[.55 0 .18],'HorizontalAlignment','center');
    end
end

% ============================================================================
function sel = voxpick(mask, vsz, x, y, zwp, bmag, BOX)
    ed = linspace(-BOX,BOX,round(2*BOX/vsz)+1);
    ix=discretize(x,ed); iy=discretize(y,ed); iz=discretize(zwp,ed);
    good = mask & ~isnan(ix)&~isnan(iy)&~isnan(iz);
    nc=numel(ed)-1;  cid=(ix-1)*nc*nc+(iy-1)*nc+iz;
    vi=find(good);  if isempty(vi), sel=[]; return; end
    [~,~,g]=unique(cid(vi));  sel=zeros(max(g),1);
    for k=1:max(g), cn=vi(g==k); [~,b]=max(bmag(cn)); sel(k)=cn(b); end
end

% ============================================================================
function draw_poles_transparent(ax, c)
    ang=c.pole_angles; islow=logical(c.pole_is_lower);
    Rn=c.R_norm*1e3; r0=c.POLE_TIP_R*1e3; L=4;   % mm
    beta=atan2(c.POLE_R,c.POLE_CONE_LEN); inc=c.upper_incline; psi0=atan2(c.R_norm_z,c.R_norm_xy);
    col=[0.72 0.74 0.78];  alp=0.14;
    for k=1:6
        th=ang(k)*pi/180;
        if islow(k), psi=-psi0; axk=[cos(th);sin(th);0];
        else, psi=+psi0; axk=[cos(inc)*cos(th);cos(inc)*sin(th);sin(inc)]; end
        tip=Rn*[cos(psi)*cos(th);cos(psi)*sin(th);sin(psi)];
        if islow(k), v=[0;0;-1]; u=cross(axk,v);
        else, u=cross(axk,[0;0;1]); u=u/norm(u); v=cross(axk,u); end
        draw_pole(ax,tip,axk,u,v,r0,beta,L,islow(k),col,alp);
    end
end

function draw_pole(ax, tip, axk, u, v, rf, beta, L, half, col, alp)
    ts=rf*(1+sin(beta)); rt=rf*cos(beta);
    psi=linspace(0,pi/2+beta,16); nt=34; t=linspace(ts,L,nt);
    axoff=[rf*(1-cos(psi)), t(2:end)]; radm=[rf*sin(psi), rt+(t(2:end)-ts)*tan(beta)];
    if half, phi=linspace(0,pi,44); else, phi=linspace(0,2*pi,64); end
    N=numel(axoff); X=zeros(N,numel(phi)); Y=X; Z=X;
    for j=1:numel(phi)
        rad=u*cos(phi(j))+v*sin(phi(j));  P=tip+axk*axoff+rad*radm;
        X(:,j)=P(1,:).'; Y(:,j)=P(2,:).'; Z(:,j)=P(3,:).';
    end
    surf(ax,X,Y,Z,'FaceColor',col,'FaceAlpha',alp,'EdgeColor','none');
    rim=tip+axk*axoff(end)+(u*cos(phi)+v*sin(phi))*radm(end); cen=tip+axk*axoff(end);
    patch(ax,[cen(1) rim(1,:)],[cen(2) rim(2,:)],[cen(3) rim(3,:)],col,'EdgeColor','none','FaceAlpha',alp);
    if half
        e0=tip+axk*axoff+u*radm; ep=tip+axk*axoff-u*radm; poly=[e0,fliplr(ep)];
        patch(ax,poly(1,:),poly(2,:),poly(3,:),col,'EdgeColor','none','FaceAlpha',alp);
    end
end
