function plot_singlepole_sideview(shapes)
% PLOT_SINGLEPOLE_SIDEVIEW  單極 y=0 (xz) 側視磁路箭頭圖 + bias 電荷位置（no_fix / +x 框）。
%   shapes (optional): cell of shapes, default {'filled','halfcut','tipcut'}.
%   +x mesh：極軸 = +x、WP 在原點、尖端頂點 (0.5,0,0)mm。view x∈[-2,7]、z∈[-3,4] mm、y=0。
%   all-source：coil 已翻向 → raw FEM 即 source（SGN=+1，不再翻 B）。
%   真實節點 grid-sample、均勻箭頭、色=|B|（不內插）。
%   標 bias 電荷位置 = (x_c, z_c) = (ℓ̂, ℓ̂·e_z)（讀 singlepole_bias_fit.mat，粉色圓點；off-axis）
%   + WP 原點 + 極錐輪廓（per shape）。tipcut 不畫 zoom inset。
%   風格① 粗體框圖。出 no_fix_dir/figures/singlepole_sideview_{filled,halfcut,tipcut}.png。

    if nargin < 1 || isempty(shapes), shapes = {'filled','halfcut','tipcut'}; end
    here   = fileparts(mfilename('fullpath'));
    nfdir  = fileparts(fileparts(here));                                                          % no_fix_dir
    FIX_MFUN = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\Calibration_using_FEM_modeling\fix_dir\code\main_function';  % load_singlepole
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');  % mt_constants
    addpath(FIX_MFUN);                                                                            % load_singlepole (reuse)
    figdir = fullfile(nfdir,'figures');  if ~exist(figdir,'dir'); mkdir(figdir); end
    DATA   = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data\coil1\singlepole';

    cnst = mt_constants();
    cnst.SPH_OFST = 0;                                  % +x mesh: WP 在原點（load_singlepole 不 z-shift）
    APEX = cnst.R_norm*1e3;                             % 尖端頂點 x [mm] = 0.5
    POLE_R = cnst.POLE_R*1e3;  POLE_LEN = cnst.POLE_CONE_LEN*1e3;
    slp  = POLE_R/POLE_LEN;                             % 錐半角斜率 = 3/15 = 0.2

    S = load(fullfile(nfdir,'data','singlepole_bias_fit.mat'));   % SHAPES, ell_um, ey, ez（bias fit）

    xlim_v = [-2, 7];   zlim_v = [-3, 4];               % view [mm]
    XV = xlim_v(2);
    nx = 40;  nz = 32;                                  % grid-sample 格數（cell ~0.22mm）

    %% ==== Pass 1：載入 3 形狀 + grid-sample，收集共用色階 CAP ====
    %%   同類比較圖共用色階（figure-style 規則）：3 形狀用同一 colorbar，跨形狀可比、
    %%   最大值數字不誤導（halfcut 近表面峰值真的較高，非 auto-scale 假象）。
    xe=linspace(xlim_v(1),xlim_v(2),nx+1);  ze=linspace(zlim_v(1),zlim_v(2),nz+1);
    D = struct('shp',{},'Xg',{},'Zg',{},'Bx_g',{},'Bz_g',{},'bsum_mT',{},'xmm',{},'zmm',{},'ymm',{},'Bx',{},'Bz',{});
    for s = 1:numel(shapes)
        shp = shapes{s};
        [Pw, BmT] = load_singlepole(fullfile(DATA, shp), cnst, 'all');   % WP 框 [m]、raw [mT]
        xmm = Pw(:,1)*1e3;  ymm = Pw(:,2)*1e3;  zmm = Pw(:,3)*1e3;
        Bx = BmT(:,1)/1e3;  By = BmT(:,2)/1e3;  Bz = BmT(:,3)/1e3;       % Tesla, SGN=+1（不翻，raw 已 source）
        inwin = xmm>=xlim_v(1)&xmm<=xlim_v(2)&zmm>=zlim_v(1)&zmm<=zlim_v(2)&abs(ymm)<0.8;
        xs=xmm(inwin); ys=ymm(inwin); zs=zmm(inwin);  bxs=Bx(inwin); bys=By(inwin); bzs=Bz(inwin);
        ix=discretize(xs,xe); iz=discretize(zs,ze);  ok=~isnan(ix)&~isnan(iz);   % grid-sample 每格挑最近 y=0
        cellid=(iz(ok)-1)*nx+ix(ok);
        xv=xs(ok); zv=zs(ok); bxv=bxs(ok); byv=bys(ok); bzv=bzs(ok); ayv=abs(ys(ok));
        gc=findgroups(cellid);  gi=(1:numel(cellid))';
        pick=accumarray(gc,gi,[],@(r) r(find(ayv(r)==min(ayv(r)),1)));
        D(s).shp=shp;  D(s).Xg=xv(pick); D(s).Zg=zv(pick);
        D(s).Bx_g=bxv(pick); D(s).Bz_g=bzv(pick);
        D(s).bsum_mT=sqrt(bxv(pick).^2+byv(pick).^2+bzv(pick).^2)*1e3;
        D(s).xmm=xmm; D(s).zmm=zmm; D(s).ymm=ymm; D(s).Bx=Bx; D(s).Bz=Bz;      % tipcut inset 用（完整節點）
        fprintf('== %s ==  arrows: %d | |B| max=%.4g mT\n', shp, numel(pick), max(D(s).bsum_mT));
    end
    CAP = max(arrayfun(@(d) prctile(d.bsum_mT, 92), D));    % 3 形狀共用色階上限
    fprintf('common CAP = %.4g mT (shared over %d shapes)\n', CAP, numel(shapes));

    %% ==== Pass 2：逐形狀用共用 CAP 渲染 ====
    arrow_len = 0.9*(xlim_v(2)-xlim_v(1))/nx;
    n_bins=28; edges_b=linspace(0,CAP,n_bins+1); cmap_n=turbo(n_bins);
    for s = 1:numel(shapes)
        shp = D(s).shp;
        Xg=D(s).Xg; Zg=D(s).Zg; Bx_g=D(s).Bx_g; Bz_g=D(s).Bz_g;  bsum_mT=D(s).bsum_mT;
        xmm=D(s).xmm; zmm=D(s).zmm; ymm=D(s).ymm; Bx=D(s).Bx; Bz=D(s).Bz;     % for tipcut inset
        Bip = hypot(Bx_g, Bz_g);  Bip(Bip<eps)=eps;
        bx_q = Bx_g./Bip*arrow_len;  bz_q = Bz_g./Bip*arrow_len;
        mag_c = min(bsum_mT, CAP);

        fig = figure('Position',[60 60 1000 780],'Color','w');
        ax=axes(fig); hold(ax,'on');
        for k=1:n_bins
            in = mag_c>=edges_b(k) & mag_c<edges_b(k+1);
            if k==n_bins, in = in | (mag_c>=edges_b(end)); end
            if any(in)
                quiver(ax, Xg(in),Zg(in), bx_q(in),bz_q(in), 0, 'Color',cmap_n(k,:), 'LineWidth',1.0, 'MaxHeadSize',0.6);
            end
        end

        %% 極錐輪廓（+x 框，per shape）：apex (APEX,0)，斜率 ±slp，延伸到 x=XV
        switch shp
          case 'halfcut'     % 切掉 z>0 半：平頂 z=0 + 下斜邊
            plot(ax,[APEX,XV],[0,0],                 'k-','LineWidth',2);           % 平頂
            plot(ax,[APEX,XV],[0,-(XV-APEX)*slp],    'k-','LineWidth',2);           % 下斜邊
          case 'tipcut'      % 尖端 x=APEX+0.01 切平：垂直平頭面(極小) + 上下斜邊
            Xc=APEX+0.01;  dz=(Xc-APEX)*slp;
            plot(ax,[Xc,Xc],[-dz,dz],               'k-','LineWidth',2);           % 平頭面
            plot(ax,[Xc,XV],[dz,(XV-APEX)*slp],     'k-','LineWidth',2);           % 上斜邊
            plot(ax,[Xc,XV],[-dz,-(XV-APEX)*slp],   'k-','LineWidth',2);           % 下斜邊
          otherwise          % filled：全錐
            plot(ax,[APEX,XV],[0,(XV-APEX)*slp],    'k-','LineWidth',2);           % 上斜邊
            plot(ax,[APEX,XV],[0,-(XV-APEX)*slp],   'k-','LineWidth',2);           % 下斜邊
        end
        plot(ax,0,0,'k+','MarkerSize',15,'LineWidth',2.4);

        %% bias 電荷位置 = (ℓ̂, ℓ̂·e_z) [mm]（粉色圓點、無標籤；off-axis）
        si = find(strcmp(S.SHAPES, shp), 1);
        cx = S.ell_um(si)/1000;  cz = S.ell_um(si)*S.ez(si)/1000;
        plot(ax,cx,cz,'o','MarkerSize',7,'MarkerFaceColor',[1 0.45 0.75],'MarkerEdgeColor',[0.80 0.20 0.50],'LineWidth',1.0);

        %% axes（風格①）
        colormap(ax,turbo);  clim(ax,[0 CAP]);
        axis(ax,'equal');  box(ax,'on');  grid(ax,'off');  xlim(ax,xlim_v);  ylim(ax,zlim_v);
        set(ax,'FontSize',16,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]);
        xt=get(ax,'XTick'); set(ax,'XTick',xt(1:2:end));
        yt=get(ax,'YTick'); set(ax,'YTick',yt(1:2:end));
        xlabel(ax,'x (mm)','FontWeight','bold');  ylabel(ax,'z (mm)','FontWeight','bold');
        cb=colorbar(ax);  cb.FontSize=16; cb.FontWeight='bold';
        cb.Ticks = cb.Ticks(1:2:end);
        cb.Label.String='|B| (mT)'; cb.Label.FontWeight='bold'; cb.Label.FontSize=16;
        ax.Toolbar.Visible='off';

        %% ===== tipcut 尖端 zoom：主圖左上角 inset（純粹加圖；不另出獨立 zoom 圖）=====
        if strcmp(shp,'tipcut')
            Rt = cnst.POLE_TIP_R*1e3;
            xin=[0.45 0.62];  zin=[-0.08 0.08];
            rectangle(ax,'Position',[xin(1) zin(1) diff(xin) diff(zin)],'EdgeColor','k','LineWidth',1.2);   % 主圖標 zoom 框
            axi = axes(fig,'Position',[0.150 0.605 0.30 0.32],'Color','w');
            draw_tip_zoom(axi, xmm,zmm,ymm,Bx,Bz, APEX,slp,Rt);
            set(axi,'XTick',[],'YTick',[],'LineWidth',1.5);        % 純畫圖：無刻度/軸標/title，只留外框
        end

        out = fullfile(figdir, sprintf('singlepole_sideview_%s.png', shp));
        exportgraphics(fig, out, 'Resolution',600);
        fprintf('  saved %s\n', out);
        close(fig);
    end
end

function draw_tip_zoom(axt, xmm,zmm,ymm,Bx,Bz, APEX,slp,Rt)
% 尖端 zoom 內容：細格 grid-sample 灰箭頭 + 完整尖端(dashed 灰) + tipcut(solid 黑) 輪廓。
%   圓角球心 Cx=APEX+Rt、半徑 Rt；平切 XCUT=0.510mm；圓角↔錐相切點 Ptan。畫在 axt（inset）。
    Cx=APEX+Rt;  bt=atan(slp);
    XCUT=0.510;  zf=sqrt(Rt^2-(Cx-XCUT)^2);
    th_tan=atan2(cos(bt),-sin(bt));  Ptan=[Cx-Rt*sin(bt), Rt*cos(bt)];
    xin=[0.45 0.62];  zin=[-0.08 0.08];
    hold(axt,'on');
    % 場箭頭（細格 grid-sample 真實節點）
    inw=xmm>=xin(1)&xmm<=xin(2)&zmm>=zin(1)&zmm<=zin(2)&abs(ymm)<0.15;
    xw=xmm(inw); zw=zmm(inw); bxw=Bx(inw); bzw=Bz(inw); aw=abs(ymm(inw));
    nxi=12; xei=linspace(xin(1),xin(2),nxi+1); zei=linspace(zin(1),zin(2),nxi+1);
    ixi=discretize(xw,xei); izi=discretize(zw,zei); oki=~isnan(ixi)&~isnan(izi);
    cid=(izi(oki)-1)*nxi+ixi(oki); xq=xw(oki); zq=zw(oki); bxq=bxw(oki); bzq=bzw(oki); aq=aw(oki);
    if ~isempty(cid)
        g=findgroups(cid); ii=(1:numel(cid))';
        pk=accumarray(g,ii,[],@(r) r(find(aq(r)==min(aq(r)),1)));
        al=0.9*diff(xin)/nxi; bp=hypot(bxq(pk),bzq(pk)); bp(bp<eps)=eps;
        quiver(axt,xq(pk),zq(pk),bxq(pk)./bp*al,bzq(pk)./bp*al,0,'Color',[.35 .35 .35],'LineWidth',0.8,'MaxHeadSize',0.5);
    end
    % 完整尖端(dashed 灰)：圓角 apex->切點 + 錐
    tf=linspace(pi,th_tan,40);
    plot(axt,Cx+Rt*cos(tf),Rt*sin(tf),'--','Color',[.55 .55 .55],'LineWidth',1.6);
    plot(axt,Cx+Rt*cos(tf),-Rt*sin(tf),'--','Color',[.55 .55 .55],'LineWidth',1.6);
    plot(axt,[Ptan(1),xin(2)],[Ptan(2),Ptan(2)+(xin(2)-Ptan(1))*slp],'--','Color',[.55 .55 .55],'LineWidth',1.6);
    plot(axt,[Ptan(1),xin(2)],[-Ptan(2),-(Ptan(2)+(xin(2)-Ptan(1))*slp)],'--','Color',[.55 .55 .55],'LineWidth',1.6);
    % tipcut(solid 黑)：平切面 + 圓角(切點->平切) + 錐
    plot(axt,[XCUT,XCUT],[-zf,zf],'k-','LineWidth',2);
    th_f=atan2(zf,XCUT-Cx); tc=linspace(th_f,th_tan,25);
    plot(axt,Cx+Rt*cos(tc),Rt*sin(tc),'k-','LineWidth',2);
    plot(axt,Cx+Rt*cos(tc),-Rt*sin(tc),'k-','LineWidth',2);
    plot(axt,[Ptan(1),xin(2)],[Ptan(2),Ptan(2)+(xin(2)-Ptan(1))*slp],'k-','LineWidth',2);
    plot(axt,[Ptan(1),xin(2)],[-Ptan(2),-(Ptan(2)+(xin(2)-Ptan(1))*slp)],'k-','LineWidth',2);
    axis(axt,'equal'); box(axt,'on'); xlim(axt,xin); ylim(axt,zin);
end
