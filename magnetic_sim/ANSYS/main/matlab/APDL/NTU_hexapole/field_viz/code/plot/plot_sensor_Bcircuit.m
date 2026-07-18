function plot_sensor_Bcircuit(PREVIEW, ZWIN, XWIN, SHOWSTEEL, MODEL)
% PLOT_SENSOR_BCIRCUIT  NTU 單極 sensor 附近磁路箭頭圖（x-z @ y=-0.99，含 S1+S2）
%   比照長飛 plot_P2sensor_Braw_P1exc.m 風格：真實 FEM 節點箭頭（不內插）、
%   per-arrow |B| log 色階(turbo)、sensor 綠 disc + 紅 n+、鋼件(薄板)截面輪廓、粗體框圖①。
%   資料 = NTU 1A×50匝 solve 的全節點場（data/singlepole/singlepole_{coord,bfield}_all.dat）。
%   pole frame（NTU WP=原點、無 SPH_OFST 位移）。raw 物理場（單極單線圈、不翻號）。
%
%   用法： plot_sensor_Bcircuit                          % 存 figures/sensor_Bcircuit_xz.png（z:-2~2）
%          plot_sensor_Bcircuit(false,[-2 11.25],[10 27.24]) % 高 z 版 → figures/sensor_Bcircuit_xz_tall.png
%          plot_sensor_Bcircuit(true)                    % PREVIEW 暫存
    if nargin<1 || isempty(PREVIEW), PREVIEW=false;    end
    if nargin<2 || isempty(ZWIN),    ZWIN=[-2 2];      end   % [MODIFIED] z 視野可調（預設同現況）
    if nargin<3 || isempty(XWIN),    XWIN=[10 24.74];  end   % [MODIFIED] x 視野可調（預設同現況）
    if nargin<4 || isempty(SHOWSTEEL),SHOWSTEEL=false;  end   % [MODIFIED] 疊導柱+yoke 鋼件輪廓（預設關）
    if nargin<5 || isempty(MODEL),   MODEL='upper_assembly'; end   % [MODIFIED] 資料/幾何模型：'upper_assembly'|'singlepole'
    isSP = strcmp(MODEL,'singlepole');                       % 單極（1 極板 + 1 短導柱、無 yoke）

    %% ---- config ----
    DPI   = 200;
    SLAB  = 0.30;                 % y 薄片半厚 [mm]（繞 sensor 的 y=-0.99）
    SLAB3 = 12.0;                 % [ADDED] gap-fill 內插 3D Delaunay 取樣 y 半寬 [mm]（放大 hull 填遠場粗網格白塊）
    YC    = -0.99;                % 切面 y [mm]
    WIN   = [XWIN ZWIN];          % [xlo xhi zlo zhi] mm（upper_assembly 視野；x/z 由 XWIN/ZWIN 參數化）
    ARROW = 0.09;                 % 箭頭長 [mm]
    CELL  = 0.09;                 % 格點抽稀邊長 [mm]（每格留 |B| 最大真實節點；非內插）
    NLEN  = 0.38;                 % n+ 箭頭長 [mm]
    nb    = 24;                   % 色階 bin 數

    %% ---- sensor 幾何（pole frame, mm）----
    S1 = [13.15, -0.99, -0.41];  n1 = [0 0 -1];   % n+ -> -z
    S2 = [13.15, -0.99,  0.66];  n2 = [0 0 +1];   % n+ -> +z

    %% ---- 薄板截面（此切面 y=-0.99：z∈[0,0.25]，x 橫跨視窗）----
    PLATE_Z = [0 0.25];

    %% ---- paths + 載資料（重用 backup import_ansys_data，處理 Fortran 定寬負號）----
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
    rr = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\NTU_hexapole\data';
    d  = import_ansys_data(fullfile(rr,MODEL,'coil1'),'all',MODEL);   % SI [m],[T]（folder+prefix 都用 MODEL）
    xw = d.x*1e3;  yw = d.y*1e3;  zw = d.z*1e3;                             % pole frame [mm]
    fprintf('載入 %d 節點 ｜ |B| max = %.4f T\n', numel(d.x), max(sqrt(d.bx.^2+d.by.^2+d.bz.^2)));

    %% ---- 選真實節點（y 薄片 + 視窗盒；不內插）----
    m = abs(yw-YC)<SLAB & xw>=WIN(1) & xw<=WIN(2) & zw>=WIN(3) & zw<=WIN(4);
    fprintf('slab |y-(%.2f)|<%.2fmm 視窗內真實節點 = %d\n', YC, SLAB, nnz(m));
    Xx=xw(m); Zz=zw(m); Bx=d.bx(m); By=d.by(m); Bz=d.bz(m);
    Bmag=sqrt(Bx.^2+By.^2+Bz.^2)*1e3;                                       % |B| [mT]（unit-reference）
    Bip=hypot(Bx,Bz); Bip(Bip<eps)=eps;
    ux=Bx./Bip*ARROW; uz=Bz./Bip*ARROW;                                    % 單位方向 × 箭頭長

    % 格點抽稀：每格留 |B| 最大之真實節點（峰值不掉、仍是節點原值、非內插）
    if CELL>0
        [~,~,gid]=unique(round([Xx Zz]/CELL),'rows');
        iu=accumarray(gid,(1:numel(Bmag))',[],@(r) r(find(Bmag(r)==max(Bmag(r)),1)));
        Xx=Xx(iu); Zz=Zz(iu); ux=ux(iu); uz=uz(iu); Bmag=Bmag(iu);
    end

    %% ---- gap-fill 重心內插（填真實節點稀疏的空洞：板內 + 板附近粗空氣）----
    %  薄板板內無真實節點、板↔sensor 球間空氣粗（節點稀）→ 對「附近無真實節點」的格點內插補箭頭。
    %  單一局部 3D Delaunay（全節點）：板內格點落鋼件 tet（鋼場）、空氣空洞格點落跨界 tet（過渡場）。
    %  B(p)=Σλ_iB_i（=FEM 線性 shape fn）。板外真實節點密處仍用真實節點、不重複。
    % [MODIFIED] sel3 box 放大（x/z 邊界 +MG，y 半寬 SLAB3）→ 3D hull 涵蓋整個視野；
    %   hull 外的極角落再以「最近節點」fallback 補上（仍屬 gap-fill、不留白）。
    MG=8;                                                             % sel3 x/z 邊界外擴 [mm]（納入視野外圍節點撐大 hull）
    sel3=(xw>WIN(1)-MG&xw<WIN(2)+MG)&(abs(yw-YC)<SLAB3)&(zw>WIN(3)-MG&zw<WIN(4)+MG);
    Pall=[xw(sel3) yw(sel3) zw(sel3)]; Ball=[d.bx(sel3) d.by(sel3) d.bz(sel3)];
    Xi=[]; Zi=[]; uxi=[]; uzi=[]; Bmi=[];
    if size(Pall,1)>=20
        TRg=delaunayTriangulation(Pall);
        [gx,gz]=meshgrid(WIN(1)+0.03:0.12:WIN(2)-0.03, WIN(3)+0.04:0.075:WIN(4)-0.04);
        gp=[gx(:) YC+0*gx(:) gz(:)];
        dmin=min(sqrt((gp(:,1)-Xx.').^2+(gp(:,3)-Zz.').^2),[],2);      % 到最近真實顯示節點距離
        gp=gp(dmin>0.10,:);                                           % 只填空洞（0.10mm 內無真實節點）
        if ~isempty(gp)
            ti=pointLocation(TRg,gp); ok=~isnan(ti);
            Bg=zeros(size(gp,1),3);
            if any(ok)                                               % hull 內：重心內插（=FEM 線性 shape fn）
                bc=cartesianToBarycentric(TRg,ti(ok),gp(ok,:));
                cn=TRg.ConnectivityList(ti(ok),:);
                for c=1:4, Bg(ok,:)=Bg(ok,:)+bc(:,c).*Ball(cn(:,c),:); end
            end
            if any(~ok)                                             % [ADDED] hull 外極角落：最近節點 fallback
                vi=nearestNeighbor(TRg, gp(~ok,:));
                Bg(~ok,:)=Ball(vi,:);
            end
            Xi=gp(:,1); Zi=gp(:,3); Bmi=sqrt(sum(Bg.^2,2))*1e3;      % mT
            Bipi=hypot(Bg(:,1),Bg(:,3)); Bipi(Bipi<eps)=eps;
            uxi=Bg(:,1)./Bipi*ARROW; uzi=Bg(:,3)./Bipi*ARROW;
        end
    end
    fprintf('gap-fill interp arrows = %d\n', numel(Xi));

    %% ---- 圖 ----
    % [MODIFIED] 圖框尺寸自適應 WIN 長寬比（axis equal 下不留大量左右白邊）；
    %   校準：預設(x14.74,z4)→1760x600（重現原圖）；高 z(z13.25)→1760x1479。
    xspan=WIN(2)-WIN(1); zspan=WIN(4)-WIN(3);
    PXMM=95; figW=round(xspan*PXMM)+360; figH=round(zspan*PXMM)+220;
    fig=figure('Position',[40 40 figW figH],'Color','w','Visible','off');
    ax=axes(fig); hold(ax,'on');

    % 薄板截面：**只黑框輪廓、不填色**（讓磁極內磁路箭頭看得見）
    %   [MODIFIED] singlepole 極板 y=YC 右緣止於 x=25.40（far apex r5 弧 @ (20.5,0)）；upper 維持全窗。
    PXR = WIN(2);  if isSP, PXR = min(WIN(2), 25.40); end
    plot([WIN(1) PXR PXR WIN(1) WIN(1)], ...
         [PLATE_Z(1) PLATE_Z(1) PLATE_Z(2) PLATE_Z(2) PLATE_Z(1)],'k-','LineWidth',2);

    % [ADDED] 鋼件輪廓（y=YC 切面 2D silhouette；來自 mesh deck 幾何）
    if SHOWSTEEL
        hw=sqrt(2.5^2-YC^2); pl=20.5-hw; pr=20.5+hw;               % 導柱 @ (20.5,0) r2.5 → 兩壁 x=18.20/22.80
        if isSP
            % 單極：只有導柱（下 r2.5 z0.25-11.25 + 上 r2.49 z11.25-17.25，step 0.01mm 忽略）；**無 yoke**
            plot([pl pr pr pl pl],[0.25 0.25 17.25 17.25 0.25],'k-','LineWidth',2);   % 柱身 z0.25→17.25（ylim clip 頂）
            text(20.5,8,'post','FontSize',12,'FontWeight','bold','Color',[.2 .2 .2],'HorizontalAlignment','center');
        else
            % upper_assembly：導柱 + yoke 融合 silhouette（yoke z[11.25,16.75]、cutout 內壁 x=15.5、右緣 clip）
            xr=WIN(2);
            sx=[xr    pr    pr   pl   pl    15.5  15.5  xr  ];
            sz=[11.25 11.25 0.25 0.25 11.25 11.25 16.75 16.75];
            plot(sx,sz,'k-','LineWidth',2);
            text(20.5, 7,'post','FontSize',12,'FontWeight','bold','Color',[.2 .2 .2],'HorizontalAlignment','center');
            text(24.5,14,'yoke','FontSize',12,'FontWeight','bold','Color',[.2 .2 .2],'HorizontalAlignment','center');
        end
    end

    % per-arrow log 色階（|B| 跨數量級；含板內內插點以統一色階）
    Ca=log10(max([Bmag;Bmi],1e-12));                                       % 合併真實+內插 定色階範圍
    Cv=log10(max(Bmag,1e-12)); cmap=turbo(nb); ed=linspace(min(Ca),max(Ca),nb+1);
    for k=1:nb                                                             % 真實節點箭頭
        in=Cv>=ed(k)&Cv<ed(k+1); if k==nb, in=in|(Cv>=ed(end)); end
        if any(in)
            quiver(Xx(in),Zz(in),ux(in),uz(in),0,'Color',cmap(k,:),'LineWidth',1.3,'MaxHeadSize',2.0);
        end
    end
    if ~isempty(Xi)                                                        % 板內內插箭頭（同色階、同風格）
        Cvi=log10(max(Bmi,1e-12));
        for k=1:nb
            in=Cvi>=ed(k)&Cvi<ed(k+1); if k==nb, in=in|(Cvi>=ed(end)); end
            if any(in)
                quiver(Xi(in),Zi(in),uxi(in),uzi(in),0,'Color',cmap(k,:),'LineWidth',1.3,'MaxHeadSize',2.0);
            end
        end
    end
    colormap(ax,turbo); clim(ax,[min(Ca) max(Ca)]);
    cb=colorbar(ax); ylabel(cb,'|B| (3D) (mT)');
    lo=floor(min(Ca)); hi=floor(max(Ca)); tk=lo:hi;
    if ~isempty(tk)&&(max(Ca)-tk(end))<0.15, tk(end)=[]; end
    labs=arrayfun(@(t)sprintf('10^{%d}',t),tk,'UniformOutput',false);
    tk=[tk, max(Ca)]; labs=[labs,{sprintf('%.3g',10^max(Ca))}];
    cb.Ticks=tk; cb.TickLabels=labs;

    % sensor：綠 disc + disc 平面短線(⊥ n+) + 紅 n+ 箭頭 + 標籤
    for S=[struct('p',S1,'n',n1,'lab','S1'), struct('p',S2,'n',n2,'lab','S2')]
        sx=S.p(1); sz=S.p(3); nx=S.n(1); nz=S.n(3);
        tdir=[-nz; nx];                                                    % disc 平面方向(⊥ n+)
        dd=[sx;sz]-0.16*tdir; de=[sx;sz]+0.16*tdir;
        plot([dd(1) de(1)],[dd(2) de(2)],'-','Color',[.1 .5 .15],'LineWidth',3);
        plot(sx,sz,'o','MarkerSize',9,'MarkerFaceColor',[.1 .7 .25],'MarkerEdgeColor','k','LineWidth',1);
        nt=[sx;sz]+NLEN*[nx;nz];
        plot([sx nt(1)],[sz nt(2)],'-','Color','w','LineWidth',5.5);
        plot([sx nt(1)],[sz nt(2)],'-','Color',[.92 .15 .15],'LineWidth',3.6);
        plot(nt(1),nt(2),'^','MarkerSize',12,'MarkerFaceColor',[.92 .15 .15],'MarkerEdgeColor','w','LineWidth',1.0);
        text(sx+0.10, sz, S.lab,'FontSize',14,'FontWeight','bold','Color','k');
        text(nt(1)+0.05, nt(2),'n_+','Color',[.92 .15 .15],'FontSize',15,'FontWeight','bold','Interpreter','tex');
    end
    text(mean(WIN(1:2)), PLATE_Z(2)+0.02, 'pole plate','FontSize',12,'FontWeight','bold', ...
         'Color',[.3 .3 .3],'HorizontalAlignment','center','VerticalAlignment','bottom');

    hold(ax,'off'); axis(ax,'equal'); xlim(ax,WIN(1:2)); ylim(ax,WIN(3:4));
    xlabel(ax,'x (mm)'); ylabel(ax,'z (mm)');
    if isSP, ttl='NTU single pole sensor'; else, ttl='NTU upper assembly sensor'; end   % [MODIFIED] 標題依 model
    title(ax,[ttl ' ｜ 1A x 50T ｜ raw'],'Interpreter','none');
    ax.Toolbar.Visible='off';

    %% ---- [STYLE] 粗體框圖①（figure-style.md）----
    set(ax,'FontSize',16,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]);
    box(ax,'on'); grid(ax,'off');
    xt=get(ax,'XTick'); set(ax,'XTick',xt(1:2:end));
    yt=get(ax,'YTick'); set(ax,'YTick',yt(1:2:end));
    set([get(ax,'XLabel') get(ax,'YLabel')],'FontSize',17,'FontWeight','bold');
    set(get(ax,'Title'),'FontSize',15,'FontWeight','bold');
    cb.FontSize=15; cb.FontWeight='bold'; cb.LineWidth=1.5;
    set(get(cb,'YLabel'),'FontSize',16,'FontWeight','bold');

    %% ---- 輸出 ----
    figdir='G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\NTU_hexapole\field_viz\figures';
    isTall = ~isequal(WIN(3:4),[-2 2]);                        % [MODIFIED] 非預設 z 視野→另存 _tall
    base   = 'sensor_Bcircuit'; if isSP, base=[base '_singlepole']; end   % [MODIFIED] singlepole 另存
    base   = [base '_xz'];      if isTall, base=[base '_tall']; end
    if PREVIEW, out=fullfile(tempdir,[base '_preview.png']);
    else,       out=fullfile(figdir,[base '.png']); end
    exportgraphics(fig,out,'Resolution',DPI);
    fprintf('saved: %s\n', out);
    close(fig);
end
