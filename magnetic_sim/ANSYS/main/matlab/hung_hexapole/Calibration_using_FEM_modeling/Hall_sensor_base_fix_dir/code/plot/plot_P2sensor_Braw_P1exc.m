function plot_P2sensor_Braw_P1exc(MODE, PREVIEW, VARIANT, FLIP, WIN, FILLEMPTY)
%   WIN (optional) = [xlo xhi zlo zhi] mm (WP frame)：覆寫節點選取盒 + 軸範圍
%                    （給自訂視窗用，如 [-4 4 0 7]）；空=用 MODE 預設方框。
%   FILLEMPTY (optional, 預設 false)：true=把「沒有真實節點的空格子」用 scatteredInterpolant
%                    ('linear','none') 補一支箭頭（只填空、real node 不動；凸包外不外插）。
% PLOT_P2SENSOR_BRAW_P1EXC  P2-sensor 局部磁路箭頭圖
%   coil1 (P1 excited)、all-source 翻 P1、真實 FEM 節點（不內插）。
%   資料源由 VARIANT 選：
%     'sensor_spheres'（預設）：sensor 加密網格。
%     'standard'              ：Long2016 verbatim baseline 粗網格。
%   三種呈現模式（MODE 參數）：
%     'circuit'：±2mm 廣角側視（同原始參考圖風格）、y 薄片 + 0.16mm 格點抽樣。
%     'pole'   ：±2.7mm 全磁路 WP→P2→sensor、對數色階。
%     'zoom'   ：±0.55mm 放大、0.045mm 格點抽樣 → 看得到 sensor 處場被細緻解析（箭頭變密）。
%     'nodes'  ：緊貼 sensor，直接畫 extract_Vmat「取樣圓柱」(半徑0.15mm、沿 n+ 高0.1mm)內的
%                真實節點，一節點一箭頭 → 箭頭數 = 實際被平均的節點數（加密版 P2 = 168）。
%                注意 baseline 粗網格圓柱內 0 節點，'nodes' 模式只在 'sensor_spheres' 適用。
%   箭頭=單位方向（看角度），顏色=|B|(3D, Tesla)。
%
%   用法： plot_P2sensor_Braw_P1exc('circuit', false, 'standard') % 存 figures/..._circuit_standard.png
%          plot_P2sensor_Braw_P1exc('pole',    false, 'standard') % 存 figures/..._pole_standard.png
%          plot_P2sensor_Braw_P1exc('nodes',   false)             % 加密版 figures/..._nodes.png
%   PREVIEW=true 時輸出到暫存 png（定案前不落 figures/）。

    if nargin < 1 || isempty(MODE),    MODE    = 'nodes';          end
    if nargin < 2 || isempty(PREVIEW), PREVIEW = true;             end
    if nargin < 3 || isempty(VARIANT), VARIANT = 'sensor_spheres'; end
    if nargin < 4 || isempty(FLIP),    FLIP    = 'flip';            end  % 'flip'=all-source（翻 P1）/ 'raw'=原始物理場
    if nargin < 5,                     WIN     = [];                end  % []=用 MODE 預設方框 / [xlo xhi zlo zhi]=自訂視窗
    if nargin < 6 || isempty(FILLEMPTY), FILLEMPTY = false;         end  % true=只內插填空格子
    if ~ismember(VARIANT,{'sensor_spheres','standard'}) && isempty(regexp(VARIANT,'^gap\d+um_mueq$','once'))
        error('VARIANT 必須是 sensor_spheres / standard / gap<NN>um_mueq');
    end
    if ~ismember(FLIP,{'flip','raw'})
        error('FLIP 必須是 flip（all-source）/ raw（原始）');
    end
    if strcmp(MODE,'nodes') && strcmp(VARIANT,'standard')
        error('baseline 粗網格 sensor 圓柱內 0 節點，nodes 模式只適用 sensor_spheres');
    end

    %% ---- config ----
    DPI       = 200;
    SENSOR_R  = 0.15e-3;  % 取樣圓柱半徑 [m]（= extract_Vmat）
    AXIAL_TOL = 0.10e-3;  % 取樣圓柱高 [m]（沿 n+，= extract_Vmat）
    switch MODE
        case 'circuit', SLAB_MM=0.30; HW_MM=2.0;  ARROW_MM=0.18;  CELL=0.16;  NLEN=0.85;
        case 'zoom',    SLAB_MM=0.20; HW_MM=0.55; ARROW_MM=0.055; CELL=0.045; NLEN=0.45;
        case 'pole',    SLAB_MM=0.30; HW_MM=2.7;  ARROW_MM=0.16;  CELL=0.17;  NLEN=0.60;
        case 'why',     SLAB_MM=0.25; HW_MM=2.3;  ARROW_MM=0.13;  CELL=0.13;  NLEN=0.50;
        case 'nodes',   SLAB_MM=0.20; HW_MM=0.25; ARROW_MM=0.035; CELL=0;     NLEN=0.18;
        otherwise, error('MODE 必須是 circuit / zoom / pole / why / nodes');
    end

    %% ---- paths ----
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
    addpath(['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
             'long2016_hexapole_halfcut\Calibration_using_FEM_modeling\Hall_sensor_base_fix_dir\code\function']);
    cnst = mt_constants();
    rr   = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data';

    %% ---- P2 sensor 幾何（i=2）----
    [sensor_pos, sensor_n] = build_sensor_geometry(cnst);   % WP 框 [m]
    sp = sensor_pos(:,2);  ni = sensor_n(:,2);              % P2 sensor 中心/法線 (3×1, [m])
    sx = sp(1)*1e3;  sz = sp(3)*1e3;                        % P2 sensor (x,z) [mm]
    nx = ni(1);  nz = ni(3);  nn = hypot(nx,nz); nx=nx/nn; nz=nz/nn;   % n+ 投到 x-z 後正規化
    fprintf('P2 sensor (x,z) = (%.3f, %.3f) mm ; n+(x,z)=(%.3f,%.3f)\n', sx,sz,nx,nz);

    %% ---- 視窗中心（決定選取盒與 xlim/ylim）----
    %  pole：WP 原點與 sensor 中點（看磁路 WP→P2→sensor）；why：轉換點↔sensor 中點；circuit/zoom：sensor；nodes：節點雲質心
    SOFF_FLIP = 3.0;   % 'why'：匯進↔漏出 sign-flip 約在沿錐面距尖 ~3mm（來自 diag_P2_Bn_map）
    if strcmp(MODE,'pole')
        vx=(0+sx)/2; vz=(0+sz)/2;
    elseif strcmp(MODE,'why')
        incw=cnst.upper_incline; betaw=atan2(3.0,15.0);
        Twh =[cnst.R_norm_xy*cosd(180), cnst.R_norm_z]*1e3;            % P2 tip (x,z) [mm]
        rotw=@(v,a)[cos(a)*v(1)-sin(a)*v(2); sin(a)*v(1)+cos(a)*v(2)];
        flank_dir = rotw([-cos(incw);sin(incw)], betaw);              % 上側 flank 方向(x-z, unit)
        Pflip = Twh(:) + SOFF_FLIP*flank_dir;                         % sign-flip 轉換點
        vx=(Pflip(1)+sx)/2; vz=(Pflip(2)+sz)/2;
    else
        vx=sx; vz=sz;
    end

    %% ---- 載 coil1 (P1 excited)，套 all-source（P1 下極→negate 使尖端射出/source）----
    %  per charge-model-source-convention：每顆激發極一律當 source。負號 ≡ 反向繞線重跑，
    %  仍是真實 FEM 場（節點原值、不內插）。
    d = import_ansys_data(fullfile(rr,'coil1',VARIANT),'all','coil1');
    if strcmp(FLIP,'flip')                                         % all-source：翻 P1（下極）→ 每極當 source
        d.bx = -d.bx;  d.by = -d.by;  d.bz = -d.bz;
    end                                                            % 'raw'：保留原始物理場（sensor 處磁通匯進 P2）
    fprintf('資料源 = coil1/%s ｜ %s ｜ 載入 %d 節點 ｜ |B| max = %.4f T\n', ...
            VARIANT, FLIP, numel(d.x), max(sqrt(d.bx.^2+d.by.^2+d.bz.^2)));   % 核指紋（baseline ~1.0 T）
    xw = d.x*1e3;  yw = d.y*1e3;  zw = (d.z - cnst.SPH_OFST)*1e3;   % WP 框 [mm]

    %% ---- 選節點 ----
    if strcmp(MODE,'nodes')
        % 直接套 extract_Vmat 的取樣圓柱（3D）：rho<=R 且 0<=axial<=H（沿 n+）
        r3    = [d.x, d.y, d.z - cnst.SPH_OFST] - sp.';            % 各節點相對 sensor 中心 [m]
        axial = r3 * ni;                                           % 沿 n+ 軸向 [m]
        rho   = vecnorm(r3 - axial*ni.', 2, 2);                    % 徑向 [m]
        m = (rho <= SENSOR_R) & (axial >= 0) & (axial <= AXIAL_TOL);
        fprintf('取樣圓柱內真實節點數 = %d\n', nnz(m));
    else
        % 廣角/放大/pole/why：y 薄片 + zoom 盒（以 vx,vz 為中心）；WIN 給定時改用自訂矩形
        if ~isempty(WIN)
            m = abs(yw) < SLAB_MM & xw>=WIN(1) & xw<=WIN(2) & zw>=WIN(3) & zw<=WIN(4);
        else
            m = abs(yw) < SLAB_MM & abs(xw - vx) < HW_MM & abs(zw - vz) < HW_MM;
        end
        if strcmp(MODE,'why')      % why：只留空氣節點（看空氣磁路、去掉鋼體內部場的雜訊）
            airsel = filter_iron_nodes(d.x,d.y,d.z,cnst,struct('visualize',false));
            isair = false(numel(d.x),1); isair(airsel)=true;
            m = m & isair;
        end
        fprintf('zoom 盒內 |y|<%.2fmm 真實節點數 = %d\n', SLAB_MM, nnz(m));
    end
    nSel = nnz(m);
    Xx = xw(m); Zz = zw(m);
    Bx = d.bx(m); By = d.by(m); Bz = d.bz(m);
    Bmag = sqrt(Bx.^2 + By.^2 + Bz.^2);              % |B| 3D [T]
    Bip  = hypot(Bx, Bz);  Bip(Bip<eps)=eps;         % in-plane mag
    ux = Bx./Bip * ARROW_MM;  uz = Bz./Bip * ARROW_MM;   % 單位方向 × 箭頭長
    if CELL > 0
        % 格點抽樣降密度（每格留「|B| 最大」的真實節點，峰值不被抽掉、colorbar 達真實峰值；仍是節點原值、非內插）
        [~,~,gid] = unique(round([Xx Zz]/CELL),'rows');
        iu = accumarray(gid,(1:numel(Bmag))',[],@(r) r(find(Bmag(r)==max(Bmag(r)),1)));
        Xx=Xx(iu); Zz=Zz(iu); ux=ux(iu); uz=uz(iu); Bmag=Bmag(iu);
        Bx=Bx(iu); By=By(iu); Bz=Bz(iu);                 % 留分量當內插源（占用格＝唯一、無重複點）
    end

    %% ---- [FILLEMPTY] 只內插「沒有真實節點的空格子」（real node 不動；凸包外不外插）----
    if FILLEMPTY && CELL>0 && numel(Xx)>=3
        if ~isempty(WIN), WINx=WIN(1:2); WINz=WIN(3:4);
        else, WINx=[vx-HW_MM vx+HW_MM]; WINz=[vz-HW_MM vz+HW_MM]; end
        Fbx=scatteredInterpolant(Xx,Zz,Bx,'linear','none');
        Fby=scatteredInterpolant(Xx,Zz,By,'linear','none');
        Fbz=scatteredInterpolant(Xx,Zz,Bz,'linear','none');
        ixr=floor(WINx(1)/CELL):ceil(WINx(2)/CELL);
        izr=floor(WINz(1)/CELL):ceil(WINz(2)/CELL);
        [IX,IZ]=meshgrid(ixr,izr); CX=IX(:)*CELL; CZ=IZ(:)*CELL;
        inw = CX>=WINx(1)&CX<=WINx(2)&CZ>=WINz(1)&CZ<=WINz(2);
        CX=CX(inw); CZ=CZ(inw);
        occ=round([Xx Zz]/CELL);                          % 已占用格
        isocc=ismember(round([CX CZ]/CELL),occ,'rows');
        CX=CX(~isocc); CZ=CZ(~isocc);                     % 只留空格
        bxi=Fbx(CX,CZ); byi=Fby(CX,CZ); bzi=Fbz(CX,CZ);
        ok=~isnan(bxi)&~isnan(bzi)&~isnan(byi);           % 凸包外 NaN → 不補（遠場真空維持空白）
        CX=CX(ok);CZ=CZ(ok);bxi=bxi(ok);byi=byi(ok);bzi=bzi(ok);
        Bmi=sqrt(bxi.^2+byi.^2+bzi.^2); bipi=hypot(bxi,bzi); bipi(bipi<eps)=eps;
        fprintf('FILLEMPTY: 補了 %d 個空格子（內插），真實節點 %d 個保持不變\n', numel(CX), numel(Xx));
        Xx=[Xx;CX]; Zz=[Zz;CZ];
        ux=[ux; bxi./bipi*ARROW_MM]; uz=[uz; bzi./bipi*ARROW_MM];
        Bmag=[Bmag; Bmi];
    end

    %% ---- 圖 ----
    fig = figure('Position',[60 60 1000 950],'Color','w');
    ax = axes(fig); hold(ax,'on');
    % 上色：pole 用對數色階（場強跨多個數量級，看 WP→sensor 漸變）；其餘線性 |B|
    useLog = ismember(MODE,{'pole','why'});
    if useLog, Cv = log10(max(Bmag,1e-12)); else, Cv = Bmag; end
    nb=24; ed=linspace(min(Cv),max(Cv),nb+1); cmap=turbo(nb);
    for k=1:nb
        in = Cv>=ed(k) & Cv<ed(k+1); if k==nb, in=in|(Cv>=ed(end)); end
        if any(in)
            quiver(Xx(in),Zz(in),ux(in),uz(in),0,'Color',cmap(k,:), ...
                   'LineWidth',1.3,'MaxHeadSize',2.0);
        end
    end
    colormap(turbo); clim([min(Cv) max(Cv)]);
    cb=colorbar;
    if useLog
        ylabel(cb,'|B| (3D) (T)');
        lo=floor(min(Cv)); hi=floor(max(Cv)); tk=lo:hi;          % integer-decade ticks within range
        if ~isempty(tk) && (max(Cv)-tk(end))<0.15, tk(end)=[]; end  % drop top decade if it crowds the max label
        labs=arrayfun(@(t)sprintf('10^{%d}',t),tk,'UniformOutput',false);
        cmaxReal = 10^max(Cv);                                   % true displayed |B|max [T]
        tk=[tk, max(Cv)]; labs=[labs, {sprintf('%.3g',cmaxReal)}]; % [ADDED] label the real max at colorbar top
        cb.Ticks=tk; cb.TickLabels=labs;
    else
        ylabel(cb,'|B| (3D) (T)');
    end

    %% ---- P2 cone 局部外框（廣角/放大模式）----
    if ~strcmp(MODE,'nodes')
        inc=cnst.upper_incline; beta=atan2(3.0,15.0);
        T=[cnst.R_norm_xy*cosd(180), cnst.R_norm_z]*1e3;             % P2 tip (x,z) [mm]
        Lsl=sqrt((cnst.POLE_CONE_LEN*1e3)^2+(cnst.POLE_R*1e3)^2);
        axp=[-cos(inc); sin(inc)];                                   % P2 pole axis 在 x-z
        rot=@(v,a)[cos(a)*v(1)-sin(a)*v(2); sin(a)*v(1)+cos(a)*v(2)];
        Ct=T(:)+Lsl*rot(axp,+beta); Cb=T(:)+Lsl*rot(axp,-beta);
        plot([T(1) Ct(1)],[T(2) Ct(2)],'k-','LineWidth',2);
        plot([T(1) Cb(1)],[T(2) Cb(2)],'k-','LineWidth',2);
        % P1 cone — 下極「半切」：整錐軸在 z≈-0.289(WP)，BLOCK 切掉軸以上半(朝 WP 側)、留下半
        %   → 輪廓 = 平切上緣(水平、在軸高 z=-0.289、面朝 WP) + 下斜邊(-beta)，非對稱整錐
        T1  = [cnst.pole_tip_x(1), cnst.pole_tip_z_wp(1)]*1e3;        % P1 tip (x,z) WP frame [mm] ≈ (0.408,-0.289)
        top1= T1(:)+Lsl*[1;0];                                       % 半切平面（水平 +x，在軸高）
        low1= T1(:)+Lsl*rot([1;0],-beta);                           % 下斜邊（-beta，往右下）
        plot([T1(1) top1(1)],[T1(2) top1(2)],'k-','LineWidth',2);
        plot([T1(1) low1(1)],[T1(2) low1(2)],'k-','LineWidth',2);
    end

    %% ---- WP 中心標記（pole 模式：點明磁路起點）----
    if strcmp(MODE,'pole')
        plot(0,0,'p','MarkerSize',17,'MarkerFaceColor',[1 .85 .1],'MarkerEdgeColor','k','LineWidth',1);
        text(0.05*HW_MM,-0.05*HW_MM,'WP','FontSize',14,'FontWeight','bold','Interpreter','tex');
    end

    %% ---- 取樣圓柱外框（nodes 模式）----
    if strcmp(MODE,'nodes')
        adir=[nx;nz]; tdir=[-nz;nx];                                % 軸向 / 切向（x-z）
        L=AXIAL_TOL*1e3; R=SENSOR_R*1e3;                            % [mm]
        base=[sx;sz]; topc=base+L*adir;
        c1=base-R*tdir; c2=base+R*tdir; c3=topc+R*tdir; c4=topc-R*tdir;
        plot([c1(1) c2(1) c3(1) c4(1) c1(1)],[c1(2) c2(2) c3(2) c4(2) c1(2)], ...
             '--','Color',[.55 .1 .55],'LineWidth',2);              % 取樣圓柱框（紫虛線）
    end

    %% ---- sensor + n+ + disc ----
    plot(sx,sz,'o','MarkerSize',8,'MarkerFaceColor',[.1 .7 .25],'MarkerEdgeColor','k');
    nt=[sx;sz]+NLEN*[nx;nz];                                        % n+ 箭頭
    plot([sx nt(1)],[sz nt(2)],'-','Color','w','LineWidth',5.5);             % 白邊襯底
    plot([sx nt(1)],[sz nt(2)],'-','Color',[.92 .15 .15],'LineWidth',3.6);
    plot(nt(1),nt(2),'^','MarkerSize',13,'MarkerFaceColor',[.92 .15 .15],'MarkerEdgeColor','w','LineWidth',1.0);
    text(nt(1)+0.04*HW_MM,nt(2),'n_+','Color',[.92 .15 .15],'FontSize',15,'FontWeight','bold','Interpreter','tex');
    % disc 平面（⊥ n+）短線
    tdir=[-nz; nx]; dd=[sx;sz]-0.15*tdir; de=[sx;sz]+0.15*tdir;
    plot([dd(1) de(1)],[dd(2) de(2)],'-','Color',[.1 .5 .15],'LineWidth',3);

    %% ---- why 模式：標出「匯進→漏出」轉換 + 兩區（解釋為何 sensor 那塊朝外）----
    if strcmp(MODE,'why')
        plot(Pflip(1),Pflip(2),'d','MarkerSize',12,'MarkerFaceColor','w', ...
             'MarkerEdgeColor','k','LineWidth',1.8);                  % sign-flip 轉換點
        text(Pflip(1)+0.06*HW_MM,Pflip(2),' 匯進↔漏出 轉換 (沿面~3mm)', ...
             'FontSize',11,'FontWeight','bold','Interpreter','tex');
        Pin = Twh(:) + 1.7*flank_dir;                                % 近尖端側（匯進區）
        text(Pin(1),Pin(2),'磁通匯進 P2','FontSize',13,'FontWeight','bold', ...
             'Color',[.1 .45 .1],'HorizontalAlignment','center','Interpreter','tex');
        text(sx+0.18*HW_MM,sz+0.18*HW_MM,'sensor 在漏出段→朝外', ...
             'FontSize',13,'FontWeight','bold','Color',[.7 .1 .1],'Interpreter','tex');
    end

    hold(ax,'off'); axis(ax,'equal');
    % nodes 取節點雲質心（縮留白）；其餘 vx,vz 已在上面設好
    if strcmp(MODE,'nodes'), vx=mean(Xx); vz=mean(Zz); end
    if ~isempty(WIN), xlim(WIN(1:2)); ylim(WIN(3:4)); else, xlim([vx-HW_MM vx+HW_MM]); ylim([vz-HW_MM vz+HW_MM]); end
    xlabel('x (mm) (WP frame)'); ylabel('z (mm) (WP frame)');
    if strcmp(VARIANT,'standard'), srcTag='baseline';
    elseif ~isempty(regexp(VARIANT,'^gap\d+um_mueq$','once'))
        gtok=regexp(VARIANT,'gap(\d+)um','tokens'); srcTag=sprintf('gap %sµm', gtok{1}{1});
    else, srcTag='sensor_spheres'; end
    if strcmp(FLIP,'flip'), signTag='all-source'; else, signTag='raw 原始'; end
    if strcmp(MODE,'nodes')
        title(sprintf('P2 sensor 取樣圓柱內 %d 真實節點 ｜ P1 激發 (%s) ｜ %s', nSel, signTag, srcTag), 'Interpreter','none');
    elseif strcmp(MODE,'why')
        title({sprintf('為何 P2 sensor 那塊磁路朝外 ｜ P1 激發 (%s) ｜ %s', signTag, srcTag), ...
               '沿 P2 cone：近尖端匯進 P2，過了~3mm 轉成漏出/回 yoke，sensor(4.57mm) 落在漏出段'}, 'Interpreter','none');
    else
        title(sprintf('P2 sensor ｜ P1 激發 (%s) ｜ %s', signTag, srcTag), 'Interpreter','none');
    end
    ax.Toolbar.Visible = 'off';                                    % 匯出不帶 axes toolbar

    %% ---- [STYLE] 選項①：粗體框圖（main/rules/figure-style.md）----
    set(ax,'FontSize',16,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]);
    box(ax,'on'); grid(ax,'off');                                  % 外框框出 + 無背景網格
    xt=get(ax,'XTick'); set(ax,'XTick',xt(1:2:end));               % x tick 減半
    yt=get(ax,'YTick'); set(ax,'YTick',yt(1:2:end));               % y tick 減半
    set([get(ax,'XLabel') get(ax,'YLabel')],'FontSize',17,'FontWeight','bold');
    set(get(ax,'Title'),'FontSize',15,'FontWeight','bold');
    cb.FontSize=15; cb.FontWeight='bold'; cb.LineWidth=1.5;        % colorbar 字加大加粗
    set(get(cb,'YLabel'),'FontSize',16,'FontWeight','bold');
    ct=get(cb,'Ticks');                                            % colorbar tick 減半（保留頂端真實峰值）
    if numel(ct)>=4
        cl=get(cb,'TickLabels');
        keep=1:2:numel(ct); if keep(end)~=numel(ct), keep=[keep numel(ct)]; end
        set(cb,'Ticks',ct(keep));
        if iscell(cl) && numel(cl)==numel(ct), set(cb,'TickLabels',cl(keep)); end
    end

    %% ---- 輸出 ----
    if strcmp(VARIANT,'standard'), vsuf='_standard';
    elseif ~isempty(regexp(VARIANT,'^gap\d+um_mueq$','once')), vsuf=['_' VARIANT];
    else, vsuf=''; end   % baseline / gap 加後綴，不蓋加密版
    if strcmp(FLIP,'raw'), fsuf='_raw'; else, fsuf=''; end                % raw 版加後綴，不蓋 all-source 版
    if ~isempty(WIN), wsuf=sprintf('_x%gto%g_z%gto%g',WIN(1),WIN(2),WIN(3),WIN(4)); else, wsuf=''; end  % 自訂視窗加後綴
    if FILLEMPTY, isuf='_fill'; else, isuf=''; end                        % 內插填空版加後綴，不蓋 raw 版
    if PREVIEW
        out = fullfile(tempdir,['p2sensor_preview_' MODE vsuf fsuf wsuf isuf '.png']);
    else
        out = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
               'long2016_hexapole_halfcut\Calibration_using_FEM_modeling\Hall_sensor_base_fix_dir\figures\' ...
               'P2sensor_Braw_P1exc_' MODE vsuf fsuf wsuf isuf '.png'];
    end
    resOut = DPI; if PREVIEW, resOut = 150; end                    % preview 用 150 DPI（<2000px 可目視）
    exportgraphics(fig,out,'Resolution',resOut);
    fprintf('saved: %s\n', out);
end
