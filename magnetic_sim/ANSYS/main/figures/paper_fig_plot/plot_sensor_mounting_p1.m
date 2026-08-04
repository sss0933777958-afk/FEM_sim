function plot_sensor_mounting_p1(SOFF, FACE, WITHFIELD, SRC)
%   [ADDED] SOFF = 藍線長度（tip→foot 沿貼附面直線距, mm）；省略 = 4.572（定案值）。
%   非預設值時輸出檔名自動加 _soff<值>mm 後綴，不覆蓋原圖。
%   [ADDED 2026-08-04] FACE = sensor 貼哪個面：
%     'cone'（預設）= 下錐面，n+ 朝下出鋼（定案位置）
%     'flat'        = **平切上表面**，沿水平錐軸走 SOFF、n+ = +z 朝上出鋼
%                     → 另存 sensor_mounting_tip40_P1_flat.png，不覆蓋原圖
% plot_sensor_mounting_p1 -- Section3_A paper figure: P1 (下磁極,半切) Hall-sensor
% mounting side view (x-z), tip40µm geometry (幾何取自 mt_constants，非硬寫)。
% Clean paper style (同 plot_sensor_mounting_p2)：每軸 3 根等距 tick、tick 朝外、
%   無尺寸數字、無軸標題、無字標(tip/sensor/P1 pole)、感測器畫實際圓柱(R0.15×H0.10mm)邊視、
%   n+ 紅箭頭(小頭、無字)、字體 36。藍線 tip→foot + 橘線 foot→sensor(直接相連)。
%   輸出 → figures/paper_fig/Section3_A/sensor_mounting_tip40_P1.png
    here   = fileparts(mfilename('fullpath'));
    figdir = fullfile(fileparts(here),'paper_fig','Section3_A');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');   % mt_constants

    % ---- tip40 geometry (P1: 0°、下極、半切、錐軸水平) ----
    cnst   = mt_constants();
    beta   = atan2(cnst.POLE_R, cnst.POLE_CONE_LEN);        % 半錐角 β ≈ 11.31°
    Lsl    = hypot(cnst.POLE_CONE_LEN, cnst.POLE_R)*1e3;    % 錐面斜長 ≈ 15.30 mm
    rf     = cnst.POLE_TIP_R*1e3;                           % 尖端倒圓半徑 [mm]（40µm = 0.04）
    if nargin < 1 || isempty(SOFF), SOFF = 4.572; end       % [MODIFIED] 沿貼附面直線距(藍線長)可調
    if nargin < 2 || isempty(FACE), FACE = 'cone'; end      % [ADDED] 'cone' | 'flat'
    % [ADDED 2026-08-04] WITHFIELD=true：疊上 APDL graded coil1（P1 自激發、1A）的磁路箭頭側視。
    %   取樣 = |y| <= YSLAB mm 內的**真實 FEM 節點**，每格挑最接近 y=0 的那顆（不內插）。
    %   水平軸延伸到 XMAX；箭頭 turbo 依 |B|(mT) 分 bin。→ 另存 *_field.png
    if nargin < 3 || isempty(WITHFIELD), WITHFIELD = false; end
    % [ADDED 2026-08-04] SRC = 場資料來源：'apdl'（graded 真實 FEM 節點）| 'maxwell'（.fld 規則格點）
    if nargin < 4 || isempty(SRC), SRC = 'apdl'; end
    YSLAB = 4.0;    XMAX = 6.0;    CELL = 0.07;             % apdl：y 取樣半寬 / x 上界 / 抽稀格邊 [mm]
    NRND  = 4200;                                            % maxwell：y=0 平面隨機位置內插的取樣點數
    AIR = 0.41;                                             % 離面 [mm]
    dir = @(el,az)[cos(el)*cos(az); sin(el)];
    rot = @(v,a)[cos(a)*v(1)-sin(a)*v(2); sin(a)*v(1)+cos(a)*v(2)];

    az = 0;
    T   = [cnst.pole_tip_x(1); cnst.pole_tip_z_wp(1)]*1e3;
    axc = dir(0, az);                                       % 錐軸水平 [1;0]
    if strcmpi(FACE,'flat')                                 % [ADDED] 平切上表面：面內方向 = 水平錐軸
        e2 = axc;                                           % 沿平切面（= 錐軸）
        nh = [0; 1];                                        % n+（朝上出鋼；鋼在平切面下方）
    else
        e2 = dir(-beta, az);                                % 下錐面 slant
        nh = dir(-beta - pi/2, az);                         % n+（朝下出鋼）
    end

    foot   = T + SOFF*e2;
    sensor = foot + AIR*nh;

    % [ADDED 2026-08-04] WITHFIELD：兩顆 sensor 都畫（平切上表面 + 下錐面），各自法線
    SP = [T + SOFF*axc          + AIR*[0;1], ...                 % 平切上表面（n+ = +z）
          T + SOFF*dir(-beta,az) + AIR*dir(-beta-pi/2,az)];      % 下錐面（n+ 朝下出鋼）
    SN = [[0;1], dir(-beta-pi/2,az)];

    % ---- 磁極截面輪廓（半切下極：平頂 + 鈍尖倒圓 + 下錐面）----
    C   = T + rf*axc;  axa = atan2(axc(2),axc(1));  angf = axa + pi;  hw = pi/2 - beta;
    fdB = rot(axc,-beta);
    tha = linspace(angf, angf+hw, 40);  arc = C + rf*[cos(tha); sin(tha)];
    low = arc(:,end) + Lsl*fdB;   top = T + Lsl*axc;
    poly = [top, T, arc, arc(:,end), low];

    %% ---- 繪圖 ----
    % ---- [MOVED UP] 視窗先算好（磁路取樣與 xlim/ylim 共用）----
    NLEN = 1.0;  a1 = sensor + NLEN*nh;                     % n+ 箭頭尖端
    pad  = 0.7;
    allx = [0 T(1) sensor(1) foot(1) a1(1)];
    allz = [0 T(2) sensor(2) a1(2)];
    if WITHFIELD                                            % [ADDED] 兩顆 sensor + 兩支箭頭都要進框
        A2 = SP + NLEN*SN;
        allx = [allx SP(1,:) A2(1,:)];
        allz = [allz SP(2,:) A2(2,:)];
    end
    XL   = [min(allx)-pad, max(allx)+pad];
    if WITHFIELD, XL(2) = XMAX; end                         % [ADDED] 水平軸延伸到 6 mm
    if strcmpi(FACE,'flat')      % 箭頭改朝上 → 下方要靠磁極下錐緣撐開，否則框被壓扁
        allz(end+1) = T(2) - (XL(2)-T(1))*tan(beta);        % 右框邊處的錐體下緣
    end
    ZL = [min(allz)-pad, max(allz)+pad];

    RED = [.92 .15 .15];  GRN = [.10 .70 .25];
    if WITHFIELD, fpos = [40 60 1560 900]; else, fpos = [60 30 980 980]; end
    figure('Position',fpos,'Color','w'); ax = axes; hold on;
    if WITHFIELD    % [MODIFIED] 輪廓改用 circuit_side / p1_maxwell_side 樣式（淡填 + 深灰藍邊）
        patch('XData',poly(1,:),'YData',poly(2,:),'FaceColor',[.82 .84 .88], ...
              'EdgeColor',[.28 .30 .36],'LineWidth',2.2,'FaceAlpha',0.30);
    else
        patch('XData',poly(1,:),'YData',poly(2,:),'FaceColor',[.80 .80 .84], ...
              'EdgeColor','k','LineWidth',2.5,'FaceAlpha',1);
    end

    % ---- [ADDED] 磁路箭頭（畫在磁極 patch 之上、標註之下；鋼內外皆顯示）----
    if WITHFIELD
        % [MODIFIED 2026-08-04] maxwell 走「快取原始 y=0 平面格點 → 每次現場抽樣」路徑：
        %   純均勻隨機會有團塊/空洞（Poisson clumping）→ 改 **jittered grid**（每格一點、格內隨機），
        %   覆蓋均勻且位置不規則。原始平面資料很小（~4k 點），改抽樣不必重讀 2GB。
        if strcmpi(SRC,'maxwell')
            zoff = -cnst.SPH_OFST*1e3;
            cs = fullfile(here,'mount_field_P1_maxwell_src.mat');
            if exist(cs,'file')
                L = load(cs,'SP2','BS','XLc','ZLc');
                if isequal(L.XLc,XL) && isequal(L.ZLc,ZL), SP2 = L.SP2;  BS = L.BS; end
                clear L;
            end
            if ~exist('SP2','var')
                MW = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\Maxwell';
                addpath(fullfile(MW,'function'));
                mcfg = model_config('long2016_hexapole_halfcut','tip40um');
                d = import_maxwell_fld(fullfile(mcfg.fld_dir, mcfg.fld_files_voltage{1}));
                fx = d.x*1e3;  fy = d.y*1e3;  fz = d.z*1e3 + zoff;
                pad2 = 0.25;                                  % 多留一圈，避免邊界外插成 NaN
                in = abs(fy) < 0.05 & fx>=XL(1)-pad2 & fx<=XL(2)+pad2 ...
                                    & fz>=ZL(1)-pad2 & fz<=ZL(2)+pad2;
                SP2 = [fx(in), fz(in)];
                BS  = mcfg.s_source(1) * [d.bx(in), d.by(in), d.bz(in)];   % 已含全 source 號誌
                XLc = XL;  ZLc = ZL;  save(cs,'SP2','BS','XLc','ZLc');
                clear d fx fy fz in;
                fprintf('y=0 平面原始格點 %d 個 → 已存 %s\n', size(SP2,1), cs);
            end
            Ix = scatteredInterpolant(SP2, BS(:,1), 'linear','none');
            Iy = scatteredInterpolant(SP2, BS(:,2), 'linear','none');
            Iz = scatteredInterpolant(SP2, BS(:,3), 'linear','none');
            % jittered grid：每格一點、格內隨機（JIT = 抖動佔格邊比例）
            cs2 = sqrt(diff(XL)*diff(ZL)/NRND);
            nxg = max(1,round(diff(XL)/cs2));   nzg = max(1,round(diff(ZL)/cs2));
            hx  = diff(XL)/nxg;                 hz  = diff(ZL)/nzg;
            [gx,gz] = meshgrid(0:nxg-1, 0:nzg-1);
            rng(0);   JIT = 0.9;
            Xs = XL(1) + (gx(:)+0.5 + JIT*(rand(numel(gx),1)-0.5))*hx;
            Zs = ZL(1) + (gz(:)+0.5 + JIT*(rand(numel(gz),1)-0.5))*hz;
            Bv = [Ix(Xs,Zs), Iy(Xs,Zs), Iz(Xs,Zs)];
            ok2 = all(isfinite(Bv),2) & vecnorm(Bv,2,2) > 1e-7;
            Xs = Xs(ok2);  Zs = Zs(ok2);  Bv = Bv(ok2,:);
            Bx = Bv(:,1);  Bz = Bv(:,3);  Bm = vecnorm(Bv,2,2)*1e3;
            fprintf('磁路取樣[maxwell]: jittered grid %dx%d（格 %.3fmm、抖動 %.0f%%）→ %d 支箭頭, |B| max %.1f mT\n', ...
                    nxg, nzg, hx, 100*JIT, numel(Xs), max(Bm));
        end
        % 抽稀後的繪圖資料快取（.dat 一顆 656k 節點、讀很久；純改樣式不需重載）
        cf = fullfile(here, sprintf('mount_field_P1_%s_%s.mat', lower(FACE), lower(SRC)));
        if ~strcmpi(SRC,'maxwell') && exist(cf,'file')
            L = load(cf,'Xs','Zs','Bx','Bz','Bm','XLc','ZLc');
            if isequal(L.XLc,XL) && isequal(L.ZLc,ZL)
                Xs=L.Xs; Zs=L.Zs; Bx=L.Bx; Bz=L.Bz; Bm=L.Bm;
                fprintf('磁路取樣: 由快取載入 %s（%d 支箭頭）\n', cf, numel(Xs));
            else
                fprintf('快取視窗不符 → 重算\n');  clear L;
            end
        end
        if ~exist('Xs','var')
            zoff = -cnst.SPH_OFST*1e3;                      % raw z → WP frame
            begin_apdl = true;  %#ok<NASGU>   （maxwell 已於前段完成取樣，這裡只走 apdl）
            if true
                CAL = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling';
                addpath(fullfile(CAL,'function'));  addpath(fullfile(CAL,'common_path'));
                d = import_ansys_data(ansys_path('long2016_hexapole_halfcut','data','graded','coil1'), 'all', 'coil1');
                sgn = 1 - 2*cnst.pole_is_lower(1);          % 全 source（下極 -1）
                fx = d.x*1e3;  fy = d.y*1e3;  fz = d.z*1e3 + zoff;
                in = fx>=XL(1) & fx<=XL(2) & fz>=ZL(1) & fz<=ZL(2) & abs(fy)<=YSLAB;
                gi = find(in);
                nx = round(diff(XL)/CELL);  nz = round(diff(ZL)/CELL);
                ix = discretize(fx(gi), linspace(XL(1),XL(2),nx+1));
                iz = discretize(fz(gi), linspace(ZL(1),ZL(2),nz+1));
                ok = ~isnan(ix) & ~isnan(iz);  gi = gi(ok);  ix = ix(ok);  iz = iz(ok);
                cid = (ix-1)*nz + iz;  ay = abs(fy(gi));
                [uc,~,g] = unique(cid);  sel = zeros(numel(uc),1);
                for k = 1:numel(uc), q = find(g==k); [~,b] = min(ay(q)); sel(k) = gi(q(b)); end
                sel = sel(d.bsum(sel) > 1e-4);
                Xs = fx(sel);  Zs = fz(sel);
                Bx = sgn*d.bx(sel);  Bz = sgn*d.bz(sel);  Bm = d.bsum(sel)*1e3;   % mT
                XLc = XL;  ZLc = ZL;  save(cf,'Xs','Zs','Bx','Bz','Bm','XLc','ZLc');
                fprintf('磁路取樣[apdl]: |y|<=%.1fmm、格 %.3fmm、每格挑最近 y=0 節點 → %d 支箭頭, |B| max %.1f mT（已存快取）\n', ...
                        YSLAB, CELL, numel(Xs), max(Bm));
            end
        end
        amax = 0.085;  bxz = hypot(Bx,Bz);  bxz(bxz==0) = 1e-12;
        scl = amax .* (Bm./max(Bm)).^0.25 ./ bxz;
        CLIM = ceil(max(Bm)/50)*50;  nb = 28;
        edg = linspace(0,CLIM,nb+1);  cmap = turbo(nb);  lw = linspace(0.5,2.2,nb);
        for k = 1:nb
            if k < nb, m = Bm>=edg(k) & Bm<edg(k+1); else, m = Bm>=edg(k); end
            if any(m)
                quiver(ax, Xs(m), Zs(m), Bx(m).*scl(m), Bz(m).*scl(m), 0, ...
                       'Color',cmap(k,:), 'LineWidth',lw(k), 'MaxHeadSize',0.35);
            end
        end
        colormap(ax,turbo);  clim(ax,[0 CLIM]);
    end

    Rs = 0.15;  Hs = 0.10;  hl = 0.20;  hw2 = 0.10;    % sensor 圓柱邊視尺寸 / 箭頭頭尺寸
    if WITHFIELD
        % [MODIFIED 2026-08-04] 只留 WP 十字 + **兩顆 sensor（平切上表面、下錐面）**與各自紅色 n+ 箭頭。
        %   拿掉：WP→tip 灰線、藍線(4.572mm)、橘線(0.41mm)、tip 黑方塊。
        plot(0,0,'+','MarkerSize',22,'Color','k','LineWidth',3);
        for k = 1:size(SP,2)
            sk = SP(:,k);  nk = SN(:,k);  tk = [-nk(2); nk(1)];
            rect = [sk+Rs*tk, sk-Rs*tk, sk+Hs*nk-Rs*tk, sk+Hs*nk+Rs*tk];
            patch('XData',rect(1,:),'YData',rect(2,:),'FaceColor',GRN,'EdgeColor','k','LineWidth',1.4);
            ak = sk + NLEN*nk;  bkk = ak - hl*nk;
            plot([sk(1) bkk(1)],[sk(2) bkk(2)],'-','Color',RED,'LineWidth',4.2);
            patch([ak(1) bkk(1)+hw2*tk(1) bkk(1)-hw2*tk(1)], ...
                  [ak(2) bkk(2)+hw2*tk(2) bkk(2)-hw2*tk(2)], RED,'EdgeColor',RED);
        end
    else
    % WP（十字，無字）+ WP→tip 連線
    plot([0 T(1)],[0 T(2)],'-','Color',[.5 .5 .5],'LineWidth',1.6);
    plot(0,0,'+','MarkerSize',22,'Color','k','LineWidth',3);

    % 藍線 tip→foot（沿錐面 4.572mm）+ 橘線 foot→sensor（0.41mm，直接相連、無數字）
    plot([T(1) foot(1)],[T(2) foot(2)],'-','Color',[.15 .35 .75],'LineWidth',3);
    plot([foot(1) sensor(1)],[foot(2) sensor(2)],'-','Color',[.75 .35 .1],'LineWidth',3);

    % tip（黑方塊，無字）
    plot(T(1),T(2),'s','MarkerSize',11,'MarkerFaceColor','k','MarkerEdgeColor','k');

    % sensor：實際圓柱邊視（2R=0.30 × H=0.10，軸沿 n+），綠填黑邊（無字）
    tdir = [-nh(2); nh(1)];
    base = sensor;  topc = sensor + Hs*nh;
    rect = [base+Rs*tdir, base-Rs*tdir, topc-Rs*tdir, topc+Rs*tdir];
    patch('XData',rect(1,:),'YData',rect(2,:),'FaceColor',GRN,'EdgeColor','k','LineWidth',1.4);

    % n+：紅箭頭（桿 + 對齊 nh 的小三角頭，無字）
    bk = a1 - hl*nh;
    plot([sensor(1) bk(1)],[sensor(2) bk(2)],'-','Color',RED,'LineWidth',4.2);
    patch([a1(1) bk(1)+hw2*tdir(1) bk(1)-hw2*tdir(1)], ...
          [a1(2) bk(2)+hw2*tdir(2) bk(2)-hw2*tdir(2)], RED,'EdgeColor',RED);
    end

    % ---- 軸：等比、box、每軸 3 根等距 tick、tick 朝外、字體 36、無軸標題 ----
    %   只框 tip+sensor 區(不含磁極 15mm 遠端)→ 磁極平頂延伸出框、自動裁切(同參考圖)。
    xlim(XL);  ylim(ZL);
    daspect([1 1 1]);  box on;  grid off;
    set(ax,'XTick',ticks3(xlim),'YTick',ticks3(ylim));
    set(ax,'FontSize',36,'FontWeight','bold','LineWidth',3.5,'TickLength',[.02 .02],'TickDir','out');
    if WITHFIELD                                             % [ADDED] |b| colorbar（同 3D 圖記號）
        cb = colorbar(ax);  cb.Label.Interpreter = 'latex';
        cb.Label.String = '$\mathbf{\|b\|\;(mT)}$';  cb.Label.FontSize = 36;
        cb.FontSize = 36;  cb.FontWeight = 'bold';
    end
    ax.Toolbar.Visible = 'off';

    sfx = ''; if abs(SOFF-4.572) > 1e-9, sfx = sprintf('_soff%gmm', SOFF); end   % [ADDED] 非預設值另存
    if strcmpi(FACE,'flat'), sfx = ['_flat' sfx]; end                            % [ADDED] 平切面版另存
    if WITHFIELD,            sfx = [sfx '_field'];  end                          % [ADDED] 疊磁路版另存
    if WITHFIELD && strcmpi(SRC,'maxwell'), sfx = [sfx '_maxwell']; end          % [ADDED] 場來源後綴
    out = fullfile(figdir, sprintf('sensor_mounting_tip40_P1%s.png', sfx));
    exportgraphics(gcf,out,'Resolution',200);
    fprintf('saved %s\n',out);
end

% ---- 每軸 3 根等距、nice-step tick（含在範圍內）----
function tk = ticks3(lm)
    s0 = (lm(2)-lm(1))/3;
    p  = 10^floor(log10(s0));  c = [1 2 3 4 5 10]*p;
    [~,i] = min(abs(c-s0));
    % [ADDED 2026-08-04] 步長太大時 3 根 tick 會有一根落在框外被裁掉（只剩 2 根數字）→
    %   往下一格 nice 值降，直到三根都在範圍內。既有圖三根本來就都在框內，輸出不變。
    for k = i:-1:1
        s = c(k);  ctr = round((lm(1)+lm(2))/2/s)*s;  tk = ctr + [-1 0 1]*s;
        if tk(1) >= lm(1) && tk(3) <= lm(2), return; end
    end
    s = c(1);  ctr = round((lm(1)+lm(2))/2/s)*s;  tk = ctr + [-1 0 1]*s;   % 保險
end
