function plot_sensor_mounting_p2(SOFF, WITHFIELD, CMODE)
%   [ADDED] SOFF = 藍線長度（tip→foot 沿錐面直線距, mm）；省略 = 4.572（定案值）。
%   非預設值時輸出檔名自動加 _soff<值>mm 後綴，不覆蓋原圖。
%   [ADDED 2026-08-05] WITHFIELD=true：在**同一視野範圍**疊上 **P1 激發**（graded coil1、1A）
%     的磁路箭頭側視 → 另存 *_fieldP1.png（不覆蓋原圖）。
%     取樣 = |y| <= YSLAB mm 內的**真實 FEM 節點**，每格挑最接近 y=0 的那顆（**不內插**，
%     per plot-real-nodes）。號誌用激發極(P1，下極)的 all-source sign。
%     ⚠ colorbar 為**本圖自身範圍**（使用者 2026-08-05 指定）——P1 激發時此視窗內 ‖b‖ 僅 ~0.06 mT，
%       與其他磁路圖的 0–150 mT 色階**不可並列比較**。
% plot_sensor_mounting_p2 -- Section3_A paper figure: P2 Hall-sensor mounting side
% view (x-z), tip40µm geometry (幾何取自 mt_constants，非硬寫)。Clean paper style:
%   每軸 3 根等距 tick、tick 朝外、無尺寸數字、無軸標題、無 "P2 pole"/"sensor" 字、
%   感測器畫成實際圓柱(R0.15×H0.10mm)邊視、n+ 紅箭頭(小頭、無字)、字體 36。
%   藍線 tip→foot(沿錐面 4.572mm) + 橘線 foot→sensor(0.41mm，直接相連；tip40 不含 gap 段)。
%   輸出 → figures/paper_fig/Section3_A/sensor_mounting_tip40_P2.png
    here   = fileparts(fileparts(mfilename('fullpath')));
    figdir = fullfile(fileparts(here),'paper_fig','Section3_A');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    % [MODIFIED 2026-08-08] 不再 addpath backup（規則 no-backup-data）；改用 live config +
    %   utils/pole_sensor_geometry（sensor 幾何唯一來源）。幾何改用 CAD STEP 實測的真實錐體。
    CAL = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling';
    addpath(fullfile(CAL,'function'), fullfile(CAL,'utils'));
    cnst = model_config('long2016_hexapole_halfcut', 'tip40um');

    % ---- tip40 geometry (P2, WP frame, mm) ----
    IP = 2;                                                 % P2（幾何全部由 pole_sensor_geometry 供給）
    if nargin < 1 || isempty(SOFF), SOFF = 4.572; end       % [MODIFIED] 沿錐面直線距(藍線長)可調
    if nargin < 2 || isempty(WITHFIELD), WITHFIELD = false; end   % [ADDED]
    % [ADDED 2026-08-05] CMODE = 色階模式（視窗內 ‖b‖ 跨 4 個數量級，單一線性色階顯示不出結構）
    %   'max'（原始：0~視窗最大值）| 'pct'（0~98 分位，超出者 clamp）| 'log'（log10 色階）
    if nargin < 3 || isempty(CMODE), CMODE = 'max'; end   % [MODIFIED] 使用者偏好線性全範圍
    YSLAB = 4.0;   CELL = 0.07;                             % [ADDED] apdl 取樣：y 半寬 / 抽稀格邊 [mm]
    AIR = 0.41;                                             % 離面 [mm]
    rot = @(v,a)[cos(a)*v(1)-sin(a)*v(2); sin(a)*v(1)+cos(a)*v(2)];

    % [MODIFIED 2026-08-08] sensor 位置/法線一律由 pole_sensor_geometry 供給（不再自己算）。
    %   P2 的子午面 = y=0 → 3D 向量取 (x,z) 即得本圖的 2D 座標。
    [spos, snor, geo] = pole_sensor_geometry(cnst, struct('soff_upper', SOFF*1e-3));
    to2 = @(v) v([1 3])*1e3;                                % 3D[m] → 2D(x,z)[mm]
    beta = geo.beta(IP);                                    % 上極**真實**半錐角 11.0138°（CAD）
    rf   = geo.r_tip(IP)*1e3;                               % 尖端倒圓半徑 [mm]
    T    = to2(geo.tip(:,IP));
    axc  = geo.axis([1 3],IP);
    nh   = snor([1 3],IP);                                  % n+（出鋼）
    sensor = to2(spos(:,IP));                               % 圓柱底面中心（離錐面精確 0.41mm）
    foot   = sensor - AIR*nh;                               % 錐面上落腳點

    % ---- 磁極截面輪廓（tip + 錐軸 + ±β 錐面 + 40µm 鈍尖倒圓；用真實錐體）----
    % 弧自極尖掃到切點（半張角 90°−β），再沿母線延伸到錐底：
    %   Lsl = (真實錐長 − 倒圓軸向推進)/cos β
    Lsl = (geo.cone_len(IP) - geo.t_tan(IP))*1e3/cos(beta);
    C   = T + rf*axc;  axa = atan2(axc(2),axc(1));  angf = axa + pi;  hw = pi/2 - beta;
    fdA = rot(axc,+beta);  fdB = rot(axc,-beta);
    tha = linspace(angf-hw,angf+hw,60);  arc = C + rf*[cos(tha); sin(tha)];
    Cu  = arc(:,1) + Lsl*fdA;  Cd = arc(:,end) + Lsl*fdB;
    poly = [Cu, arc, Cd];

    %% ---- 視野範圍（先算，供場取樣與軸範圍共用；與原圖逐字相同的算法）----
    % [MODIFIED 2026-08-05] 原本在檔尾才算 xlim/ylim；WITHFIELD 需要先知道視窗才能取樣。
    NLEN = 1.0;  a1 = sensor + NLEN*nh;                     % n+ 箭頭尖端（原在下方，移上來）
    allx = [0 T(1) sensor(1) foot(1) a1(1)];
    allz = [0 T(2) sensor(2) a1(2)];
    pad  = 0.7;  XL = [min(allx)-pad max(allx)+pad];  ZL = [min(allz)-pad max(allz)+pad];

    %% ---- 繪圖 ----
    RED = [.92 .15 .15];  GRN = [.10 .70 .25];
    if WITHFIELD, fpos=[60 30 1240 980]; else, fpos=[60 30 980 980]; end   % [MODIFIED] 留 colorbar 空間
    figure('Position',fpos,'Color','w'); ax = axes; hold on;
    % [MODIFIED 2026-08-05] WITHFIELD 時輪廓改用 p1_maxwell_side / circuit_side 樣式
    %   （淡填 FaceAlpha 0.30 + 深灰藍邊），讓底下的磁路箭頭看得見；無場版維持原本實心黑邊。
    if WITHFIELD
        patch('XData',poly(1,:),'YData',poly(2,:),'FaceColor',[.82 .84 .88], ...
              'EdgeColor',[.28 .30 .36],'LineWidth',2.2,'FaceAlpha',0.30);
    else
        patch('XData',poly(1,:),'YData',poly(2,:),'FaceColor',[.80 .80 .84], ...
              'EdgeColor','k','LineWidth',2.5,'FaceAlpha',1);
    end

    %% ---- [ADDED 2026-08-05] P1 激發的磁路箭頭（真實 FEM 節點、不內插）----
    if WITHFIELD
        cf = fullfile(here, 'data', 'mount_field_P2fromP1_apdl.mat');
        if exist(cf,'file')
            L = load(cf,'Xs','Zs','Bx','Bz','Bm','XLc','ZLc');
            if isequal(L.XLc,XL) && isequal(L.ZLc,ZL)
                Xs=L.Xs; Zs=L.Zs; Bx=L.Bx; Bz=L.Bz; Bm=L.Bm;
                fprintf('磁路取樣: 由快取載入 %s（%d 支箭頭）\n', cf, numel(Xs));
            else
                fprintf('快取視窗不符 → 重算\n');
            end
        end
        if ~exist('Xs','var')
            CAL = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling';
            addpath(fullfile(CAL,'function'));  addpath(fullfile(CAL,'common_path'));
            zoff = -cnst.SPH_OFST*1e3;                       % raw z → WP frame
            d   = import_ansys_data(ansys_path('long2016_hexapole_halfcut','data','graded','coil1'), 'all', 'coil1');
            sgn = 1 - 2*cnst.pole_is_lower(1);               % 激發極 = P1（下極）→ all-source 取 -1
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
            sel = sel(d.bsum(sel) > 0);
            Xs = fx(sel);  Zs = fz(sel);
            Bx = sgn*d.bx(sel);  Bz = sgn*d.bz(sel);  Bm = d.bsum(sel)*1e3;   % mT
            XLc = XL;  ZLc = ZL;  save(cf,'Xs','Zs','Bx','Bz','Bm','XLc','ZLc');
            fprintf('磁路取樣[apdl,P1激發]: |y|<=%.1fmm、格 %.3fmm、每格挑最近 y=0 節點 → %d 支箭頭, |b| max %.4f mT（已存快取）\n', ...
                    YSLAB, CELL, numel(Xs), max(Bm));
        end
        % [MODIFIED 2026-08-05] 色階：視窗同時含 P1 尖端(~1070mT)與 P2 感測區(~0.06mT)，跨 4 個數量級
        nb = 28;  cmap = turbo(nb);  lw = linspace(0.5,2.2,nb);
        switch lower(CMODE)
            case 'max',  CLIM = max(Bm);              val = Bm;              cbl = '$\mathbf{\|b\|\;(mT)}$';
            case 'pct',  CLIM = prctile(Bm,98);       val = min(Bm,CLIM);    cbl = '$\mathbf{\|b\|\;(mT)}$';
            case 'log',  val = log10(max(Bm,1e-4));   CLIM = [];             cbl = '$\mathbf{\log_{10}\|b\|\;(mT)}$';
            otherwise,   error('CMODE 必為 max|pct|log');
        end
        if strcmpi(CMODE,'log'), lo = prctile(val,1); hi = max(val); else, lo = 0; hi = CLIM; end
        amax = 0.085;  bxz = hypot(Bx,Bz);  bxz(bxz==0) = 1e-12;
        u = (val-lo)/max(hi-lo,eps);  u = min(max(u,0),1);          % 0~1 正規化（決定色與長度）
        scl  = amax .* (0.25 + 0.75*u) ./ bxz;                       % 長度隨正規化值線性、不再被極值壓扁
        edg  = linspace(lo,hi,nb+1);
        for k = 1:nb
            if k < nb, m = val>=edg(k) & val<edg(k+1); else, m = val>=edg(k); end
            if any(m)
                quiver(ax, Xs(m), Zs(m), Bx(m).*scl(m), Bz(m).*scl(m), 0, ...
                       'Color',cmap(k,:), 'LineWidth',lw(k), 'MaxHeadSize',0.35);
            end
        end
        colormap(ax,turbo);  clim(ax,[lo hi]);  CBLAB = cbl;
    end

    % WP（十字，無字）+ WP→tip 連線
    plot([0 T(1)],[0 T(2)],'-','Color',[.5 .5 .5],'LineWidth',1.6);
    plot(0,0,'+','MarkerSize',22,'Color','k','LineWidth',3);

    % 藍線 tip→foot（沿錐面 4.572mm）+ 橘線 foot→sensor（0.41mm，直接相連、無數字）
    plot([T(1) foot(1)],[T(2) foot(2)],'-','Color',[.15 .35 .75],'LineWidth',3);
    plot([foot(1) sensor(1)],[foot(2) sensor(2)],'-','Color',[.75 .35 .1],'LineWidth',3);

    % tip（黑方塊，無字）
    plot(T(1),T(2),'s','MarkerSize',11,'MarkerFaceColor','k','MarkerEdgeColor','k');

    % sensor：實際圓柱邊視（2R=0.30 × H=0.10，軸沿 n+），綠填黑邊（無字）
    Rs = 0.15;  Hs = 0.10;  tdir = [-nh(2); nh(1)];
    base = sensor;  top = sensor + Hs*nh;
    rect = [base+Rs*tdir, base-Rs*tdir, top-Rs*tdir, top+Rs*tdir];
    patch('XData',rect(1,:),'YData',rect(2,:),'FaceColor',GRN,'EdgeColor','k','LineWidth',1.4);

    % n+：紅箭頭（桿 + 對齊 nh 的小三角頭，無字）   [MODIFIED] NLEN/a1 已於上方視野段先算
    hl = 0.20;  hw2 = 0.10;  bk = a1 - hl*nh;
    plot([sensor(1) bk(1)],[sensor(2) bk(2)],'-','Color',RED,'LineWidth',4.2);
    patch([a1(1) bk(1)+hw2*tdir(1) bk(1)-hw2*tdir(1)], ...
          [a1(2) bk(2)+hw2*tdir(2) bk(2)-hw2*tdir(2)], RED,'EdgeColor',RED);

    % ---- 軸：等比、box、每軸 3 根等距 tick、tick 朝外、字體 36、無軸標題 ----
    xlim(XL);  ylim(ZL);                                    % [MODIFIED] 用上方先算好的視窗（與場取樣一致）
    daspect([1 1 1]);  box on;  grid off;
    set(ax,'XTick',ticks3(xlim),'YTick',ticks3(ylim));
    set(ax,'FontSize',36,'FontWeight','bold','LineWidth',3.5,'TickLength',[.02 .02],'TickDir','out');
    ax.Toolbar.Visible = 'off';

    if WITHFIELD                                            % [ADDED] ‖b‖ colorbar（figure-style 標準記號）
        % [MODIFIED 2026-08-05] 先把座標軸往左縮，留出 colorbar + 標題的寬度，避免標題被裁掉。
        drawnow;  ap = get(ax,'Position');  set(ax,'Position',[ap(1) ap(2) ap(3)*0.74 ap(4)]);
        ap = get(ax,'Position');
        cb = colorbar(ax);  cb.Units = 'normalized';
        cb.Position = [ap(1)+ap(3)+0.025, ap(2), 0.030, ap(4)];
        cb.Label.Interpreter = 'latex';
        cb.Label.String = CBLAB;  cb.Label.FontSize = 36;
        cb.FontSize = 36;  cb.FontWeight = 'bold';
    end

    sfx = ''; if abs(SOFF-4.572) > 1e-9, sfx = sprintf('_soff%gmm', SOFF); end   % [ADDED] 非預設值另存
    if WITHFIELD, sfx = [sfx '_fieldP1_' lower(CMODE)]; end                                    % [ADDED] P1 激發磁路版另存
    out = fullfile(figdir, sprintf('sensor_mounting_tip40_P2%s.png', sfx));
    exportgraphics(gcf,out,'Resolution',200);
    fprintf('saved %s\n',out);
end

% ---- 每軸 3 根等距、nice-step tick（含在範圍內）----
function tk = ticks3(lm)
    s0 = (lm(2)-lm(1))/3;
    p  = 10^floor(log10(s0));  c = [1 2 3 4 5 10]*p;
    [~,i] = min(abs(c-s0));  s = c(i);
    ctr = round((lm(1)+lm(2))/2/s)*s;
    tk  = ctr + [-1 0 1]*s;
end
