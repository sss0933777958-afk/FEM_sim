function plot_sensorBcircuit_P1exc(pole_i, I_scale)
% PLOT_SENSORBCIRCUIT_P1EXC  任一極 sensor 局部磁路箭頭圖（P1 激發、all-source）。
%   pole_i  = sensor 極 (paper index 1..6)，預設 4 (P4)。
%   I_scale = 激發電流 (A)，預設 1。2/3 -> 載 coil1/standard_<I>A 線性放大資料
%             (|B| x I_scale；箭頭為單位方向不變，只 colorbar 平移)。
%   切面 = 該極「pole_axis × n+」局部平面（過 sensor），含真實 FEM 節點薄片（不內插）。
%   顯示方向 = 全域 +z 為上（world-up），故磁極依真實傾角呈現（如 P2 那張看得出上斜）。
%   資料：coil1 (P1 excited) FEM，套 all-source（P1 下極→整場 negate→source）。
%   箭頭=單位方向，顏色=|B|(3D, log)。輸出 figures/<Pn>sensor_Braw_P1exc_<I>A.png。

    if nargin < 1 || isempty(pole_i),  pole_i  = 4; end
    if nargin < 2 || isempty(I_scale), I_scale = 1; end

    %% ---- config ----
    SLAB_MM = 0.30;  HW_MM = 4.0;  ARROW_MM = 0.18;  DPI = 200;

    %% ---- paths ----
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
    addpath(['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
             'long2016_hexapole_halfcut\Calibration_using_FEM_modeling\Hall_sensor_base_no_fix_dir\code\function']);
    cnst = mt_constants();
    rr = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data';
    figdir = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
              'long2016_hexapole_halfcut\Calibration_using_FEM_modeling\Hall_sensor_base_no_fix_dir\figures'];

    %% ---- sensor 幾何 + 切面 + world-up 顯示基底 ----
    [sensor_pos, sensor_n] = build_sensor_geometry(cnst);
    c  = sensor_pos(:,pole_i);                         % 切面中心 = sensor [m]
    n  = sensor_n(:,pole_i);  n = n/norm(n);           % n+
    ax = cnst.pole_axis(:,pole_i);  ax = ax/norm(ax);  % 極軸 (tip->base)
    nrm = cross(ax, n);  nrm = nrm/norm(nrm);          % 切面法線（平面含 ax 與 n+）
    vup = [0;0;1] - ([0;0;1].'*nrm)*nrm; vup = vup/norm(vup);   % 平面內、對齊全域 +z
    vh  = cross(nrm, vup); vh = vh/norm(vh);                    % 平面內水平（極在 -x 側，對齊參考圖）
    fprintf('%s：n+·vup=%.2f（n+ 在顯示框的傾斜）\n', cnst.pole_labels{pole_i}, n.'*vup);

    %% ---- coil1 (P1 excited) + all-source negate ----
    if I_scale == 1, dfold = 'standard'; else, dfold = sprintf('standard_%dA',I_scale); end
    d = import_ansys_data(fullfile(rr,'coil1',dfold),'all','coil1');
    B = -[d.bx, d.by, d.bz];                           % all-source（P1 下極）；資料已含 xI_scale
    P = [d.x, d.y, d.z - cnst.SPH_OFST];               % WP 框 [m]

    %% ---- 真實節點薄片（切面內）+ 投到 WP-centred 顯示框 ----
    dn  = (P - c.')*nrm;                               % 距切面（過 sensor c）距離
    H   = P*vh*1e3;  V = P*vup*1e3;                    % WP 為原點的顯示座標 [mm]
    sc  = [c.'*vh; c.'*vup]*1e3;                        % sensor 顯示座標 [mm]
    m   = abs(dn) < SLAB_MM*1e-3 & abs(H) < HW_MM & abs(V) < HW_MM;
    fprintf('切面薄片內真實節點數 = %d\n', nnz(m));
    Hh = H(m); Vv = V(m);
    Bh = B(m,:)*vh;  Bv = B(m,:)*vup;  Bmag = sqrt(sum(B(m,:).^2,2));
    Bip = hypot(Bh,Bv); Bip(Bip<eps)=eps;
    aH = Bh./Bip*ARROW_MM;  aV = Bv./Bip*ARROW_MM;
    % 格點抽樣降密度（每格留「|B| 最大」的真實節點，峰值不被抽掉、colorbar 達真實峰值；仍是節點原值、非內插）
    CELL = 0.16;
    [~,~,gid] = unique(round([Hh Vv]/CELL),'rows');
    iu = accumarray(gid,(1:numel(Bmag))',[],@(r) r(find(Bmag(r)==max(Bmag(r)),1)));
    Hh=Hh(iu); Vv=Vv(iu); aH=aH(iu); aV=aV(iu); Bmag=Bmag(iu);

    %% ---- 圖 ----
    fig = figure('Position',[60 60 1000 950],'Color','w'); hold on;
    % LOG colour mapping (match reference): floor at cmin, log-spaced turbo bins
    cmax = max(Bmag);  cmin = max(min(Bmag(Bmag>0)), cmax*1e-4);
    nb=24; ed=logspace(log10(cmin),log10(cmax),nb+1); cm=turbo(nb);
    for k=1:nb
        if k==1,  in = Bmag<ed(k+1);
        elseif k==nb, in = Bmag>=ed(k);
        else, in = Bmag>=ed(k) & Bmag<ed(k+1); end
        if any(in), quiver(Hh(in),Vv(in),aH(in),aV(in),0,'Color',cm(k,:),'LineWidth',1.3,'MaxHeadSize',2.0); end
    end
    colormap(turbo); clim([cmin cmax]); set(gca,'ColorScale','log');
    cb=colorbar; ylabel(cb,'|B| (3D) [Tesla]');
    % [ADDED] integer-decade ticks + explicit real-max tick at colorbar top
    decs=ceil(log10(cmin)):floor(log10(cmax));
    if ~isempty(decs) && (log10(cmax)-decs(end))<0.15, decs(end)=[]; end  % drop top decade if it crowds the max
    tk=[10.^decs, cmax]; labs=[arrayfun(@(d)sprintf('10^{%d}',d),decs,'UniformOutput',false), {sprintf('%.3g',cmax)}];
    cb.Ticks=tk; cb.TickLabels=labs;

    pj = @(p3)[p3.'*vh; p3.'*vup]*1e3;                 % 3D 點(WP框) → 顯示 (H,V) [mm]
    d2 = @(v3)[v3.'*vh; v3.'*vup];                      % 3D 方向 → 顯示 (h,v)

    % ---- WP working-point star at display origin ----
    plot(0,0,'p','MarkerSize',18,'MarkerFaceColor',[1 .85 .1],'MarkerEdgeColor','k','LineWidth',1.0);
    text(0+0.12,0-0.18,'WP','FontSize',13,'FontWeight','bold');

    % sensor (at sc) + n+ + disc
    plot(sc(1),sc(2),'o','MarkerSize',8,'MarkerFaceColor',[.1 .7 .25],'MarkerEdgeColor','k');
    nd = d2(n); nd = nd/norm(nd)*0.85;
    plot(sc(1)+[0 nd(1)],sc(2)+[0 nd(2)],'-','Color','w','LineWidth',5.5);               % 白邊襯底
    plot(sc(1)+[0 nd(1)],sc(2)+[0 nd(2)],'-','Color',[.92 .15 .15],'LineWidth',3.6);
    plot(sc(1)+nd(1),sc(2)+nd(2),'^','MarkerSize',13,'MarkerFaceColor',[.92 .15 .15],'MarkerEdgeColor','w','LineWidth',1.0);
    text(sc(1)+nd(1)+0.08,sc(2)+nd(2),'n_+','Color',[.92 .15 .15],'FontSize',15,'FontWeight','bold','Interpreter','tex');
    e1 = ax - (ax.'*n)*n; e1=e1/norm(e1);              % 平面內 ⊥ n+（disc 切向）
    td = d2(e1); td=td/norm(td)*0.15;
    plot(sc(1)+[-td(1) td(1)],sc(2)+[-td(2) td(2)],'-','Color',[.1 .5 .15],'LineWidth',3);

    % ---- 磁極 cone 真實輪廓（3D 算、投到顯示框）----
    beta = atan2(3.0,15.0);
    Lsl  = sqrt(cnst.POLE_CONE_LEN^2 + cnst.POLE_R^2);   % slant 長 [m]
    tip3 = [cnst.pole_tip_x(pole_i); cnst.pole_tip_y(pole_i); cnst.pole_tip_z_wp(pole_i)];
    kxa  = cross(nrm,ax);                                 % ⊥ ax，平面內
    dT   = ax*cos(beta) + kxa*sin(beta);                  % 上斜邊方向 (Rodrigues, ax⊥nrm)
    dB   = ax*cos(beta) - kxa*sin(beta);                  % 下斜邊方向
    Tp = pj(tip3); Ct = pj(tip3 + Lsl*dT); Cb = pj(tip3 + Lsl*dB);
    plot([Tp(1) Ct(1)],[Tp(2) Ct(2)],'k-','LineWidth',1.8);
    plot([Tp(1) Cb(1)],[Tp(2) Cb(2)],'k-','LineWidth',1.8);
    plot([Ct(1) Cb(1)],[Ct(2) Cb(2)],'k-','LineWidth',1.8);
    plot(Tp(1),Tp(2),'ks','MarkerSize',8,'MarkerFaceColor','k');
    if abs(Tp(1))<HW_MM*1.3 && abs(Tp(2))<HW_MM*1.3, text(Tp(1)+0.08,Tp(2),'tip','FontSize',11); end

    hold off; axis equal; grid on;
    % window: include WP (0,0) and sensor (sc), with margin (reference layout)
    padx=1.2; padz=1.0;
    xlim([min(0,sc(1))-padx, max(0,sc(1))+padx]);
    ylim([min(0,sc(2))-padz, max(0,sc(2))+padz]);
    xlabel('x [mm] (WP frame)'); ylabel('z [mm] (WP frame)');
    title(sprintf('%s sensor ｜ P1 %gA', cnst.pole_labels{pole_i}, I_scale), 'Interpreter','none');

    if ~exist(figdir,'dir'); mkdir(figdir); end
    out = fullfile(figdir, sprintf('%ssensor_Braw_P1exc_%gA.png', cnst.pole_labels{pole_i}, I_scale));
    exportgraphics(fig,out,'Resolution',DPI);
    fprintf('saved: %s\n', out);
end
